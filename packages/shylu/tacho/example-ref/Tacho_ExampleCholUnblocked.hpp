#ifndef __TACHO_EXAMPLE_CHOL_UNBLOCKED_HPP__
#define __TACHO_EXAMPLE_CHOL_UNBLOCKED_HPP__

#include <Kokkos_Core.hpp>
#include <impl/Kokkos_Timer.hpp>

#include "Tacho_Util.hpp"

#include "Tacho_CrsMatrixBase.hpp"
#include "Tacho_CrsMatrixView.hpp"
#include "Tacho_CrsRowView.hpp"
#include "Tacho_CrsMatrixTools.hpp"

#include "Tacho_MatrixMarket.hpp"

#include "Tacho_GraphTools.hpp"

#include "Tacho_GraphTools_Scotch.hpp"
#include "Tacho_GraphTools_CAMD.hpp"

#include "Tacho_SymbolicFactorization.hpp"

#include "Tacho_TaskView.hpp"
#include "Tacho_TaskFactory.hpp"

#include "Tacho_Gemm.hpp"
#include "Tacho_Herk.hpp"
#include "Tacho_Trsm.hpp"
#include "Tacho_Chol.hpp"

namespace Tacho {

  template<typename DeviceSpaceType>
  int exampleCholUnblocked(const std::string file_input,
                           const int treecut,
                           const int prunecut,
                           const int fill_level,
                           const int rows_per_team,
                           const bool verbose) {
    typedef typename
      Kokkos::Impl::is_space<DeviceSpaceType>::host_mirror_space::execution_space HostSpaceType ;

    const bool detail = false;
    std::cout << "DeviceSpace::  "; DeviceSpaceType::print_configuration(std::cout, detail);
    std::cout << "HostSpace::    ";   HostSpaceType::print_configuration(std::cout, detail);

    typedef CrsMatrixBase<value_type,ordinal_type,size_type,HostSpaceType> CrsMatrixBaseHostType;
    typedef GraphTools<ordinal_type,size_type,HostSpaceType> GraphToolsHostType;

    typedef GraphTools_Scotch<ordinal_type,size_type,HostSpaceType> GraphToolsHostType_Scotch;
    typedef GraphTools_CAMD<ordinal_type,size_type,HostSpaceType> GraphToolsHostType_CAMD;

    typedef SymbolicFactorization<CrsMatrixBaseHostType> SymbolicFactorizationType;

    typedef Kokkos::Experimental::TaskPolicy<DeviceSpaceType> PolicyType;

    typedef CrsMatrixBase<value_type,ordinal_type,size_type,DeviceSpaceType> CrsMatrixBaseDeviceType;

    typedef CrsMatrixView<CrsMatrixBaseDeviceType> CrsMatrixViewDeviceType;
    typedef TaskView<CrsMatrixViewDeviceType> CrsTaskViewDeviceType;

    int r_val = 0;

    Kokkos::Impl::Timer timer;

    CrsMatrixBaseHostType AA_host("AA_host");
    timer.reset();
    {
      std::ifstream in;
      in.open(file_input);
      if (!in.good()) {
        std::cout << "Failed in open the file: " << file_input << std::endl;
        return -1;
      }
      MatrixMarket::read(AA_host, in);
    }
    double t_read = timer.seconds();

    if (verbose)
      AA_host.showMe(std::cout) << std::endl;

    typename GraphToolsHostType::size_type_array rptr("Graph::RowPtrArray", AA_host.NumRows() + 1);
    typename GraphToolsHostType::ordinal_type_array cidx("Graph::ColIndexArray", AA_host.NumNonZeros());

    timer.reset();
    GraphToolsHostType::getGraph(rptr, cidx, AA_host);
    double t_graph = timer.seconds();

    GraphToolsHostType_Scotch S;
    S.setGraph(AA_host.NumRows(), rptr, cidx);
    S.setSeed(0);
    S.setTreeLevel();
    S.setStrategy( SCOTCH_STRATSPEED
                   | SCOTCH_STRATLEVELMAX  
                   | SCOTCH_STRATLEVELMIN  
                   | SCOTCH_STRATLEAFSIMPLE
                   | SCOTCH_STRATSEPASIMPLE
                   );

    timer.reset();
    S.computeOrdering(treecut);
    double t_scotch = timer.seconds();

    S.pruneTree(prunecut);
    if (verbose)
      S.showMe(std::cout) << std::endl;

    CrsMatrixBaseHostType BB_host("BB_host");
    BB_host.createConfTo(AA_host);

    CrsMatrixTools::copy(BB_host,
                         S.PermVector(),
                         S.InvPermVector(),
                         AA_host);

    if (verbose)
      BB_host.showMe(std::cout) << std::endl;

    timer.reset();
    GraphToolsHostType::getGraph(rptr, cidx, BB_host);
    t_graph += timer.seconds();

    GraphToolsHostType_CAMD C;
    C.setGraph(BB_host.NumRows(),
               rptr, cidx,
               S.NumBlocks(),
               S.RangeVector());

    timer.reset();
    C.computeOrdering();
    double t_camd = timer.seconds();

    if (verbose)
      C.showMe(std::cout) << std::endl;

    CrsMatrixBaseHostType CC_host("CC_host");
    CC_host.createConfTo(BB_host);

    CrsMatrixTools::copy(CC_host,
                         C.PermVector(),
                         C.InvPermVector(),
                         BB_host);

    if (verbose)
      CC_host.showMe(std::cout) << std::endl;

    CrsMatrixBaseHostType DD_host("DD_host");

    timer.reset();
    SymbolicFactorizationType::createNonZeroPattern(DD_host,
                                                    fill_level,
                                                    Uplo::Upper,
                                                    CC_host,
                                                    rows_per_team);
    double t_symbolic = timer.seconds();

    if (verbose)
      DD_host.showMe(std::cout) << std::endl;

    // ==================================================================================

    CrsMatrixBaseDeviceType AA_device("AA_device");
    AA_device.mirror(DD_host);

    const size_type max_concurrency = 10;
    const size_type max_task_size = (3*sizeof(CrsTaskViewDeviceType)+sizeof(PolicyType)+128);
    const size_type max_task_dependence = 0;
    const size_type team_size = 1;

    PolicyType policy(max_concurrency,
                      max_task_size,
                      max_task_dependence,
                      team_size);

    CrsMatrixViewDeviceType A_device(AA_device);
    Kokkos::View<typename CrsMatrixViewDeviceType::row_view_type*,DeviceSpaceType> 
      rowviews("RowViewInMatView", AA_device.NumRows());
    A_device.setRowViewArray(rowviews);
    
    timer.reset();
    int ierr = Chol<Uplo::Upper,AlgoChol::Unblocked,Variant::One>::invoke
      (policy, policy.member_single(),
       A_device);
    double t_chol = timer.seconds();
    TACHO_TEST_FOR_ABORT( ierr, "Fail to perform Cholesky (serial)");

    if (verbose) {
      DD_host.mirror(AA_device);
      DD_host.showMe(std::cout) << std::endl;
    }

    {
      const auto prec = std::cout.precision();
      std::cout.precision(4);

      std::cout << std::scientific;
      std::cout << "SymbolicFactorization:: Given matrix  dimension = " << AA_host.NumRows() << " x " << AA_host.NumCols()
                << ", " << " nnz = " << AA_host.NumNonZeros() << std::endl;
      std::cout << "SymbolicFactorization:: Upper factors dimension = " << DD_host.NumRows() << " x " << DD_host.NumCols()
                << ", " << " nnz = " << DD_host.NumNonZeros() << std::endl;

      std::cout << "SymbolicFactorization:: "
                << "read = " << t_read << " [sec], "
                << "graph generation = " << (t_graph/2.0) << " [sec] "
                << "scotch reordering = " << t_scotch << " [sec] "
                << "camd reordering = " << t_camd << " [sec] "
                << "symbolic factorization = " << t_symbolic << " [sec] "
                << "Cholesky factorization = " << t_chol << " [sec] "
                << std::endl;

      std::cout.unsetf(std::ios::scientific);
      std::cout.precision(prec);
    }

    return r_val;
  }
}

#endif
