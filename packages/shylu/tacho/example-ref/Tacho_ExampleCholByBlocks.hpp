#ifndef __TACHO_EXAMPLE_CHOL_BY_BLOCKS_HPP__
#define __TACHO_EXAMPLE_CHOL_BY_BLOCKS_HPP__

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
  int exampleCholByBlocks(const std::string file_input,
                          const int treecut,
                          const int prunecut,
                          const int fill_level,
                          const int rows_per_team,
                          const int max_concurrency,
                          const int max_task_dependence,
                          const int team_size,
                          const bool check,
                          const bool verbose) {
    typedef typename
      Kokkos::Impl::is_space<DeviceSpaceType>::host_mirror_space::execution_space HostSpaceType ;

    const bool detail = false;
    std::cout << "DeviceSpace::  "; DeviceSpaceType::print_configuration(std::cout, detail);
    std::cout << "HostSpace::    ";   HostSpaceType::print_configuration(std::cout, detail);

    // for simple test, let's use host space only here, for device it needs mirroring.

    typedef CrsMatrixBase<value_type,ordinal_type,size_type,HostSpaceType> CrsMatrixBaseHostType;
    typedef CrsMatrixView<CrsMatrixBaseHostType> CrsMatrixViewHostType;

    typedef GraphTools<ordinal_type,size_type,HostSpaceType> GraphToolsHostType;

    typedef GraphTools_Scotch<ordinal_type,size_type,HostSpaceType> GraphToolsHostType_Scotch;
    typedef GraphTools_CAMD<ordinal_type,size_type,HostSpaceType> GraphToolsHostType_CAMD;

    typedef SymbolicFactorization<CrsMatrixBaseHostType> SymbolicFactorizationType;

    typedef Kokkos::Experimental::TaskPolicy<DeviceSpaceType> PolicyType;

    typedef TaskView<CrsMatrixViewHostType> CrsTaskViewHostType;
    typedef CrsMatrixBase<CrsTaskViewHostType,ordinal_type,size_type,HostSpaceType> CrsHierBaseHostType;
    typedef CrsMatrixView<CrsHierBaseHostType> CrsHierViewHostType;
    typedef TaskView<CrsHierViewHostType> CrsTaskHierViewHostType;

    int r_val = 0;
    
    Kokkos::Impl::Timer timer;

    ///
    /// Read from matrix market
    ///
    ///     input  - file
    ///     output - AA_host
    ///
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

    ///
    /// Create a graph structure for Scotch and CAMD (rptr, cidx)
    ///     
    ///     rptr and cidx are need to be set up for Scotch and CAMD
    ///
    typename GraphToolsHostType::size_type_array rptr("Graph::RowPtrArray", AA_host.NumRows() + 1);
    typename GraphToolsHostType::ordinal_type_array cidx("Graph::ColIndexArray", AA_host.NumNonZeros());

    ///
    /// Run Scotch
    ///
    ///     input  - rptr, cidx, A_host
    ///     output - S (perm, iperm, nblks, range, tree), AA_scotch_host (permuted)
    ///
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
    
    if (verbose)
      S.showMe(std::cout) << std::endl;

    CrsMatrixBaseHostType AA_scotch_host("AA_scotch_host");
    AA_scotch_host.createConfTo(AA_host);
    
    CrsMatrixTools::copy(AA_scotch_host,
                         S.PermVector(),
                         S.InvPermVector(),
                         AA_host);

    if (verbose)
      AA_scotch_host.showMe(std::cout) << std::endl;

    ///
    /// Run CAMD
    ///
    ///     input  - rptr, cidx, A_host
    ///     output - S (perm, iperm, nblks, range, tree), AA_scotch_host (permuted)
    ///
    timer.reset();
    GraphToolsHostType::getGraph(rptr, cidx, AA_scotch_host);
    t_graph += timer.seconds();

    GraphToolsHostType_CAMD C;
    C.setGraph(AA_scotch_host.NumRows(),
               rptr, cidx,
               S.NumBlocks(),
               S.RangeVector());

    timer.reset();
    C.computeOrdering();
    double t_camd = timer.seconds();

    if (verbose)
      C.showMe(std::cout) << std::endl;
    
    CrsMatrixBaseHostType AA_camd_host("AA_camd_host");
    AA_camd_host.createConfTo(AA_scotch_host);
    
    CrsMatrixTools::copy(AA_camd_host,
                         C.PermVector(),
                         C.InvPermVector(),
                         AA_scotch_host);
    
    if (verbose)
      AA_camd_host.showMe(std::cout) << std::endl;

    ///
    /// Symbolic factorization
    ///
    ///     input  - 
    ///     output - S (perm, iperm, nblks, range, tree), AA_scotch_host (permuted)
    ///
    CrsMatrixBaseHostType AA_factor_host("AA_factor_host");

    timer.reset();
    SymbolicFactorizationType::createNonZeroPattern(AA_factor_host,
                                                    fill_level,
                                                    Uplo::Upper,
                                                    AA_camd_host,
                                                    rows_per_team);
    double t_symbolic = timer.seconds();

    if (verbose)
      AA_factor_host.showMe(std::cout) << std::endl;

    ///
    /// Clean tempoerary matrices
    ///
    ///     input  - AA_scotch_host, AA_camd_host, C, rptr, cidx
    ///     output - none
    ///
    AA_scotch_host = CrsMatrixBaseHostType();
    AA_camd_host   = CrsMatrixBaseHostType();

    C = GraphToolsHostType_CAMD();
    rptr = typename GraphToolsHostType::size_type_array();
    cidx = typename GraphToolsHostType::ordinal_type_array();

    ///
    /// Create task policy
    ///
    ///     input  - max_task_size
    ///     output - policy
    ///
    const size_type max_task_size = (3*sizeof(CrsTaskViewHostType)+sizeof(PolicyType)+128);

    timer.reset();
    PolicyType policy(max_concurrency,
                      max_task_size,
                      max_task_dependence,
                      team_size);
    double t_policy = timer.seconds();

    ///
    /// Sequential execution
    ///
    ///     input  - AA_factor_host (matrix to be compared), rowviews
    ///     output - BB_factor_host, B_factor_host
    ///
    double t_chol_serial = 0;
    CrsMatrixBaseHostType BB_factor_host("BB_factor_host");
    if (check) {
      BB_factor_host.createConfTo(AA_factor_host);
      CrsMatrixTools::copy(BB_factor_host, AA_factor_host);
      
      CrsTaskViewHostType B_factor_host(BB_factor_host);
      Kokkos::View<typename CrsTaskViewHostType::row_view_type*,HostSpaceType>
        rowviews("RowViewInMatView", B_factor_host.NumRows());
      B_factor_host.setRowViewArray(rowviews);
      
      timer.reset();
      {
        auto future = policy.proc_create_team(Chol<Uplo::Upper,AlgoChol::Unblocked,Variant::One>
                                              ::createTaskFunctor(policy, B_factor_host));
        policy.spawn(future);
        Kokkos::Experimental::wait(policy);
        TACHO_TEST_FOR_ABORT( future.get(), "Fail to perform Cholesky (serial)");
      }
      t_chol_serial = timer.seconds();

      if (verbose)
        BB_factor_host.showMe(std::cout) << std::endl;
    }

    ///
    /// Task parallel execution
    ///
    ///    input  - AA_factor_host, rowviews
    ///    output - HA_factor_host, AA_factor_host, B_factor_host
    ///
    double t_hier = 0, t_blocks = 0, t_chol_parallel = 0;
    CrsHierBaseHostType HA_factor_host("HA_factor_host");
    {
      timer.reset();
      S.pruneTree(prunecut);
      CrsMatrixTools::createHierMatrix(HA_factor_host, 
                                       AA_factor_host, 
                                       S.NumBlocks(),
                                       S.RangeVector(),
                                       S.TreeVector());
      t_hier = timer.seconds();    
      
      timer.reset();
      size_type nblocks = HA_factor_host.NumNonZeros();
      
      Kokkos::View<ordinal_type*,HostSpaceType> 
        ap_rowview_blocks("NumRowViewInBlocks", nblocks + 1);
      
      ap_rowview_blocks(0) = 0;
      for (ordinal_type k=0;k<nblocks;++k) 
        ap_rowview_blocks(k+1) = ap_rowview_blocks(k) + HA_factor_host.Value(k).NumRows();
      
      Kokkos::View<typename CrsMatrixViewHostType::row_view_type*,HostSpaceType>
        rowview_blocks("RowViewInBlocks", ap_rowview_blocks(nblocks));
      
      Kokkos::parallel_for(Kokkos::RangePolicy<HostSpaceType>(0, nblocks),
                           [&](const ordinal_type k) {
                             const ordinal_type begin = ap_rowview_blocks(k);
                             const ordinal_type end   = ap_rowview_blocks(k+1);
                             HA_factor_host.Value(k).setRowViewArray
                               (Kokkos::subview(rowview_blocks, Kokkos::pair<ordinal_type,ordinal_type>(begin, end)));
                           } );
      CrsMatrixTools::filterEmptyBlocks(HA_factor_host);
      t_blocks = timer.seconds();    

      {
        size_type nblocks_filtered = HA_factor_host.NumNonZeros(), nnz_blocks = 0;
        for (size_type k=0;k<nblocks_filtered; ++k) 
          nnz_blocks += HA_factor_host.Value(k).NumNonZeros();
        
        TACHO_TEST_FOR_ABORT( nnz_blocks != AA_factor_host.NumNonZeros(),
                              "nnz counted in blocks is different from nnz in the base matrix.");
      }
      
      CrsTaskHierViewHostType H_factor_host(HA_factor_host);
      timer.reset();
      {
        auto future = policy.proc_create_team(Chol<Uplo::Upper,AlgoChol::ByBlocks,Variant::One>
                                              ::createTaskFunctor(policy, H_factor_host));
        policy.spawn(future);
        Kokkos::Experimental::wait(policy);
        TACHO_TEST_FOR_ABORT( future.get(), "Fail to perform Cholesky (serial)");
      }
      t_chol_parallel = timer.seconds();
      
      if (verbose)
        AA_factor_host.showMe(std::cout) << std::endl;
    }

    if (check) {
      double diff = 0, norm = 0;    
      TACHO_TEST_FOR_ABORT( BB_factor_host.NumNonZeros() != AA_factor_host.NumNonZeros(),
                            "nnz used in serial is not same as nnz used in parallel");

      const size_type nnz = AA_factor_host.NumNonZeros();
      for (size_type k=0;k<nnz;++k) {
        norm += Util::abs(BB_factor_host.Value(k));
        diff += Util::abs(AA_factor_host.Value(k) - BB_factor_host.Value(k));
      }
      std::cout << std::scientific;
      std::cout << "CholByBlocks:: check with serial execution " << std::endl
                << "               diff = " << diff << ", norm = " << norm << ", rel err = " << (diff/norm) << std::endl;
      std::cout.unsetf(std::ios::scientific);      
    }

    {
      const auto prec = std::cout.precision();
      std::cout.precision(4);

      std::cout << std::scientific;
      std::cout << "CholByBlocks:: Given    matrix = " << AA_host.NumRows() << " x " << AA_host.NumCols() << ", nnz = " << AA_host.NumNonZeros() << std::endl;
      std::cout << "CholByBlocks:: Factored matrix = " << AA_factor_host.NumRows() << " x " << AA_factor_host.NumCols() << ", nnz = " << AA_factor_host.NumNonZeros() << std::endl;
      std::cout << "CholByBlocks:: Hier     matrix = " << HA_factor_host.NumRows() << " x " << HA_factor_host.NumCols() << ", nnz = " << HA_factor_host.NumNonZeros() << std::endl;

      std::cout << "CholByBlocks:: "
                << "read = " << t_read << " [sec], "
                << "graph generation = " << (t_graph/2.0) << " [sec] "
                << "scotch reordering = " << t_scotch << " [sec] "
                << "camd reordering = " << t_camd << " [sec] " 
                << std::endl
                << "CholByBlocks:: "
                << "symbolic factorization = " << t_symbolic << " [sec] "
                << std::endl
                << "CholByBlocks:: "
                << "policy creation = " << t_policy << " [sec] "
                << "hier creation = " << t_hier << " [sec] "
                << "block specification = " << t_blocks << " [sec] "
                << std::endl
                << "CholByBlocks:: "
                << "Chol Parallel = " << t_chol_parallel << " [sec] ";
      if (check) 
        std::cout << "Chol Serial = " << (check ? t_chol_serial : -1) << " [sec] "
                  << "speed-up = " << (t_chol_serial/t_chol_parallel) << " [sec] ";
      
      std::cout << std::endl;
      
      std::cout.unsetf(std::ios::scientific);
      std::cout.precision(prec);
    }

    return r_val;
  }
}

#endif
