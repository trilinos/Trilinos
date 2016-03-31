#ifndef __TACHO_EXAMPLE_HIER_MATRIX_HPP__
#define __TACHO_EXAMPLE_HIER_MATRIX_HPP__

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

namespace Tacho {

  template<typename DeviceSpaceType>
  int exampleHierMatrix(const std::string file_input,
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

    // for simple test, let's use host space only here

    typedef CrsMatrixBase<value_type,ordinal_type,size_type,HostSpaceType> CrsMatrixBaseHostType;
    typedef CrsMatrixView<CrsMatrixBaseHostType> CrsMatrixViewHostType;

    typedef GraphTools<ordinal_type,size_type,HostSpaceType> GraphToolsHostType;

    typedef GraphTools_Scotch<ordinal_type,size_type,HostSpaceType> GraphToolsHostType_Scotch;
    typedef GraphTools_CAMD<ordinal_type,size_type,HostSpaceType> GraphToolsHostType_CAMD;

    typedef SymbolicFactorization<CrsMatrixBaseHostType> SymbolicFactorizationType;

    typedef CrsMatrixBase<CrsMatrixViewHostType,ordinal_type,size_type,HostSpaceType> CrsHierBaseHostType;

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

    timer.reset();
    CrsHierBaseHostType HA_host("HA_host");
    CrsMatrixTools::createHierMatrix(HA_host, DD_host, 
                                     S.NumBlocks(),
                                     S.RangeVector(),
                                     S.TreeVector());
    double t_hier = timer.seconds();    

    // block preprocessing
    timer.reset();
    {
      Kokkos::View<ordinal_type*,HostSpaceType> nrowview_blocks("NumRowViewInBlocks", HA_host.NumNonZeros() + 1);
      
      nrowview_blocks(0) = 0;
      for (ordinal_type i=0;i<HA_host.NumNonZeros();++i) 
        nrowview_blocks(i+1) = nrowview_blocks(i) + HA_host.Value(i).NumRows();
      
      Kokkos::View<typename CrsMatrixViewHostType::row_view_type*,HostSpaceType>
        rowview_blocks("RowViewInBlocks", nrowview_blocks(HA_host.NumNonZeros()));

      typedef Kokkos::pair<ordinal_type,ordinal_type> range_type;
      Kokkos::parallel_for(Kokkos::RangePolicy<HostSpaceType>(0, HA_host.NumNonZeros()),
                           [&](const size_type i) {
                             const auto begin = nrowview_blocks(i);
                             const auto end   = nrowview_blocks(i+1);
                             HA_host.Value(i).setRowViewArray
                               (Kokkos::subview(rowview_blocks, range_type(begin, end)));
                           } );

      CrsMatrixTools::filterEmptyBlocks(HA_host);
    }
    double t_block = timer.seconds();    

    size_type nnz_blocks = 0;
    for (size_type i=0; i<HA_host.NumNonZeros(); ++i) {
      nnz_blocks += HA_host.Value(i).NumNonZeros();
      HA_host.Value(i).showMe(std::cout) << std::endl;
    }

    if (nnz_blocks != DD_host.NumNonZeros()) 
      std::cout << "Error, nnz blocks = " << nnz_blocks << " vs flat nnz = " << DD_host.NumNonZeros() 
                << std::endl;
    else 
      std::cout << "Pass, nnz blocks match to flat nnz" << std::endl;      

    {
      const auto prec = std::cout.precision();
      std::cout.precision(4);

      std::cout << std::scientific;
      std::cout << "HierMatrix:: Given matrix  dimension = " << AA_host.NumRows() << " x " << AA_host.NumCols()
                << ", " << " nnz = " << AA_host.NumNonZeros() << std::endl;
      std::cout << "HierMatrix:: Upper factors dimension = " << DD_host.NumRows() << " x " << DD_host.NumCols()
                << ", " << " nnz = " << DD_host.NumNonZeros() << std::endl;
      std::cout << "HierMatrix:: Hier matrix   dimension = " << HA_host.NumRows() << " x " << HA_host.NumCols()
                << ", " << " nnz = " << HA_host.NumNonZeros() << std::endl;

      std::cout << "HierMatrix:: "
                << "read = " << t_read << " [sec], "
                << "graph generation = " << (t_graph/2.0) << " [sec] "
                << "scotch reordering = " << t_scotch << " [sec] "
                << "camd reordering = " << t_camd << " [sec] " << std::endl
                << "HierMatrix:: "
                << "symbolic factorization = " << t_symbolic << " [sec] "
                << "hier creation = " << t_hier << " [sec] "
                << "block specification = " << t_block << " [sec] "
                << std::endl;

      std::cout.unsetf(std::ios::scientific);
      std::cout.precision(prec);
    }

    return r_val;
  }
}

#endif
