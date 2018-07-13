#ifndef __TACHO_EXAMPLE_INCOMPLETE_SYMBOLIC_FACTORIZATION_HPP__
#define __TACHO_EXAMPLE_INCOMPLETE_SYMBOLIC_FACTORIZATION_HPP__

#include <Kokkos_Core.hpp>
#include <impl/Kokkos_Timer.hpp>

#include "Tacho_Util.hpp"

#include "Tacho_CrsMatrixBase.hpp"
#include "Tacho_CrsMatrixTools.hpp"

#include "Tacho_MatrixMarket.hpp"

#include "Tacho_GraphTools.hpp"

#include "Tacho_GraphTools_Scotch.hpp"
#include "Tacho_GraphTools_CAMD.hpp"

#include "Tacho_IncompleteSymbolicFactorization.hpp"

namespace Tacho {
  
  template<typename DeviceSpaceType>
  int exampleIncompleteSymbolicFactorization(const std::string file_input,
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

    typedef IncompleteSymbolicFactorization<CrsMatrixBaseHostType> IncompleteSymbolicFactorizationType;

    int r_val = 0;

    Kokkos::Impl::Timer timer;
    
    CrsMatrixBaseHostType AA("AA");
    timer.reset();
    {
      std::ifstream in;
      in.open(file_input);
      if (!in.good()) {
        std::cout << "Failed in open the file: " << file_input << std::endl;
        return -1;
      }
      MatrixMarket::read(AA, in);
    }
    double t_read = timer.seconds();

    if (verbose)
      AA.showMe(std::cout) << std::endl;
    
    typename GraphToolsHostType::size_type_array rptr("Graph::RowPtrArray", AA.NumRows() + 1);
    typename GraphToolsHostType::ordinal_type_array cidx("Graph::ColIndexArray", AA.NumNonZeros());

    timer.reset();
    GraphToolsHostType::getGraph(rptr, cidx, AA);
    double t_graph = timer.seconds();

    GraphToolsHostType_Scotch S;
    S.setGraph(AA.NumRows(), rptr, cidx);
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

    std::cout << "Using tree level ( " << (S.TreeLevel() - treecut)<< " ), Scotch # of ranges = " << S.NumBlocks() << std::endl;
    S.pruneTree(prunecut);
    std::cout << "After prunecut level ( " << prunecut << " ),  # of ranges = " << S.NumBlocks()  << std::endl;

    if (verbose)
      S.showMe(std::cout) << std::endl;

    CrsMatrixBaseHostType BB("BB");
    BB.createConfTo(AA);

    CrsMatrixTools::copy(BB, 
                         S.PermVector(),
                         S.InvPermVector(),
                         AA);

    if (verbose)
      BB.showMe(std::cout) << std::endl;

    timer.reset();
    GraphToolsHostType::getGraph(rptr, cidx, BB);
    t_graph += timer.seconds();

    GraphToolsHostType_CAMD C;
    C.setGraph(BB.NumRows(), 
               rptr, cidx, 
               S.NumBlocks(), 
               S.RangeVector());

    timer.reset();
    C.computeOrdering();
    double t_camd = timer.seconds();

    if (verbose)
      C.showMe(std::cout) << std::endl;

    CrsMatrixBaseHostType CC("CC");
    CC.createConfTo(BB);

    CrsMatrixTools::copy(CC, 
                         C.PermVector(),
                         C.InvPermVector(),
                         BB);

    if (verbose)
      CC.showMe(std::cout) << std::endl;

    CrsMatrixBaseHostType DD("DD");    
    
    timer.reset();
    IncompleteSymbolicFactorizationType::createNonZeroPattern(DD, 
                                                    fill_level, 
                                                    Uplo::Upper,
                                                    CC,
                                                    rows_per_team);
    double t_symbolic = timer.seconds();

    if (verbose)
      DD.showMe(std::cout) << std::endl;
    
    {
      const auto prec = std::cout.precision();
      std::cout.precision(4);

      std::cout << std::scientific;
      std::cout << "SymbolicFactorization:: Given matrix  dimension = " << AA.NumRows() << " x " << AA.NumCols()
                << ", " << " nnz = " << AA.NumNonZeros() << std::endl;
      std::cout << "SymbolicFactorization:: Upper factors dimension = " << DD.NumRows() << " x " << DD.NumCols()
                << ", " << " nnz = " << DD.NumNonZeros() << std::endl;
      
      std::cout << "SymbolicFactorization:: " 
                << "read = " << t_read << " [sec], "
                << "graph generation = " << (t_graph/2.0) << " [sec] "
                << "scotch reordering = " << t_scotch << " [sec] "
                << "camd reordering = " << t_camd << " [sec] "
                << "symbolic factorization = " << t_symbolic << " [sec] "
                << std::endl;
      
      std::cout.unsetf(std::ios::scientific);
      std::cout.precision(prec);
    }

    return r_val;
  }
}

#endif
