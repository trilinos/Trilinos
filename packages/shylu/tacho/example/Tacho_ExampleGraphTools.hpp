#ifndef __TACHO_EXAMPLE_GRAPH_TOOLS_HPP__
#define __TACHO_EXAMPLE_GRAPH_TOOLS_HPP__

#include <Kokkos_Core.hpp>
#include <impl/Kokkos_Timer.hpp>

#include "ShyLUTacho_config.h"

#include "Tacho_Util.hpp"

#include "Tacho_CrsMatrixBase.hpp"
#include "Tacho_CrsMatrixTools.hpp"

#include "Tacho_MatrixMarket.hpp"

#include "Tacho_GraphTools.hpp"
//#include "Tacho_GraphTools_Scotch.hpp"
//#include "Tacho_GraphTools_CAMD.hpp"

namespace Tacho {
  
  template<typename DeviceSpaceType>
  int exampleGraphTools(const std::string file_input,
                        const bool verbose) {
    typedef typename
      Kokkos::Impl::is_space<DeviceSpaceType>::host_mirror_space::execution_space HostSpaceType ;

    const bool detail = false;
    std::cout << "DeviceSpace::  "; DeviceSpaceType::print_configuration(std::cout, detail);
    std::cout << "HostSpace::    ";   HostSpaceType::print_configuration(std::cout, detail);

    typedef CrsMatrixBase<value_type,ordinal_type,size_type,HostSpaceType> CrsMatrixBaseHostType;
    typedef GraphTools<ordinal_type,size_type,HostSpaceType> GraphToolsHostType;

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
    
    typename GraphToolsHostType::size_type_array rptr("Graph::RowPtrArray", AA.NumRows() + 1);
    typename GraphToolsHostType::ordinal_type_array cidx("Graph::ColIndexArray", AA.NumNonZeros());

    timer.reset();
    GraphToolsHostType::getGraph(rptr, cidx, AA);
    double t_graph = timer.seconds();

    {
      const auto prec = std::cout.precision();
      std::cout.precision(4);

      std::cout << std::scientific
                << "GraphTools:: dimension = " << AA.NumRows() << " x " << AA.NumCols()
                << ", " << " nnz = " << AA.NumNonZeros() << ", "
                << "read = " << t_read << " [sec], "
                << "graph generation = " << t_graph << " [sec] "
                << std::endl;

      std::cout.unsetf(std::ios::scientific);
      std::cout.precision(prec);
    }

    return r_val;
  }
}

#endif
