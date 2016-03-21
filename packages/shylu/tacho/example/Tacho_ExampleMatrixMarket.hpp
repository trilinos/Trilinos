#ifndef __TACHO_EXAMPLE_MATRIX_MARKET_HPP__
#define __TACHO_EXAMPLE_MATRIX_MARKET_HPP__

#include <Kokkos_Core.hpp>
#include <impl/Kokkos_Timer.hpp>

#include "ShyLUTacho_config.h"

#include "Tacho_Util.hpp"
#include "Tacho_CrsMatrixBase.hpp"
#include "Tacho_CrsMatrixTools.hpp"
#include "Tacho_MatrixMarket.hpp"

namespace Tacho {
  
  template<typename DeviceSpaceType>
  int exampleMatrixMarket(const std::string file_input,
                          const bool verbose) {
    typedef typename
      Kokkos::Impl::is_space<DeviceSpaceType>::host_mirror_space::execution_space HostSpaceType ;
    
    typedef CrsMatrixBase<value_type,ordinal_type,size_type,HostSpaceType>   CrsMatrixBaseHostType;
    typedef CrsMatrixBase<value_type,ordinal_type,size_type,DeviceSpaceType> CrsMatrixBaseDeviceType;

    int r_val = 0;

    Kokkos::Impl::Timer timer;
    
    CrsMatrixBaseHostType AA("AA");
    {
      timer.reset();

      std::ifstream in;
      in.open(file_input);
      if (!in.good()) {
        std::cout << "Failed in open the file: " << file_input << std::endl;
        return -1;
      }

      MatrixMarket::import(AA, in);
    }

    std::cout << "CrsMatrixBase:: test matrices "
              <<":: m = " << AA.NumRows() << " , n = " << AA.NumCols() << " , nnz = " << AA.NumNonZeros() << std::endl;

    AA.showMe(std::cout) << std::endl;

    return r_val;
  }
}

#endif
