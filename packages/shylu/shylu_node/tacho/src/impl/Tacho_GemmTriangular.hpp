#ifndef __TACHO_GEMM_TRIANGULAR_HPP__
#define __TACHO_GEMM_TRIANGULAR_HPP__

/// \file Tacho_GemmTriangular.hpp
/// \brief Update the upper/lower triangular of C, C must be a square matrix.
/// \author Kyungjoo Kim (kyukim@sandia.gov)

#include "Tacho_Util.hpp"

namespace Tacho {

  ///
  /// GemmTriangular:
  ///
  
  /// various implementation for different uplo and algo parameters
  template<typename ArgTransA, typename ArgTransB, typename ArgUploC, typename ArgAlgo>
  struct GemmTriangular;

}

#endif
