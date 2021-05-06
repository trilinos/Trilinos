#ifndef __TACHO_APPLY_PERMUTATION_HPP__
#define __TACHO_APPLY_PERMUTATION_HPP__

/// \file Tacho_ApplyPermutation.hpp
/// \brief Front interface to apply permutation
/// \author Kyungjoo Kim (kyukim@sandia.gov)

#include "Tacho_Util.hpp"

namespace Tacho {

  ///
  /// Apply Permutation
  ///

  /// various implementation for different uplo and algo parameters
  template<typename ArgSide, typename ArgTrans, typename ArgAlgo>
  struct ApplyPermutation;
}

#endif
