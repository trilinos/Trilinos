#ifndef __TACHO_APPLY_PIVOTS_HPP__
#define __TACHO_APPLY_PIVOTS_HPP__

/// \file Tacho_ApplyPivots.hpp
/// \brief Front interface to apply pivots 
/// \author Kyungjoo Kim (kyukim@sandia.gov)

#include "Tacho_Util.hpp"

namespace Tacho {

  ///
  /// Apply Pivots
  ///

  /// various implementation for different uplo and algo parameters
  template<typename ArgPivotMode, typename ArgSide, typename ArgDirect, typename ArgAlgo>
  struct ApplyPivots;
}

#endif
