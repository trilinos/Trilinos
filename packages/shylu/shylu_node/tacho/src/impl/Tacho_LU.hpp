#ifndef __TACHO_LU_HPP__
#define __TACHO_LU_HPP__

/// \file Tacho_LU.hpp
/// \brief Front interface for LU dense factorization
/// \author Kyungjoo Kim (kyukim@sandia.gov)

#include "Tacho_Util.hpp"

namespace Tacho {

///
/// LU:
///
///

/// various implementation for different uplo and algo parameters
template <typename ArgAlgo> struct LU;

struct LU_Algorithm {
  using type = ActiveAlgorithm::type;
};

} // namespace Tacho

#endif
