#ifndef __TACHO_CHOL_HPP__
#define __TACHO_CHOL_HPP__

/// \file Tacho_Chol.hpp
/// \brief Front interface for Cholesky dense factorization
/// \author Kyungjoo Kim (kyukim@sandia.gov)

#include "Tacho_Util.hpp"

namespace Tacho {

///
/// Chol:
///
///

/// various implementation for different uplo and algo parameters
template <typename ArgUplo, typename ArgAlgo> struct Chol;

struct CholAlgorithm {
  using type = ActiveAlgorithm::type;
};

} // namespace Tacho

#endif
