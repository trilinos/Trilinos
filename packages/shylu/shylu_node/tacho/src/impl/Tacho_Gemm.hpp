#ifndef __TACHO_GEMM_HPP__
#define __TACHO_GEMM_HPP__

/// \file Tacho_Herk.hpp
/// \brief Front interface for Herk operators
/// \author Kyungjoo Kim (kyukim@sandia.gov)

#include "Tacho_Util.hpp"

namespace Tacho {

///
/// Gemm:
///

/// various implementation for different uplo and algo parameters
template <typename ArgTransA, typename ArgTransB, typename ArgAlgo> struct Gemm;

struct GemmAlgorithm {
  using type = ActiveAlgorithm::type;
};

} // namespace Tacho

#endif
