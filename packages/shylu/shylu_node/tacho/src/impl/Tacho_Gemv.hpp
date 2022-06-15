#ifndef __TACHO_GEMV_HPP__
#define __TACHO_GEMV_HPP__

/// \file Tacho_Gemv.hpp
/// \brief Front interface for Herk operators
/// \author Kyungjoo Kim (kyukim@sandia.gov)

#include "Tacho_Util.hpp"

namespace Tacho {

///
/// Gemm:
///

/// various implementation for different uplo and algo parameters
template <typename ArgTrans, typename ArgAlgo> struct Gemv;

struct GemvAlgorithm {
  using type = ActiveAlgorithm::type;
};

} // namespace Tacho

#endif
