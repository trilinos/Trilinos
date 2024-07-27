// clang-format off
// @HEADER
// *****************************************************************************
//                            Tacho package
//
// Copyright 2022 NTESS and the Tacho contributors.
// SPDX-License-Identifier: BSD-2-Clause
// *****************************************************************************
// @HEADER
// clang-format on
#ifndef __TACHO_TRSM_HPP__
#define __TACHO_TRSM_HPP__

/// \file Tacho_Trsm.hpp
/// \author Kyungjoo Kim (kyukim@sandia.gov)

#include "Tacho_Util.hpp"

namespace Tacho {

///
/// Trsm:
///

/// various implementation for different uplo and algo parameters
template <typename ArgSide, typename ArgUplo, typename ArgTrans, typename ArgAlgo> struct Trsm;

struct TrsmAlgorithm {
  using type = ActiveAlgorithm<runsOnCudaOrHIP()>::type;
};

struct TrsmAlgorithm_Team {
#if defined(KOKKOS_ENABLE_OPENMP)
  using type = ActiveHostAlgorithm<runsWithOMP()>::type;
#else
  using type = ActiveAlgorithm<runsOnCudaOrHIP()>::type;
#endif
};

} // namespace Tacho

#endif
