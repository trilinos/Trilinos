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
  using type = ActiveAlgorithm<runsOnCudaOrHIP()>::type;
};

struct GemmAlgorithm_Team {
#if defined(KOKKOS_ENABLE_OPENMP)
  using type = ActiveHostAlgorithm<runsWithOMP()>::type;
#else
  using type = ActiveAlgorithm<runsOnCudaOrHIP()>::type;
#endif
};

} // namespace Tacho

#endif
