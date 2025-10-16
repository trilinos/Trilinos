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
#ifndef __TACHO_TRMV_HPP__
#define __TACHO_TRMV_HPP__

/// \file Tacho_Trmv.hpp
/// \brief Front interface for Herk operators
/// \author Kyungjoo Kim (kyukim@sandia.gov)

#include "Tacho_Util.hpp"

namespace Tacho {

///
/// Gemm:
///

/// various implementation for different uplo and algo parameters
template <typename ArgUplo, typename ArgTrans, typename ArgAlgo> struct Trmv;

struct TrmvAlgorithm {
  using type = ActiveAlgorithm<runsOnCudaOrHIP()>::type;
};

struct TrmvAlgorithm_Team {
#if defined(KOKKOS_ENABLE_CUDA) || defined(KOKKOS_ENABLE_HIP)
  using type = ActiveAlgorithm<runsOnCudaOrHIP()>::type;
#else
  using type = ActiveHostAlgorithm<runsWithOMP()>::type;
#endif
};
} // namespace Tacho

#endif
