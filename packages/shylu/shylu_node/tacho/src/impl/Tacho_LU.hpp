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
  using type = ActiveAlgorithm<runsOnCudaOrHIP()>::type;
};

struct LU_Algorithm_Team {
#if defined(KOKKOS_ENABLE_OPENMP)
  using type = ActiveHostAlgorithm<runsWithOMP()>::type;
#else
  using type = ActiveAlgorithm<runsOnCudaOrHIP()>::type;
#endif
};

} // namespace Tacho

#endif
