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
#ifndef __TACHO_LDL_HPP__
#define __TACHO_LDL_HPP__

/// \file Tacho_LDL.hpp
/// \brief Front interface for LDL^t dense factorization
/// \author Kyungjoo Kim (kyukim@sandia.gov)

#include "Tacho_Util.hpp"

namespace Tacho {

///
/// LDL:
///
///

/// various implementation for different uplo and algo parameters
template <typename ArgUplo, typename ArgAlgo> struct LDL;

struct LDL_Algorithm {
  using type = ActiveAlgorithm<runsOnCudaOrHIP()>::type;
};

struct LDL_Algorithm_Team {
#if defined(KOKKOS_ENABLE_OPENMP)
  using type = ActiveHostAlgorithm<runsWithOMP()>::type;
#else
  using type = ActiveAlgorithm<runsOnCudaOrHIP()>::type;
#endif
};
} // namespace Tacho

#endif
