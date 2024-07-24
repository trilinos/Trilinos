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
  using type = ActiveAlgorithm<runsOnCudaOrHIP()>::type;
};

struct CholAlgorithm_Team {
#if defined(KOKKOS_ENABLE_OPENMP)
  using type = ActiveHostAlgorithm<runsWithOMP()>::type;
#else
  using type = ActiveAlgorithm<runsOnCudaOrHIP()>::type;
#endif
};
} // namespace Tacho

#endif
