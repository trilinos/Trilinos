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
#ifndef __TACHO_GEMM_TRIANGULAR_HPP__
#define __TACHO_GEMM_TRIANGULAR_HPP__

/// \file Tacho_GemmTriangular.hpp
/// \brief Update the upper/lower triangular of C, C must be a square matrix.
/// \author Kyungjoo Kim (kyukim@sandia.gov)

#include "Tacho_Util.hpp"

namespace Tacho {

///
/// GemmTriangular:
///

/// various implementation for different uplo and algo parameters
template <typename ArgTransA, typename ArgTransB, typename ArgUploC, typename ArgAlgo> struct GemmTriangular;

} // namespace Tacho

#endif
