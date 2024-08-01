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
#ifndef __TACHO_SCALE_2X2_BLOCK_INVERSE_DIAGONALS_HPP__
#define __TACHO_SCALE_2X2_BLOCK_INVERSE_DIAGONALS_HPP__

/// \file Tacho_Scale2x2_BlockInverseDiagonals.hpp
/// \brief Scale matrix with given diagonal blocks
/// \author Kyungjoo Kim (kyukim@sandia.gov)

#include "Tacho_Util.hpp"

namespace Tacho {

///
/// Scale2x2_BlockInverseDiagonals
///

/// various implementation for different uplo and algo parameters
template <typename ArgSide, typename ArgAlgo> struct Scale2x2_BlockInverseDiagonals;
} // namespace Tacho

#endif
