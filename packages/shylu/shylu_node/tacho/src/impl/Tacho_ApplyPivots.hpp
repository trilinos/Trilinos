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
#ifndef __TACHO_APPLY_PIVOTS_HPP__
#define __TACHO_APPLY_PIVOTS_HPP__

/// \file Tacho_ApplyPivots.hpp
/// \brief Front interface to apply pivots
/// \author Kyungjoo Kim (kyukim@sandia.gov)

#include "Tacho_Util.hpp"

namespace Tacho {

///
/// Apply Pivots
///

/// various implementation for different uplo and algo parameters
template <typename ArgPivotMode, typename ArgSide, typename ArgDirect, typename ArgAlgo> struct ApplyPivots;

} // namespace Tacho

#endif
