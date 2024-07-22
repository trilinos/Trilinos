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
#ifndef __TACHO_SYMMETRIZE_HPP__
#define __TACHO_SYMMETRIZE_HPP__

/// \file Tacho_Symmetrize.hpp
/// \brief Front interface for Symmetrize
/// \author Kyungjoo Kim (kyukim@sandia.gov)

#include "Tacho_Util.hpp"

namespace Tacho {

///
/// Symmetrize
///

/// various implementation for different uplo and algo parameters
template <typename ArgUplo, typename ArgAlgo> struct Symmetrize;
} // namespace Tacho

#endif
