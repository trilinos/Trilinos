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
#ifndef __TACHO_COPY_HPP__
#define __TACHO_COPY_HPP__

/// \file Tacho_Copy.hpp
/// \brief Front interface for Copy
/// \author Kyungjoo Kim (kyukim@sandia.gov)

#include "Tacho_Util.hpp"

namespace Tacho {

///
/// Copy
///

/// various implementation for different uplo and algo parameters
template <typename ArgAlgo> struct Copy;
} // namespace Tacho

#endif
