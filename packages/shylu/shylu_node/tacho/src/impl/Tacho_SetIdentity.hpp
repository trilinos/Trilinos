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
#ifndef __TACHO_SET_IDENTITY_HPP__
#define __TACHO_SET_IDENTITY_HPP__

/// \file Tacho_SetIdentity.hpp
/// \brief Front interface for SetIdentity
/// \author Kyungjoo Kim (kyukim@sandia.gov)

#include "Tacho_Util.hpp"

namespace Tacho {

///
/// SetIdentity
///

/// various implementation for different uplo and algo parameters
template <typename ArgAlgo> struct SetIdentity;

} // namespace Tacho

#endif
