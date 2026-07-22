// @HEADER
// *****************************************************************************
//        Phalanx: A Partial Differential Equation Field Evaluation 
//       Kernel for Flexible Management of Complex Dependency Chains
//
// Copyright 2008 NTESS and the Phalanx contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef PHX_DIMENSION_HPP
#define PHX_DIMENSION_HPP
#include "Phalanx_ExtentTraits.hpp"

struct DIM {};
struct POINT {};
struct QP {};
struct BASIS {};
struct CELL {};

namespace PHX {
  template<> struct is_extent<DIM> : std::true_type {};
  template<> struct is_extent<POINT> : std::true_type {};
  template<> struct is_extent<QP> : std::true_type {};
  template<> struct is_extent<BASIS> : std::true_type {};
  template<> struct is_extent<CELL> : std::true_type {};
}
#endif
