// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef PANZER_DIMENSION_HPP
#define PANZER_DIMENSION_HPP

#include "Phalanx_ExtentTraits.hpp"

namespace panzer {
  /*
  //! Spatial Dimension Tag
  SHARDS_ARRAY_DIM_TAG_SIMPLE_DECLARATION(Dim)
  //! Integration Point Tag
  SHARDS_ARRAY_DIM_TAG_SIMPLE_DECLARATION(IP)
  //! Basis Point Tag (generalization of the NODE)
  SHARDS_ARRAY_DIM_TAG_SIMPLE_DECLARATION(BASIS)

  //! Node Point Tag
  // SHARDS_ARRAY_DIM_TAG_SIMPLE_DECLARATION(NODE)
  using NODE = BASIS;

  //! Generic Point Tag
  SHARDS_ARRAY_DIM_TAG_SIMPLE_DECLARATION(Point)
  //! Cell Tag
  SHARDS_ARRAY_DIM_TAG_SIMPLE_DECLARATION(Cell)
  //! Face Tag
  SHARDS_ARRAY_DIM_TAG_SIMPLE_DECLARATION(Face)
  //! Dummy Tag
  SHARDS_ARRAY_DIM_TAG_SIMPLE_DECLARATION(Dummy)
  //! Edge Tag
  SHARDS_ARRAY_DIM_TAG_SIMPLE_DECLARATION(Edge)  // Suzey: 06/11/201
  */

  struct Dim{};
  struct IP{};
  struct BASIS{};
  using NODE = BASIS;
  struct Point{};
  struct Cell{};
  struct Face{};
  struct Edge{};
  struct Dummy{};
}

namespace PHX {
  template<> std::string print<panzer::Dim>();
  template<> std::string print<panzer::IP>();
  template<> std::string print<panzer::BASIS>();
  template<> std::string print<panzer::Point>();
  template<> std::string print<panzer::Cell>();
  template<> std::string print<panzer::Face>();
  template<> std::string print<panzer::Edge>();
  template<> std::string print<panzer::Dummy>();
}

PHX_IS_EXTENT(panzer::Dim)
PHX_IS_EXTENT(panzer::IP)
PHX_IS_EXTENT(panzer::BASIS)
PHX_IS_EXTENT(panzer::Point)
PHX_IS_EXTENT(panzer::Cell)
PHX_IS_EXTENT(panzer::Face)
PHX_IS_EXTENT(panzer::Edge)
PHX_IS_EXTENT(panzer::Dummy)

#endif
