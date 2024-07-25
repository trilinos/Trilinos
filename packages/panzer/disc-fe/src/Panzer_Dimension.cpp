// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Panzer_Dimension.hpp"

namespace panzer {

  // SHARDS_ARRAY_DIM_TAG_SIMPLE_IMPLEMENTATION(Dim)
  // SHARDS_ARRAY_DIM_TAG_SIMPLE_IMPLEMENTATION(IP)
  // SHARDS_ARRAY_DIM_TAG_SIMPLE_IMPLEMENTATION(BASIS)
  // // SHARDS_ARRAY_DIM_TAG_SIMPLE_IMPLEMENTATION(NODE)
  // SHARDS_ARRAY_DIM_TAG_SIMPLE_IMPLEMENTATION(Point)
  // SHARDS_ARRAY_DIM_TAG_SIMPLE_IMPLEMENTATION(Cell)
  // SHARDS_ARRAY_DIM_TAG_SIMPLE_IMPLEMENTATION(Face)
  // SHARDS_ARRAY_DIM_TAG_SIMPLE_IMPLEMENTATION(Dummy)
  // SHARDS_ARRAY_DIM_TAG_SIMPLE_IMPLEMENTATION(Edge)

}

namespace PHX {
  template<> std::string print<panzer::Dim>(){return "D";}
  template<> std::string print<panzer::IP>(){return "IP";}
  template<> std::string print<panzer::BASIS>(){return "B";}
  template<> std::string print<panzer::Point>(){return "P";}
  template<> std::string print<panzer::Cell>(){return "C";}
  template<> std::string print<panzer::Face>(){return "F";}
  template<> std::string print<panzer::Edge>(){return "E";}
  template<> std::string print<panzer::Dummy>(){return "";}
}
