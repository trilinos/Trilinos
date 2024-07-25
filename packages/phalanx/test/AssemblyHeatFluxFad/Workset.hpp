// @HEADER
// *****************************************************************************
//        Phalanx: A Partial Differential Equation Field Evaluation 
//       Kernel for Flexible Management of Complex Dependency Chains
//
// Copyright 2008 NTESS and the Phalanx contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef PHX_EXAMPLE_MY_WORKSET_HPP
#define PHX_EXAMPLE_MY_WORKSET_HPP

#include "Phalanx_config.hpp" // for std::vector
#include "Cell.hpp"

struct MyWorkset {
  
  std::size_t local_offset;

  std::size_t num_cells;
  
  std::vector<MyCell>::iterator begin;

  std::vector<MyCell>::iterator end;

};

#endif
