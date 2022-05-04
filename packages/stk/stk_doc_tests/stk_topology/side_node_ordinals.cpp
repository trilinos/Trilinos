// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
// 
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
// 
//     * Redistributions in binary form must reproduce the above
//       copyright notice, this list of conditions and the following
//       disclaimer in the documentation and/or other materials provided
//       with the distribution.
// 
//     * Neither the name of NTESS nor the names of its contributors
//       may be used to endorse or promote products derived from this
//       software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
// "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
// LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
// A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
// OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
// SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
// LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
// DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
// THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
// OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
// 

#include "gtest/gtest.h"              // for Test, EXPECT_EQ, Message, TestPartResult, SuiteApiR...
#include "stk_topology/topology.hpp"  // for topology, topology::HEX_8

#ifndef __IBMCPP__

namespace {

TEST(side_node_ordinals, hex8)
{
  stk::topology hex8 = stk::topology::HEX_8;
  const unsigned nodes_per_side = 4;
  unsigned side1_node_ordinals[nodes_per_side] = {0, 0, 0, 0};

  unsigned expected_side1_node_ordinals[nodes_per_side] = {1, 2, 6, 5};

  //important usage note: side-ordinal is a 0-based number.
  //i.e., side 1 is labeled as side 2 in the exodus manual. it is the second side of the element.
  //note also that the side-node-ordinals returned are 0-based whereas the node labels in the
  //exodus manual diagrams are 1-based.

  const unsigned side1_ordinal = 1;

  //the following call to hex8.side_node_ordinals(..) fills our side1_node_ordinals array.
  hex8.side_node_ordinals(side1_ordinal, side1_node_ordinals);

  EXPECT_EQ(expected_side1_node_ordinals[0], side1_node_ordinals[0]);
  EXPECT_EQ(expected_side1_node_ordinals[1], side1_node_ordinals[1]);
  EXPECT_EQ(expected_side1_node_ordinals[2], side1_node_ordinals[2]);
  EXPECT_EQ(expected_side1_node_ordinals[3], side1_node_ordinals[3]);
}

}

#endif // __IBMCPP__
