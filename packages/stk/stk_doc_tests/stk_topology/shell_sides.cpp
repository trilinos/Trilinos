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

#include "gtest/gtest.h"
#include "stk_topology/topology.hpp"

namespace {

//begin_num_sides
TEST(stk_topology, shell_side_num_sides) {
  stk::topology shell = stk::topology::SHELL_QUAD_4;

  EXPECT_TRUE(shell.is_valid());
  EXPECT_TRUE(shell.is_shell());

  EXPECT_EQ(shell.rank(),stk::topology::ELEMENT_RANK);
  EXPECT_EQ(shell.side_rank(),stk::topology::FACE_RANK);

  EXPECT_EQ(shell.num_vertices(),4u);
  EXPECT_EQ(shell.num_edges(),4u);

  EXPECT_EQ(shell.num_faces(),2u);
  EXPECT_EQ(shell.num_sides(),6u);
  EXPECT_EQ(shell.num_sub_topology(shell.side_rank()), 2u);
  EXPECT_NE(shell.num_sub_topology(shell.side_rank()), shell.num_sides());
}
//end_num_sides

//begin_shell_side_topo
TEST(stk_topology, shell_side_topology) {
  stk::topology shell = stk::topology::SHELL_QUAD_4;

  EXPECT_TRUE(shell.is_valid());
  EXPECT_TRUE(shell.is_shell());

  EXPECT_EQ(shell.num_faces(),2u);
  EXPECT_EQ(shell.face_topology(0), stk::topology::QUAD_4);
  EXPECT_EQ(shell.face_topology(1), stk::topology::QUAD_4);

  EXPECT_EQ(shell.num_sides(),6u);
  EXPECT_EQ(shell.side_topology(0), stk::topology::QUAD_4);
  EXPECT_EQ(shell.side_topology(1), stk::topology::QUAD_4);
  EXPECT_EQ(shell.side_topology(2), stk::topology::LINE_2);
  EXPECT_EQ(shell.side_topology(3), stk::topology::LINE_2);
  EXPECT_EQ(shell.side_topology(4), stk::topology::LINE_2);
  EXPECT_EQ(shell.side_topology(5), stk::topology::LINE_2);
}
//end_shell_side_topo

}