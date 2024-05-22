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

#include <gtest/gtest.h>
#include <stk_util/stk_config.h>
#include <stk_mesh/base/MetaData.hpp>
#include <Shards_BasicTopologies.hpp>

namespace {

void test_mapping(stk::topology stkTopo, bool expectInvalid = false)
{
  shards::CellTopology cellTopo = stk::mesh::get_cell_topology(stkTopo);

  if (stkTopo == stk::topology::INVALID_TOPOLOGY || expectInvalid) {
    EXPECT_FALSE(cellTopo.isValid()) << stkTopo;
    return;
  }

  EXPECT_TRUE(cellTopo.isValid()) << stkTopo;
  EXPECT_EQ(cellTopo.getNodeCount(), stkTopo.num_nodes()) << stkTopo;
  EXPECT_EQ(cellTopo.getVertexCount(), stkTopo.num_vertices()) << stkTopo;
}

TEST( GetCellTopology, mappingTest )
{
  constexpr bool expectInvalid = true;

  test_mapping(stk::topology::NODE);
  test_mapping(stk::topology::LINE_2);
  test_mapping(stk::topology::LINE_3);
  test_mapping(stk::topology::TRI_3);
  test_mapping(stk::topology::TRI_4);
  test_mapping(stk::topology::TRI_6);
  test_mapping(stk::topology::QUAD_4);
  test_mapping(stk::topology::QUAD_6, expectInvalid);
  test_mapping(stk::topology::QUAD_8);
  test_mapping(stk::topology::QUAD_9);
  test_mapping(stk::topology::PARTICLE);
  test_mapping(stk::topology::LINE_2_1D);
  test_mapping(stk::topology::LINE_3_1D);
  test_mapping(stk::topology::BEAM_2);
  test_mapping(stk::topology::BEAM_3);
  test_mapping(stk::topology::SHELL_LINE_2);
  test_mapping(stk::topology::SHELL_LINE_3);
  test_mapping(stk::topology::SPRING_2, expectInvalid);
  test_mapping(stk::topology::SPRING_3, expectInvalid);
  test_mapping(stk::topology::TRI_3_2D);
  test_mapping(stk::topology::TRI_4_2D);
  test_mapping(stk::topology::TRI_6_2D);
  test_mapping(stk::topology::QUAD_4_2D);
  test_mapping(stk::topology::QUAD_8_2D);
  test_mapping(stk::topology::QUAD_9_2D);
  test_mapping(stk::topology::SHELL_TRI_3);
  test_mapping(stk::topology::SHELL_TRI_4, expectInvalid);
  test_mapping(stk::topology::SHELL_TRI_6);
  test_mapping(stk::topology::SHELL_QUAD_4);
  test_mapping(stk::topology::SHELL_QUAD_8);
  test_mapping(stk::topology::SHELL_QUAD_9);
  test_mapping(stk::topology::TET_4);
  test_mapping(stk::topology::TET_8);
  test_mapping(stk::topology::TET_10);
  test_mapping(stk::topology::TET_11);
  test_mapping(stk::topology::PYRAMID_5);
  test_mapping(stk::topology::PYRAMID_13);
  test_mapping(stk::topology::PYRAMID_14);
  test_mapping(stk::topology::WEDGE_6);
  test_mapping(stk::topology::WEDGE_12, expectInvalid);
  test_mapping(stk::topology::WEDGE_15);
  test_mapping(stk::topology::WEDGE_18);
  test_mapping(stk::topology::HEX_8);
  test_mapping(stk::topology::HEX_20);
  test_mapping(stk::topology::HEX_27);
}

} // empty namespace
