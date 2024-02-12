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

#include "gtest/gtest.h"              // for AssertionResult, Test, Message, SuiteApiResolver
#include "stk_topology/topology.hpp"  // for topology, topology::END_TOPOLOGY, topology::topology_t

TEST(stk_topology, isTri)
{
  for (unsigned topology = 0; topology < stk::topology::END_TOPOLOGY; ++topology) {
    bool amITri = stk::isTriangleElement((stk::topology::topology_t)topology);
    if (topology >= stk::topology::TRI_3_2D && topology <= stk::topology::TRI_6_2D) {
      EXPECT_TRUE(amITri);
    }
    else {
      EXPECT_FALSE(amITri);
    }
  }
}


TEST(stk_topology, isQuad)
{
  for (unsigned topology = 0; topology < stk::topology::END_TOPOLOGY; ++topology) {
    bool amIQuad = stk::isQuadrilateralElement((stk::topology::topology_t)topology);
    if (topology >= stk::topology::QUAD_4_2D && topology <= stk::topology::QUAD_9_2D) {
      EXPECT_TRUE(amIQuad);
    }
    else{
      EXPECT_FALSE(amIQuad);
    }
  }
}

TEST(stk_topology, isHex)
{
  for (unsigned topology = 0; topology < stk::topology::END_TOPOLOGY; ++topology) {
    bool amIHex = stk::isHexahedronElement((stk::topology::topology_t)topology);
    if (topology >= stk::topology::HEX_8 && topology <= stk::topology::HEX_27) {
      EXPECT_TRUE(amIHex);
    }
    else {
      EXPECT_FALSE(amIHex);
    }
  }
}

TEST(stk_topology, isTet)
{
  for (unsigned topology = 0; topology < stk::topology::END_TOPOLOGY; ++topology) {
    bool amITet = stk::isTetrahedronElement((stk::topology::topology_t)topology);
    if (topology >= stk::topology::TET_4 && topology <= stk::topology::TET_11) {
      EXPECT_TRUE(amITet);
    }
    else {
      EXPECT_FALSE(amITet);
    }
  }
}

