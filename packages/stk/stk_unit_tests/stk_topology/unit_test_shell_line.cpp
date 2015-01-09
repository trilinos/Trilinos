// Copyright (c) 2013, Sandia Corporation.
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
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
//     * Neither the name of Sandia Corporation nor the names of its
//       contributors may be used to endorse or promote products derived
//       from this software without specific prior written permission.
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
#include <stk_topology/topology.hpp>

TEST( stk_topology, shell_line_2)
{
  using stk::topology;

  topology t = topology::SHELL_LINE_2;


  EXPECT_TRUE(t.is_valid());
  EXPECT_FALSE(t.has_homogeneous_faces());
  EXPECT_TRUE(t.is_shell());

  EXPECT_EQ(t.rank(),topology::ELEMENT_RANK);
  EXPECT_EQ(t.side_rank(),topology::EDGE_RANK);


  EXPECT_EQ(t.dimension(),2u);
  EXPECT_EQ(t.num_nodes(),2u);
  EXPECT_EQ(t.num_vertices(),2u);
  EXPECT_EQ(t.num_edges(),2u);
  EXPECT_EQ(t.num_faces(),0u);
  EXPECT_EQ(t.num_permutations(),2u);
  EXPECT_EQ(t.num_positive_permutations(),1u);

  EXPECT_FALSE(t.defined_on_spatial_dimension(1));
  EXPECT_TRUE(t.defined_on_spatial_dimension(2));
  EXPECT_FALSE(t.defined_on_spatial_dimension(3));

  EXPECT_EQ(t.base(),topology::SHELL_LINE_2);

  const char a[] = "ab";
  {
    const char b[] = "ab";
    EXPECT_TRUE(t.equivalent(a,b).first);
    EXPECT_EQ(t.equivalent(a,b).second,0u);
  }

  {
    const char b[] = "ba";
    EXPECT_TRUE(t.equivalent(a,b).first);
    EXPECT_EQ(t.equivalent(a,b).second,1u);
  }

  {
    char edge[2];
    t.edge_nodes(a,0,edge);
    t.edge_topology().equivalent(edge,"ab");

    t.edge_nodes(a,1,edge);
    t.edge_topology().equivalent(edge,"ba");
  }
}

TEST( stk_topology, shell_line_3)
{
  using stk::topology;

  topology t = topology::SHELL_LINE_3;


  EXPECT_TRUE(t.is_valid());
  EXPECT_FALSE(t.has_homogeneous_faces());
  EXPECT_TRUE(t.is_shell());

  EXPECT_EQ(t.rank(),topology::ELEMENT_RANK);
  EXPECT_EQ(t.side_rank(),topology::EDGE_RANK);


  EXPECT_EQ(t.dimension(),2u);
  EXPECT_EQ(t.num_nodes(),3u);
  EXPECT_EQ(t.num_vertices(),2u);
  EXPECT_EQ(t.num_edges(),2u);
  EXPECT_EQ(t.num_faces(),0u);
  EXPECT_EQ(t.num_permutations(),2u);
  EXPECT_EQ(t.num_positive_permutations(),1u);

  EXPECT_FALSE(t.defined_on_spatial_dimension(1));
  EXPECT_TRUE(t.defined_on_spatial_dimension(2));
  EXPECT_FALSE(t.defined_on_spatial_dimension(3));

  EXPECT_EQ(t.base(),topology::SHELL_LINE_2);

  const char a[] = "abc";
  {
    const char b[] = "abc";
    EXPECT_TRUE(t.equivalent(a,b).first);
    EXPECT_EQ(t.equivalent(a,b).second,0u);
  }

  {
    const char b[] = "bac";
    EXPECT_TRUE(t.equivalent(a,b).first);
    EXPECT_EQ(t.equivalent(a,b).second,1u);
  }

  {
    char edge[3];
    t.edge_nodes(a,0,edge);
    t.edge_topology().equivalent(edge,"abc");

    t.edge_nodes(a,1,edge);
    t.edge_topology().equivalent(edge,"bac");
  }
}

