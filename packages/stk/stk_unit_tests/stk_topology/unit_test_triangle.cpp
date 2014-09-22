// Copyright (c) 2013, Sandia Corporation.
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Governement retains certain rights in this software.
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

TEST( stk_topology, tri_3)
{
  using stk::topology;

  topology t = topology::TRI_3;

  EXPECT_TRUE(t.is_valid());
  EXPECT_FALSE(t.has_homogeneous_faces());
  EXPECT_FALSE(t.is_shell());

  EXPECT_EQ(t.rank(),topology::FACE_RANK);
  EXPECT_EQ(t.side_rank(),topology::EDGE_RANK);

  EXPECT_EQ(t.dimension(),2u);
  EXPECT_EQ(t.num_nodes(),3u);
  EXPECT_EQ(t.num_vertices(),3u);
  EXPECT_EQ(t.num_edges(),3u);
  EXPECT_EQ(t.num_faces(),0u);
  EXPECT_EQ(t.num_permutations(),6u);
  EXPECT_EQ(t.num_positive_permutations(),3u);

  EXPECT_FALSE(t.defined_on_spatial_dimension(1));
  EXPECT_FALSE(t.defined_on_spatial_dimension(2));
  EXPECT_TRUE(t.defined_on_spatial_dimension(3));

  EXPECT_EQ(t.base(),topology::TRI_3);

  const char a[] = "abc";

  {
    const char b[] = "abc";
    EXPECT_TRUE(t.equivalent(a,b).first);
    EXPECT_EQ(t.equivalent(a,b).second,0u);
  }

  {
    const char b[] = "cab";
    EXPECT_TRUE(t.equivalent(a,b).first);
    EXPECT_EQ(t.equivalent(a,b).second,1u);
  }

  {
    const char b[] = "bca";
    EXPECT_TRUE(t.equivalent(a,b).first);
    EXPECT_EQ(t.equivalent(a,b).second,2u);
  }

  {
    const char b[] = "acb";
    EXPECT_TRUE(t.equivalent(a,b).first);
    EXPECT_EQ(t.equivalent(a,b).second,3u);
  }

  {
    const char b[] = "cba";
    EXPECT_TRUE(t.equivalent(a,b).first);
    EXPECT_EQ(t.equivalent(a,b).second,4u);
  }

  {
    const char b[] = "bac";
    EXPECT_TRUE(t.equivalent(a,b).first);
    EXPECT_EQ(t.equivalent(a,b).second,5u);
  }
}

TEST( stk_topology, tri_4)
{
  using stk::topology;

  topology t = topology::TRI_4;

  EXPECT_TRUE(t.is_valid());
  EXPECT_FALSE(t.has_homogeneous_faces());
  EXPECT_FALSE(t.is_shell());

  EXPECT_EQ(t.rank(),topology::FACE_RANK);
  EXPECT_EQ(t.side_rank(),topology::EDGE_RANK);

  EXPECT_EQ(t.dimension(),2u);
  EXPECT_EQ(t.num_nodes(),4u);
  EXPECT_EQ(t.num_vertices(),3u);
  EXPECT_EQ(t.num_edges(),3u);
  EXPECT_EQ(t.num_faces(),0u);
  EXPECT_EQ(t.num_permutations(),6u);
  EXPECT_EQ(t.num_positive_permutations(),3u);

  EXPECT_FALSE(t.defined_on_spatial_dimension(1));
  EXPECT_FALSE(t.defined_on_spatial_dimension(2));
  EXPECT_TRUE(t.defined_on_spatial_dimension(3));

  EXPECT_EQ(t.base(),topology::TRI_3);
  const char a[] = "abcd";

  {
    const char b[] = "abcd";
    EXPECT_TRUE(t.equivalent(a,b).first);
    EXPECT_EQ(t.equivalent(a,b).second,0u);
  }

  {
    const char b[] = "cabd";
    EXPECT_TRUE(t.equivalent(a,b).first);
    EXPECT_EQ(t.equivalent(a,b).second,1u);
  }

  {
    const char b[] = "bcad";
    EXPECT_TRUE(t.equivalent(a,b).first);
    EXPECT_EQ(t.equivalent(a,b).second,2u);
  }

  {
    const char b[] = "acbd";
    EXPECT_TRUE(t.equivalent(a,b).first);
    EXPECT_EQ(t.equivalent(a,b).second,3u);
  }

  {
    const char b[] = "cbad";
    EXPECT_TRUE(t.equivalent(a,b).first);
    EXPECT_EQ(t.equivalent(a,b).second,4u);
  }

  {
    const char b[] = "bacd";
    EXPECT_TRUE(t.equivalent(a,b).first);
    EXPECT_EQ(t.equivalent(a,b).second,5u);
  }
}

TEST( stk_topology, tri_6)
{
  using stk::topology;

  topology t = topology::TRI_6;

  EXPECT_TRUE(t.is_valid());
  EXPECT_FALSE(t.has_homogeneous_faces());
  EXPECT_FALSE(t.is_shell());

  EXPECT_EQ(t.rank(),topology::FACE_RANK);
  EXPECT_EQ(t.side_rank(),topology::EDGE_RANK);

  EXPECT_EQ(t.dimension(),2u);
  EXPECT_EQ(t.num_nodes(),6u);
  EXPECT_EQ(t.num_vertices(),3u);
  EXPECT_EQ(t.num_edges(),3u);
  EXPECT_EQ(t.num_faces(),0u);
  EXPECT_EQ(t.num_permutations(),6u);
  EXPECT_EQ(t.num_positive_permutations(),3u);

  EXPECT_FALSE(t.defined_on_spatial_dimension(1));
  EXPECT_FALSE(t.defined_on_spatial_dimension(2));
  EXPECT_TRUE(t.defined_on_spatial_dimension(3));

  const char a[] = "abc012";

  {
    const char b[] = "abc012";
    EXPECT_TRUE(t.equivalent(a,b).first);
    EXPECT_EQ(t.equivalent(a,b).second,0u);
  }

  {
    const char b[] = "cab201";
    EXPECT_TRUE(t.equivalent(a,b).first);
    EXPECT_EQ(t.equivalent(a,b).second,1u);
  }

  {
    const char b[] = "bca120";
    EXPECT_TRUE(t.equivalent(a,b).first);
    EXPECT_EQ(t.equivalent(a,b).second,2u);
  }

  {
    const char b[] = "acb210";
    EXPECT_TRUE(t.equivalent(a,b).first);
    EXPECT_EQ(t.equivalent(a,b).second,3u);
  }

  {
    const char b[] = "cba102";
    EXPECT_TRUE(t.equivalent(a,b).first);
    EXPECT_EQ(t.equivalent(a,b).second,4u);
  }

  {
    const char b[] = "bac021";
    EXPECT_TRUE(t.equivalent(a,b).first);
    EXPECT_EQ(t.equivalent(a,b).second,5u);
  }
}

