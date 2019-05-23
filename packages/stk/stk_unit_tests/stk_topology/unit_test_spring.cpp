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
#include <stk_ngp_test/ngp_test.hpp>
#include "topology_test_utils.hpp"

TEST(stk_topology, spring_2)
{
  stk::topology t = stk::topology::SPRING_2;

  EXPECT_TRUE(t.is_valid());
  EXPECT_FALSE(t.has_homogeneous_faces());
  EXPECT_FALSE(t.is_shell());

  EXPECT_EQ(t.rank(),stk::topology::ELEMENT_RANK);
  EXPECT_EQ(t.side_rank(),stk::topology::NODE_RANK);
  EXPECT_EQ(t.num_sides(),2u);

  EXPECT_EQ(t.dimension(),2u);
  EXPECT_EQ(t.num_nodes(),2u);
  EXPECT_EQ(t.num_vertices(),2u);
  EXPECT_EQ(t.num_edges(),0u);
  EXPECT_EQ(t.num_faces(),0u);
  EXPECT_EQ(t.num_permutations(),2u);
  EXPECT_EQ(t.num_positive_permutations(),2u);

  EXPECT_TRUE(t.defined_on_spatial_dimension(1));
  EXPECT_TRUE(t.defined_on_spatial_dimension(2));
  EXPECT_TRUE(t.defined_on_spatial_dimension(3));

  EXPECT_EQ(t.base(),stk::topology::SPRING_2);

  EXPECT_EQ(t.face_topology(0), stk::topology::INVALID_TOPOLOGY);

  std::vector<std::vector<unsigned>> gold_side_node_ordinals = { {0},
                                                                 {1} };
  check_side_node_ordinals(t, gold_side_node_ordinals);
  check_side_nodes(t, gold_side_node_ordinals);

  std::vector<std::vector<unsigned>> gold_permutation_node_ordinals = { {0, 1},
                                                                        {1, 0} };
  check_permutation_node_ordinals(t, gold_permutation_node_ordinals);
  check_permutation_nodes(t, gold_permutation_node_ordinals);

  check_equivalent(t, gold_permutation_node_ordinals);
  check_lexicographical_smallest_permutation(t, gold_permutation_node_ordinals);
}

void check_spring2_on_device()
{
  Kokkos::parallel_for(1, KOKKOS_LAMBDA(const int i)
  {
    stk::topology t = stk::topology::SPRING_2;

    NGP_EXPECT_TRUE(t.is_valid());
    NGP_EXPECT_FALSE(t.has_homogeneous_faces());
    NGP_EXPECT_FALSE(t.is_shell());

    NGP_EXPECT_EQ(t.rank(),stk::topology::ELEMENT_RANK);
    NGP_EXPECT_EQ(t.side_rank(),stk::topology::NODE_RANK);
    NGP_EXPECT_EQ(t.num_sides(),2u);

    NGP_EXPECT_EQ(t.dimension(),2u);
    NGP_EXPECT_EQ(t.num_nodes(),2u);
    NGP_EXPECT_EQ(t.num_vertices(),2u);
    NGP_EXPECT_EQ(t.num_edges(),0u);
    NGP_EXPECT_EQ(t.num_faces(),0u);
    NGP_EXPECT_EQ(t.num_permutations(),2u);
    NGP_EXPECT_EQ(t.num_positive_permutations(),2u);

    NGP_EXPECT_TRUE(t.defined_on_spatial_dimension(1));
    NGP_EXPECT_TRUE(t.defined_on_spatial_dimension(2));
    NGP_EXPECT_TRUE(t.defined_on_spatial_dimension(3));

    NGP_EXPECT_EQ(t.base(),stk::topology::SPRING_2);

    NGP_EXPECT_EQ(t.face_topology(0), stk::topology::INVALID_TOPOLOGY);

    unsigned gold_side_node_ordinals[2][1] = { {0},
                                               {1} };
    check_side_node_ordinals_ngp(t, gold_side_node_ordinals);
    check_side_nodes_ngp(t, gold_side_node_ordinals);

    unsigned gold_permutation_node_ordinals[2][2] = { {0, 1},
                                                      {1, 0} };
    check_permutation_node_ordinals_ngp(t, gold_permutation_node_ordinals);
    check_permutation_nodes_ngp(t, gold_permutation_node_ordinals);

    check_equivalent_ngp(t, gold_permutation_node_ordinals);
    check_lexicographical_smallest_permutation_ngp(t, gold_permutation_node_ordinals);
  });
}

NGP_TEST(stk_topology_ngp, spring_2)
{
  check_spring2_on_device();
}

TEST(stk_topology, spring_3)
{
  stk::topology t = stk::topology::SPRING_3;

  EXPECT_TRUE(t.is_valid());
  EXPECT_FALSE(t.has_homogeneous_faces());
  EXPECT_FALSE(t.is_shell());

  EXPECT_EQ(t.rank(),stk::topology::ELEMENT_RANK);
  EXPECT_EQ(t.side_rank(),stk::topology::NODE_RANK);
  EXPECT_EQ(t.num_sides(),2u);

  EXPECT_EQ(t.dimension(),2u);
  EXPECT_EQ(t.num_nodes(),3u);
  EXPECT_EQ(t.num_vertices(),2u);
  EXPECT_EQ(t.num_edges(),0u);
  EXPECT_EQ(t.num_faces(),0u);
  EXPECT_EQ(t.num_permutations(),2u);
  EXPECT_EQ(t.num_positive_permutations(),2u);

  EXPECT_TRUE(t.defined_on_spatial_dimension(1));
  EXPECT_TRUE(t.defined_on_spatial_dimension(2));
  EXPECT_TRUE(t.defined_on_spatial_dimension(3));

  EXPECT_EQ(t.base(),stk::topology::SPRING_2);

  EXPECT_EQ(t.face_topology(0), stk::topology::INVALID_TOPOLOGY);

  std::vector<std::vector<unsigned>> gold_side_node_ordinals = { {0},
                                                                 {1} };
  check_side_node_ordinals(t, gold_side_node_ordinals);
  check_side_nodes(t, gold_side_node_ordinals);

  std::vector<std::vector<unsigned>> gold_permutation_node_ordinals = { {0, 1, 2},
                                                                        {1, 0, 2} };

  check_permutation_node_ordinals(t, gold_permutation_node_ordinals);
  check_permutation_nodes(t, gold_permutation_node_ordinals);

  check_equivalent(t, gold_permutation_node_ordinals);
  check_lexicographical_smallest_permutation(t, gold_permutation_node_ordinals);
}

void check_spring3_on_device()
{
  Kokkos::parallel_for(1, KOKKOS_LAMBDA(const int i)
  {
    stk::topology t = stk::topology::SPRING_3;

    NGP_EXPECT_TRUE(t.is_valid());
    NGP_EXPECT_FALSE(t.has_homogeneous_faces());
    NGP_EXPECT_FALSE(t.is_shell());

    NGP_EXPECT_EQ(t.rank(),stk::topology::ELEMENT_RANK);
    NGP_EXPECT_EQ(t.side_rank(),stk::topology::NODE_RANK);
    NGP_EXPECT_EQ(t.num_sides(),2u);

    NGP_EXPECT_EQ(t.dimension(),2u);
    NGP_EXPECT_EQ(t.num_nodes(),3u);
    NGP_EXPECT_EQ(t.num_vertices(),2u);
    NGP_EXPECT_EQ(t.num_edges(),0u);
    NGP_EXPECT_EQ(t.num_faces(),0u);
    NGP_EXPECT_EQ(t.num_permutations(),2u);
    NGP_EXPECT_EQ(t.num_positive_permutations(),2u);

    NGP_EXPECT_TRUE(t.defined_on_spatial_dimension(1));
    NGP_EXPECT_TRUE(t.defined_on_spatial_dimension(2));
    NGP_EXPECT_TRUE(t.defined_on_spatial_dimension(3));

    NGP_EXPECT_EQ(t.base(),stk::topology::SPRING_2);

    NGP_EXPECT_EQ(t.face_topology(0), stk::topology::INVALID_TOPOLOGY);

    unsigned gold_side_node_ordinals[2][1] = { {0},
                                               {1} };
    check_side_node_ordinals_ngp(t, gold_side_node_ordinals);
    check_side_nodes_ngp(t, gold_side_node_ordinals);

    unsigned gold_permutation_node_ordinals[2][3] = { {0, 1, 2},
                                                      {1, 0, 2} };
    check_permutation_node_ordinals_ngp(t, gold_permutation_node_ordinals);
    check_permutation_nodes_ngp(t, gold_permutation_node_ordinals);

    check_equivalent_ngp(t, gold_permutation_node_ordinals);
    check_lexicographical_smallest_permutation_ngp(t, gold_permutation_node_ordinals);
  });
}

NGP_TEST(stk_topology_ngp, spring_3)
{
  check_spring3_on_device();
}
