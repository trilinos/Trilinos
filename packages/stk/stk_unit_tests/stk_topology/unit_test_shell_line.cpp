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

#include "Kokkos_Core.hpp"            // for parallel_for, KOKKOS_LAMBDA
#include "gtest/gtest.h"              // for AssertionResult, Message, TestPartResult, EXPECT_EQ
#include "stk_ngp_test/ngp_test.hpp"  // for NGP_EXPECT_EQ, NGP_EXPECT_FALSE, NGP_EXPECT_TRUE
#include "stk_topology/topology.hpp"  // for topology, topology::SHELL_LINE_2, topology::EDGE_RANK
#include "topology_test_utils.hpp"    // for check_edge_node_ordinals, check_edge_node_ordinals_ngp
#include <vector>                     // for vector

namespace {

std::vector<std::vector<uint8_t>> get_gold_edge_node_ordinals_shell_line2() {
  return std::vector<std::vector<uint8_t>> {
    {0, 1},
    {1, 0}
  };
}

std::vector<std::vector<uint8_t>> get_gold_permutation_node_ordinals_shell_line2() {
  return std::vector<std::vector<uint8_t>> {
    {0, 1},
    {1, 0}
  };
}

TEST(stk_topology, shell_line_2)
{
  stk::topology t = stk::topology::SHELL_LINE_2;

  EXPECT_TRUE(t.is_valid());
  EXPECT_FALSE(t.has_homogeneous_faces());
  EXPECT_TRUE(t.is_shell());

  EXPECT_EQ(t.rank(),stk::topology::ELEMENT_RANK);
  EXPECT_EQ(t.side_rank(),stk::topology::EDGE_RANK);
  EXPECT_EQ(t.num_sides(),2u);

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

  EXPECT_EQ(t.base(),stk::topology::SHELL_LINE_2);

  EXPECT_EQ(t.face_topology(0), stk::topology::INVALID_TOPOLOGY);

  check_side_node_ordinals(t, get_gold_edge_node_ordinals_shell_line2());
  check_edge_node_ordinals(t, get_gold_edge_node_ordinals_shell_line2());
  check_side_nodes(t, get_gold_edge_node_ordinals_shell_line2());
  check_edge_nodes(t, get_gold_edge_node_ordinals_shell_line2());

  check_permutation_node_ordinals(t, get_gold_permutation_node_ordinals_shell_line2());
  check_permutation_nodes(t, get_gold_permutation_node_ordinals_shell_line2());

  check_equivalent(t, get_gold_permutation_node_ordinals_shell_line2());
  check_lexicographical_smallest_permutation(t, get_gold_permutation_node_ordinals_shell_line2());
}

void check_shell_line_2_on_device()
{
  OrdinalType goldEdgeNodeOrdinals = fillGoldOrdinals(get_gold_edge_node_ordinals_shell_line2());
  OrdinalType goldPermutationNodeOrdinals = fillGoldOrdinals(get_gold_permutation_node_ordinals_shell_line2());

  Kokkos::parallel_for(stk::ngp::DeviceRangePolicy(0, 1), KOKKOS_LAMBDA(const int i)
  {
    stk::topology t = stk::topology::SHELL_LINE_2;

    NGP_EXPECT_TRUE(t.is_valid());
    NGP_EXPECT_FALSE(t.has_homogeneous_faces());
    NGP_EXPECT_TRUE(t.is_shell());

    NGP_EXPECT_EQ(t.rank(),stk::topology::ELEMENT_RANK);
    NGP_EXPECT_EQ(t.side_rank(),stk::topology::EDGE_RANK);
    NGP_EXPECT_EQ(t.num_sides(),2u);

    NGP_EXPECT_EQ(t.dimension(),2u);
    NGP_EXPECT_EQ(t.num_nodes(),2u);
    NGP_EXPECT_EQ(t.num_vertices(),2u);
    NGP_EXPECT_EQ(t.num_edges(),2u);
    NGP_EXPECT_EQ(t.num_faces(),0u);
    NGP_EXPECT_EQ(t.num_permutations(),2u);
    NGP_EXPECT_EQ(t.num_positive_permutations(),1u);

    NGP_EXPECT_FALSE(t.defined_on_spatial_dimension(1));
    NGP_EXPECT_TRUE(t.defined_on_spatial_dimension(2));
    NGP_EXPECT_FALSE(t.defined_on_spatial_dimension(3));

    NGP_EXPECT_EQ(t.base(),stk::topology::SHELL_LINE_2);

    constexpr unsigned numNodes = stk::topology_detail::topology_data<stk::topology::SHELL_LINE_2>::num_nodes;

    check_side_node_ordinals_ngp<numNodes>(t, goldEdgeNodeOrdinals);
    check_edge_node_ordinals_ngp<numNodes>(t, goldEdgeNodeOrdinals);
    check_side_nodes_ngp<numNodes>(t, goldEdgeNodeOrdinals);
    check_edge_nodes_ngp<numNodes>(t, goldEdgeNodeOrdinals);

    check_permutation_node_ordinals_ngp<numNodes>(t, goldPermutationNodeOrdinals);
    check_permutation_nodes_ngp<numNodes>(t, goldPermutationNodeOrdinals);

    check_equivalent_ngp<numNodes>(t, goldPermutationNodeOrdinals);
    check_lexicographical_smallest_permutation_ngp<numNodes>(t, goldPermutationNodeOrdinals);
  });
}

NGP_TEST(stk_topology_ngp, shell_line_2)
{
  check_shell_line_2_on_device();
}


std::vector<std::vector<uint8_t>> get_gold_edge_node_ordinals_shell_line3() {
  return std::vector<std::vector<uint8_t>> {
    {0, 1, 2},
    {1, 0, 2}
  };
}

std::vector<std::vector<uint8_t>> get_gold_permutation_node_ordinals_shell_line3() {
  return std::vector<std::vector<uint8_t>> {
    {0, 1, 2},
    {1, 0, 2}
  };
}

TEST(stk_topology, shell_line_3)
{
  stk::topology t = stk::topology::SHELL_LINE_3;

  EXPECT_TRUE(t.is_valid());
  EXPECT_FALSE(t.has_homogeneous_faces());
  EXPECT_TRUE(t.is_shell());

  EXPECT_EQ(t.rank(),stk::topology::ELEMENT_RANK);
  EXPECT_EQ(t.side_rank(),stk::topology::EDGE_RANK);
  EXPECT_EQ(t.num_sides(),2u);

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

  EXPECT_EQ(t.base(),stk::topology::SHELL_LINE_2);

  EXPECT_EQ(t.face_topology(0), stk::topology::INVALID_TOPOLOGY);

  check_side_node_ordinals(t, get_gold_edge_node_ordinals_shell_line3());
  check_edge_node_ordinals(t, get_gold_edge_node_ordinals_shell_line3());
  check_side_nodes(t, get_gold_edge_node_ordinals_shell_line3());
  check_edge_nodes(t, get_gold_edge_node_ordinals_shell_line3());

  check_permutation_node_ordinals(t, get_gold_permutation_node_ordinals_shell_line3());
  check_permutation_nodes(t, get_gold_permutation_node_ordinals_shell_line3());

  check_equivalent(t, get_gold_permutation_node_ordinals_shell_line3());
  check_lexicographical_smallest_permutation(t, get_gold_permutation_node_ordinals_shell_line3());
}

void check_shell_line_3_on_device()
{
  OrdinalType goldEdgeNodeOrdinals = fillGoldOrdinals(get_gold_edge_node_ordinals_shell_line3());
  OrdinalType goldPermutationNodeOrdinals = fillGoldOrdinals(get_gold_permutation_node_ordinals_shell_line3());

  Kokkos::parallel_for(stk::ngp::DeviceRangePolicy(0, 1), KOKKOS_LAMBDA(const int i)
  {
    stk::topology t = stk::topology::SHELL_LINE_3;

    NGP_EXPECT_TRUE(t.is_valid());
    NGP_EXPECT_FALSE(t.has_homogeneous_faces());
    NGP_EXPECT_TRUE(t.is_shell());

    NGP_EXPECT_EQ(t.rank(),stk::topology::ELEMENT_RANK);
    NGP_EXPECT_EQ(t.side_rank(),stk::topology::EDGE_RANK);
    NGP_EXPECT_EQ(t.num_sides(),2u);

    NGP_EXPECT_EQ(t.dimension(),2u);
    NGP_EXPECT_EQ(t.num_nodes(),3u);
    NGP_EXPECT_EQ(t.num_vertices(),2u);
    NGP_EXPECT_EQ(t.num_edges(),2u);
    NGP_EXPECT_EQ(t.num_faces(),0u);
    NGP_EXPECT_EQ(t.num_permutations(),2u);
    NGP_EXPECT_EQ(t.num_positive_permutations(),1u);

    NGP_EXPECT_FALSE(t.defined_on_spatial_dimension(1));
    NGP_EXPECT_TRUE(t.defined_on_spatial_dimension(2));
    NGP_EXPECT_FALSE(t.defined_on_spatial_dimension(3));

    NGP_EXPECT_EQ(t.base(),stk::topology::SHELL_LINE_2);

    constexpr unsigned numNodes = stk::topology_detail::topology_data<stk::topology::SHELL_LINE_3>::num_nodes;

    check_side_node_ordinals_ngp<numNodes>(t, goldEdgeNodeOrdinals);
    check_edge_node_ordinals_ngp<numNodes>(t, goldEdgeNodeOrdinals);
    check_side_nodes_ngp<numNodes>(t, goldEdgeNodeOrdinals);
    check_edge_nodes_ngp<numNodes>(t, goldEdgeNodeOrdinals);

    check_permutation_node_ordinals_ngp<numNodes>(t, goldPermutationNodeOrdinals);
    check_permutation_nodes_ngp<numNodes>(t, goldPermutationNodeOrdinals);

    check_equivalent_ngp<numNodes>(t, goldPermutationNodeOrdinals);
    check_lexicographical_smallest_permutation_ngp<numNodes>(t, goldPermutationNodeOrdinals);
  });
}

NGP_TEST(stk_topology_ngp, shell_line_3)
{
  check_shell_line_3_on_device();
}

}
