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
#include "stk_topology/topology.hpp"  // for topology, topology::QUAD_4, topology::EDGE_RANK
#include "topology_test_utils.hpp"    // for check_edge_node_ordinals, check_edge_node_ordinals_ngp
#include <stdio.h>                    // for printf
#include <iostream>                   // for operator<<, endl, basic_ostream, cout, ostream
#include <vector>                     // for vector

namespace {

std::vector<std::vector<uint8_t>> get_gold_edge_node_ordinals_quad4() {
  return std::vector<std::vector<uint8_t>> {
    {0, 1},
    {1, 2},
    {2, 3},
    {3, 0}
  };
}

std::vector<std::vector<uint8_t>> get_gold_permutation_node_ordinals_quad4() {
  return std::vector<std::vector<uint8_t>> {
    {0, 1, 2, 3},
    {3, 0, 1, 2},
    {2, 3, 0, 1},
    {1, 2, 3, 0},
    {0, 3, 2, 1},
    {3, 2, 1, 0},
    {2, 1, 0, 3},
    {1, 0, 3, 2}
  };
}

TEST( stk_topology, quad_4 )
{
  using stk::topology;

  topology t = topology::QUAD_4;

  EXPECT_TRUE(t.is_valid());
  EXPECT_FALSE(t.has_homogeneous_faces());
  EXPECT_FALSE(t.is_shell());

  EXPECT_EQ(t.rank(),topology::FACE_RANK);
  EXPECT_EQ(t.side_rank(),topology::EDGE_RANK);
  EXPECT_EQ(t.num_sides(),4u);

  EXPECT_EQ(t.dimension(),2u);
  EXPECT_EQ(t.num_nodes(),4u);
  EXPECT_EQ(t.num_vertices(),4u);
  EXPECT_EQ(t.num_edges(),4u);
  EXPECT_EQ(t.num_faces(),0u);
  EXPECT_EQ(t.num_permutations(),8u);
  EXPECT_EQ(t.num_positive_permutations(),4u);

  EXPECT_FALSE(t.defined_on_spatial_dimension(1));
  EXPECT_FALSE(t.defined_on_spatial_dimension(2));
  EXPECT_TRUE(t.defined_on_spatial_dimension(3));

  EXPECT_EQ(t.base(),topology::QUAD_4);

  EXPECT_EQ(t.face_topology(0), stk::topology::INVALID_TOPOLOGY);

  check_side_node_ordinals(t, get_gold_edge_node_ordinals_quad4());
  check_edge_node_ordinals(t, get_gold_edge_node_ordinals_quad4());
  check_side_nodes(t, get_gold_edge_node_ordinals_quad4());
  check_edge_nodes(t, get_gold_edge_node_ordinals_quad4());

  check_permutation_node_ordinals(t, get_gold_permutation_node_ordinals_quad4());
  check_permutation_nodes(t, get_gold_permutation_node_ordinals_quad4());

  check_equivalent(t, get_gold_permutation_node_ordinals_quad4());
  check_lexicographical_smallest_permutation(t, get_gold_permutation_node_ordinals_quad4());
}

void check_quad_4_on_device()
{
  OrdinalType goldEdgeNodeOrdinals = fillGoldOrdinals(get_gold_edge_node_ordinals_quad4());
  OrdinalType goldPermutationNodeOrdinals = fillGoldOrdinals(get_gold_permutation_node_ordinals_quad4());

  stk::topology t = stk::topology::QUAD_4;
  constexpr unsigned numNodes = stk::topology_detail::topology_data<stk::topology::QUAD_4>::num_nodes;
  EXPECT_EQ(4u, numNodes);

  Kokkos::parallel_for(stk::ngp::DeviceRangePolicy(0, 1), KOKKOS_LAMBDA(const int i)
  {
    NGP_EXPECT_TRUE(t.is_valid());
    NGP_EXPECT_FALSE(t.has_homogeneous_faces());
    NGP_EXPECT_FALSE(t.is_shell());

    NGP_EXPECT_EQ(t.rank(),stk::topology::FACE_RANK);
    NGP_EXPECT_EQ(t.side_rank(),stk::topology::EDGE_RANK);
    NGP_EXPECT_EQ(t.num_sides(),4u);

    NGP_EXPECT_EQ(t.dimension(),2u);
    NGP_EXPECT_EQ(t.num_nodes(),4u);
    NGP_EXPECT_EQ(t.num_vertices(),4u);
    NGP_EXPECT_EQ(t.num_edges(),4u);
    NGP_EXPECT_EQ(t.num_faces(),0u);
    NGP_EXPECT_EQ(t.num_permutations(),8u);
    NGP_EXPECT_EQ(t.num_positive_permutations(),4u);

    NGP_EXPECT_FALSE(t.defined_on_spatial_dimension(1));
    NGP_EXPECT_FALSE(t.defined_on_spatial_dimension(2));
    NGP_EXPECT_TRUE(t.defined_on_spatial_dimension(3));

    NGP_EXPECT_EQ(t.base(),stk::topology::QUAD_4);

    NGP_EXPECT_EQ(t.face_topology(0), stk::topology::INVALID_TOPOLOGY);

    check_side_node_ordinals_ngp<numNodes>(t, goldEdgeNodeOrdinals);
    check_edge_node_ordinals_ngp<numNodes>(t, goldEdgeNodeOrdinals);
    check_side_nodes_ngp<numNodes>(t, goldEdgeNodeOrdinals);
    check_edge_nodes_ngp<numNodes>(t, goldEdgeNodeOrdinals);
  });

  Kokkos::parallel_for(stk::ngp::DeviceRangePolicy(0, 1), KOKKOS_LAMBDA(const int i)
  {
    check_permutation_node_ordinals_ngp<numNodes>(t, goldPermutationNodeOrdinals);
    check_permutation_nodes_ngp<numNodes>(t, goldPermutationNodeOrdinals);

    check_equivalent_ngp<numNodes>(t, goldPermutationNodeOrdinals);
    check_lexicographical_smallest_permutation_ngp<numNodes>(t, goldPermutationNodeOrdinals);
  });
}

NGP_TEST(stk_topology_ngp, quad_4)
{
  check_quad_4_on_device();
}


std::vector<std::vector<uint8_t>> get_gold_edge_node_ordinals_quad6() {
  return std::vector<std::vector<uint8_t>> {
    {0, 1, 4},
    {1, 2},
    {2, 3, 5},
    {3, 0}
  };
}

std::vector<std::vector<uint8_t>> get_gold_permutation_node_ordinals_quad6() {
  return std::vector<std::vector<uint8_t>> {
    {0, 1, 2, 3,  4, 5},
    {3, 0, 1, 2,  4, 5},
    {2, 3, 0, 1,  5, 4},
    {1, 2, 3, 0,  5, 4},
    {0, 3, 2, 1,  5, 4},
    {3, 2, 1, 0,  5, 4},
    {2, 1, 0, 3,  4, 5},
    {1, 0, 3, 2,  4, 5}
  };
}

TEST( stk_topology, quad_6 )
{
  using stk::topology;

  topology t = topology::QUAD_6;

  EXPECT_TRUE(t.is_valid());
  EXPECT_FALSE(t.has_homogeneous_faces());
  EXPECT_FALSE(t.is_shell());

  EXPECT_EQ(t.rank(),topology::FACE_RANK);
  EXPECT_EQ(t.side_rank(),topology::EDGE_RANK);
  EXPECT_EQ(t.num_sides(),4u);

  EXPECT_EQ(t.dimension(),2u);
  EXPECT_EQ(t.num_nodes(),6u);
  EXPECT_EQ(t.num_vertices(),4u);
  EXPECT_EQ(t.num_edges(),4u);
  EXPECT_EQ(t.num_faces(),0u);
  EXPECT_EQ(t.num_permutations(),8u);
  EXPECT_EQ(t.num_positive_permutations(),4u);

  EXPECT_FALSE(t.defined_on_spatial_dimension(1));
  EXPECT_FALSE(t.defined_on_spatial_dimension(2));
  EXPECT_TRUE(t.defined_on_spatial_dimension(3));

  EXPECT_EQ(t.base(),topology::QUAD_4);

  EXPECT_EQ(t.face_topology(0), stk::topology::INVALID_TOPOLOGY);

  check_side_node_ordinals(t, get_gold_edge_node_ordinals_quad6());
  check_edge_node_ordinals(t, get_gold_edge_node_ordinals_quad6());
  check_side_nodes(t, get_gold_edge_node_ordinals_quad6());
  check_edge_nodes(t, get_gold_edge_node_ordinals_quad6());

  std::cout << "Reminder: we still need to enable permutation for QUAD_6" << std::endl;
  const bool enabledPermutation = false;
  if (enabledPermutation) {
    check_permutation_node_ordinals(t, get_gold_permutation_node_ordinals_quad6());
  }
  check_permutation_nodes(t, get_gold_permutation_node_ordinals_quad6());

  check_equivalent(t, get_gold_permutation_node_ordinals_quad6());
  check_lexicographical_smallest_permutation(t, get_gold_permutation_node_ordinals_quad6());
}

void check_quad_6_on_device()
{
  OrdinalType goldEdgeNodeOrdinals = fillGoldOrdinals(get_gold_edge_node_ordinals_quad6());
  OrdinalType goldPermutationNodeOrdinals = fillGoldOrdinals(get_gold_permutation_node_ordinals_quad6());

  stk::topology t = stk::topology::QUAD_6;
  constexpr unsigned numNodes = stk::topology_detail::topology_data<stk::topology::QUAD_6>::num_nodes;
  EXPECT_EQ(6u, numNodes);

  Kokkos::parallel_for(stk::ngp::DeviceRangePolicy(0, 1), KOKKOS_LAMBDA(const int i)
  {
    NGP_EXPECT_TRUE(t.is_valid());
    NGP_EXPECT_FALSE(t.has_homogeneous_faces());
    NGP_EXPECT_FALSE(t.is_shell());

    NGP_EXPECT_EQ(t.rank(),stk::topology::FACE_RANK);
    NGP_EXPECT_EQ(t.side_rank(),stk::topology::EDGE_RANK);
    NGP_EXPECT_EQ(t.num_sides(),4u);

    NGP_EXPECT_EQ(t.dimension(),2u);
    NGP_EXPECT_EQ(t.num_nodes(),6u);
    NGP_EXPECT_EQ(t.num_vertices(),4u);
    NGP_EXPECT_EQ(t.num_edges(),4u);
    NGP_EXPECT_EQ(t.num_faces(),0u);
    NGP_EXPECT_EQ(t.num_permutations(),8u);
    NGP_EXPECT_EQ(t.num_positive_permutations(),4u);

    NGP_EXPECT_FALSE(t.defined_on_spatial_dimension(1));
    NGP_EXPECT_FALSE(t.defined_on_spatial_dimension(2));
    NGP_EXPECT_TRUE(t.defined_on_spatial_dimension(3));

    NGP_EXPECT_EQ(t.base(),stk::topology::QUAD_4);

    NGP_EXPECT_EQ(t.face_topology(0), stk::topology::INVALID_TOPOLOGY);

    check_side_node_ordinals_ngp<numNodes>(t, goldEdgeNodeOrdinals);
    check_edge_node_ordinals_ngp<numNodes>(t, goldEdgeNodeOrdinals);
    check_side_nodes_ngp<numNodes>(t, goldEdgeNodeOrdinals);
    check_edge_nodes_ngp<numNodes>(t, goldEdgeNodeOrdinals);
  });

  Kokkos::parallel_for(stk::ngp::DeviceRangePolicy(0, 1), KOKKOS_LAMBDA(const int i)
  {
#if KOKKOS_VERSION < 40200
    printf("Reminder: we still need to enable permutation for QUAD_6\n");
#else
    Kokkos::printf("Reminder: we still need to enable permutation for QUAD_6\n");
#endif
    const bool enabledPermutation = false;
    if (enabledPermutation) {
      check_permutation_node_ordinals_ngp<numNodes>(t, goldPermutationNodeOrdinals);
    }
    check_permutation_nodes_ngp<numNodes>(t, goldPermutationNodeOrdinals);

    check_equivalent_ngp<numNodes>(t, goldPermutationNodeOrdinals);
    check_lexicographical_smallest_permutation_ngp<numNodes>(t, goldPermutationNodeOrdinals);
  });
}

NGP_TEST(stk_topology_ngp, quad_6)
{
  check_quad_6_on_device();
}


std::vector<std::vector<uint8_t>> get_gold_edge_node_ordinals_quad8() {
  return std::vector<std::vector<uint8_t>> {
    {0, 1, 4},
    {1, 2, 5},
    {2, 3, 6},
    {3, 0, 7}
  };
}

std::vector<std::vector<uint8_t>> get_gold_permutation_node_ordinals_quad8() {
  return std::vector<std::vector<uint8_t>> {
    {0, 1, 2, 3,  4, 5, 6, 7},
    {3, 0, 1, 2,  7, 4, 5, 6},
    {2, 3, 0, 1,  6, 7, 4, 5},
    {1, 2, 3, 0,  5, 6, 7, 4},
    {0, 3, 2, 1,  7, 6, 5, 4},
    {3, 2, 1, 0,  6, 5, 4, 7},
    {2, 1, 0, 3,  5, 4, 7, 6},
    {1, 0, 3, 2,  4, 7, 6, 5}
  };
}

TEST( stk_topology, quad_8 )
{
  using stk::topology;

  topology t = topology::QUAD_8;

  EXPECT_TRUE(t.is_valid());
  EXPECT_FALSE(t.has_homogeneous_faces());
  EXPECT_FALSE(t.is_shell());

  EXPECT_EQ(t.rank(),topology::FACE_RANK);
  EXPECT_EQ(t.side_rank(),topology::EDGE_RANK);
  EXPECT_EQ(t.num_sides(),4u);

  EXPECT_EQ(t.dimension(),2u);
  EXPECT_EQ(t.num_nodes(),8u);
  EXPECT_EQ(t.num_vertices(),4u);
  EXPECT_EQ(t.num_edges(),4u);
  EXPECT_EQ(t.num_faces(),0u);
  EXPECT_EQ(t.num_permutations(),8u);
  EXPECT_EQ(t.num_positive_permutations(),4u);

  EXPECT_FALSE(t.defined_on_spatial_dimension(1));
  EXPECT_FALSE(t.defined_on_spatial_dimension(2));
  EXPECT_TRUE(t.defined_on_spatial_dimension(3));

  EXPECT_EQ(t.base(),topology::QUAD_4);

  EXPECT_EQ(t.face_topology(0), stk::topology::INVALID_TOPOLOGY);

  check_side_node_ordinals(t, get_gold_edge_node_ordinals_quad8());
  check_edge_node_ordinals(t, get_gold_edge_node_ordinals_quad8());
  check_side_nodes(t, get_gold_edge_node_ordinals_quad8());
  check_edge_nodes(t, get_gold_edge_node_ordinals_quad8());

  check_permutation_node_ordinals(t, get_gold_permutation_node_ordinals_quad8());
  check_permutation_nodes(t, get_gold_permutation_node_ordinals_quad8());

  check_equivalent(t, get_gold_permutation_node_ordinals_quad8());
  check_lexicographical_smallest_permutation(t, get_gold_permutation_node_ordinals_quad8());
}

void check_quad_8_on_device()
{
  OrdinalType goldEdgeNodeOrdinals = fillGoldOrdinals(get_gold_edge_node_ordinals_quad8());
  OrdinalType goldPermutationNodeOrdinals = fillGoldOrdinals(get_gold_permutation_node_ordinals_quad8());

  stk::topology t = stk::topology::QUAD_8;
  constexpr unsigned numNodes = stk::topology_detail::topology_data<stk::topology::QUAD_8>::num_nodes;
  EXPECT_EQ(8u, numNodes);

  Kokkos::parallel_for(stk::ngp::DeviceRangePolicy(0, 1), KOKKOS_LAMBDA(const int i)
  {
    NGP_EXPECT_TRUE(t.is_valid());
    NGP_EXPECT_FALSE(t.has_homogeneous_faces());
    NGP_EXPECT_FALSE(t.is_shell());

    NGP_EXPECT_EQ(t.rank(),stk::topology::FACE_RANK);
    NGP_EXPECT_EQ(t.side_rank(),stk::topology::EDGE_RANK);
    NGP_EXPECT_EQ(t.num_sides(),4u);

    NGP_EXPECT_EQ(t.dimension(),2u);
    NGP_EXPECT_EQ(t.num_nodes(),8u);
    NGP_EXPECT_EQ(t.num_vertices(),4u);
    NGP_EXPECT_EQ(t.num_edges(),4u);
    NGP_EXPECT_EQ(t.num_faces(),0u);
    NGP_EXPECT_EQ(t.num_permutations(),8u);
    NGP_EXPECT_EQ(t.num_positive_permutations(),4u);

    NGP_EXPECT_FALSE(t.defined_on_spatial_dimension(1));
    NGP_EXPECT_FALSE(t.defined_on_spatial_dimension(2));
    NGP_EXPECT_TRUE(t.defined_on_spatial_dimension(3));

    NGP_EXPECT_EQ(t.base(),stk::topology::QUAD_4);

    NGP_EXPECT_EQ(t.face_topology(0), stk::topology::INVALID_TOPOLOGY);

    check_side_node_ordinals_ngp<numNodes>(t, goldEdgeNodeOrdinals);
    check_edge_node_ordinals_ngp<numNodes>(t, goldEdgeNodeOrdinals);
    check_side_nodes_ngp<numNodes>(t, goldEdgeNodeOrdinals);
    check_edge_nodes_ngp<numNodes>(t, goldEdgeNodeOrdinals);
  });

  Kokkos::parallel_for(stk::ngp::DeviceRangePolicy(0, 1), KOKKOS_LAMBDA(const int i)
  {
    check_permutation_node_ordinals_ngp<numNodes>(t, goldPermutationNodeOrdinals);
    check_permutation_nodes_ngp<numNodes>(t, goldPermutationNodeOrdinals);

    check_equivalent_ngp<numNodes>(t, goldPermutationNodeOrdinals);
    check_lexicographical_smallest_permutation_ngp<numNodes>(t, goldPermutationNodeOrdinals);
  });
}

NGP_TEST(stk_topology_ngp, quad_8)
{
  check_quad_8_on_device();
}


std::vector<std::vector<uint8_t>> get_gold_edge_node_ordinals_quad9() {
  return std::vector<std::vector<uint8_t>> {
    {0, 1, 4},
    {1, 2, 5},
    {2, 3, 6},
    {3, 0, 7}
  };
}

std::vector<std::vector<uint8_t>> get_gold_permutation_node_ordinals_quad9() {
  return std::vector<std::vector<uint8_t>> {
    {0, 1, 2, 3,  4, 5, 6, 7,  8},
    {3, 0, 1, 2,  7, 4, 5, 6,  8},
    {2, 3, 0, 1,  6, 7, 4, 5,  8},
    {1, 2, 3, 0,  5, 6, 7, 4,  8},
    {0, 3, 2, 1,  7, 6, 5, 4,  8},
    {3, 2, 1, 0,  6, 5, 4, 7,  8},
    {2, 1, 0, 3,  5, 4, 7, 6,  8},
    {1, 0, 3, 2,  4, 7, 6, 5,  8}
  };
}

TEST( stk_topology, quad_9 )
{
  using stk::topology;

  topology t = topology::QUAD_9;

  EXPECT_TRUE(t.is_valid());
  EXPECT_FALSE(t.has_homogeneous_faces());
  EXPECT_FALSE(t.is_shell());

  EXPECT_EQ(t.rank(),topology::FACE_RANK);
  EXPECT_EQ(t.side_rank(),topology::EDGE_RANK);
  EXPECT_EQ(t.num_sides(),4u);

  EXPECT_EQ(t.dimension(),2u);
  EXPECT_EQ(t.num_nodes(),9u);
  EXPECT_EQ(t.num_vertices(),4u);
  EXPECT_EQ(t.num_edges(),4u);
  EXPECT_EQ(t.num_faces(),0u);
  EXPECT_EQ(t.num_permutations(),8u);
  EXPECT_EQ(t.num_positive_permutations(),4u);

  EXPECT_FALSE(t.defined_on_spatial_dimension(1));
  EXPECT_FALSE(t.defined_on_spatial_dimension(2));
  EXPECT_TRUE(t.defined_on_spatial_dimension(3));

  EXPECT_EQ(t.base(),topology::QUAD_4);

  EXPECT_EQ(t.face_topology(0), stk::topology::INVALID_TOPOLOGY);

  check_side_node_ordinals(t, get_gold_edge_node_ordinals_quad9());
  check_edge_node_ordinals(t, get_gold_edge_node_ordinals_quad9());
  check_side_nodes(t, get_gold_edge_node_ordinals_quad9());
  check_edge_nodes(t, get_gold_edge_node_ordinals_quad9());

  check_permutation_node_ordinals(t, get_gold_permutation_node_ordinals_quad9());
  check_permutation_nodes(t, get_gold_permutation_node_ordinals_quad9());

  check_equivalent(t, get_gold_permutation_node_ordinals_quad9());
  check_lexicographical_smallest_permutation(t, get_gold_permutation_node_ordinals_quad9());
}

void check_quad_9_on_device()
{
  OrdinalType goldEdgeNodeOrdinals = fillGoldOrdinals(get_gold_edge_node_ordinals_quad9());
  OrdinalType goldPermutationNodeOrdinals = fillGoldOrdinals(get_gold_permutation_node_ordinals_quad9());

  stk::topology t = stk::topology::QUAD_9;
  constexpr unsigned numNodes = stk::topology_detail::topology_data<stk::topology::QUAD_9>::num_nodes;
  EXPECT_EQ(9u, numNodes);
 
  Kokkos::parallel_for(stk::ngp::DeviceRangePolicy(0, 1), KOKKOS_LAMBDA(const int i)
  {
    NGP_EXPECT_TRUE(t.is_valid());
    NGP_EXPECT_FALSE(t.has_homogeneous_faces());
    NGP_EXPECT_FALSE(t.is_shell());

    NGP_EXPECT_EQ(t.rank(),stk::topology::FACE_RANK);
    NGP_EXPECT_EQ(t.side_rank(),stk::topology::EDGE_RANK);
    NGP_EXPECT_EQ(t.num_sides(),4u);

    NGP_EXPECT_EQ(t.dimension(),2u);
    NGP_EXPECT_EQ(t.num_nodes(),9u);
    NGP_EXPECT_EQ(t.num_vertices(),4u);
    NGP_EXPECT_EQ(t.num_edges(),4u);
    NGP_EXPECT_EQ(t.num_faces(),0u);
    NGP_EXPECT_EQ(t.num_permutations(),8u);
    NGP_EXPECT_EQ(t.num_positive_permutations(),4u);

    NGP_EXPECT_FALSE(t.defined_on_spatial_dimension(1));
    NGP_EXPECT_FALSE(t.defined_on_spatial_dimension(2));
    NGP_EXPECT_TRUE(t.defined_on_spatial_dimension(3));

    NGP_EXPECT_EQ(t.base(),stk::topology::QUAD_4);

    NGP_EXPECT_EQ(t.face_topology(0), stk::topology::INVALID_TOPOLOGY);

    check_side_node_ordinals_ngp<numNodes>(t, goldEdgeNodeOrdinals);
    check_edge_node_ordinals_ngp<numNodes>(t, goldEdgeNodeOrdinals);
    check_side_nodes_ngp<numNodes>(t, goldEdgeNodeOrdinals);
    check_edge_nodes_ngp<numNodes>(t, goldEdgeNodeOrdinals);
  });

  Kokkos::parallel_for(stk::ngp::DeviceRangePolicy(0, 1), KOKKOS_LAMBDA(const int i)
  {
    check_permutation_node_ordinals_ngp<numNodes>(t, goldPermutationNodeOrdinals);
    check_permutation_nodes_ngp<numNodes>(t, goldPermutationNodeOrdinals);

    check_equivalent_ngp<numNodes>(t, goldPermutationNodeOrdinals);
    check_lexicographical_smallest_permutation_ngp<numNodes>(t, goldPermutationNodeOrdinals);
  });
}

NGP_TEST(stk_topology_ngp, quad_9)
{
  check_quad_9_on_device();
}

}
