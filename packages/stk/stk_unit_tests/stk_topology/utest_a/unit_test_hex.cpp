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

#include "Kokkos_Core.hpp"                   // for parallel_for, KOKKOS_LAMBDA
#include "gtest/gtest.h"                     // for AssertionResult, Message, TestPartResult
#include "stk_ngp_test/ngp_test.hpp"         // for NGP_EXPECT_EQ, NGP_EXPECT_FALSE, NGP_EXPECT_...
#include "stk_topology/topology.hpp"         // for topology, topology::QUAD_4, topology::QUAD_8
#include "stk_util/environment/CPUTime.hpp"  // for cpu_time
#include "topology_test_utils.hpp"           // for check_edge_node_ordinals, check_edge_node_or...
#include <cstddef>                           // for size_t
#include <iostream>                          // for operator<<, basic_ostream, basic_ostream<>::...
#include <vector>                            // for vector

namespace {

TEST(stk_topology, DISABLED_hex_8_defined_on_spatial_dimension_perf)
{
  stk::topology hex8 = stk::topology::HEX_8;

  const size_t numIterations = 10000000000;
  size_t result = 0;

  double startTime = stk::cpu_time();
  for (size_t i = 0; i < numIterations; ++i) {
    result += hex8.defined_on_spatial_dimension(i % 4) ? 1 : 0;
  }
  double stopTime = stk::cpu_time();
  std::cout << "result = " << result << ", time = " << stopTime - startTime << "s" << std::endl;
}

TEST(stk_topology, DISABLED_hex_8_side_node_ordinals_perf)
{
  stk::topology hex8 = stk::topology::HEX_8;

  const size_t numIterations = 1000000000;
  size_t result = 0;

  double startTime = stk::cpu_time();
  std::vector<uint8_t> sideNodeOrdinals(4);
  const unsigned numSides = hex8.num_sides();
  for (size_t i = 0; i < numIterations; ++i) {
    for (unsigned side = 0; side < numSides; ++side) {
      hex8.side_node_ordinals(side, sideNodeOrdinals.data());
      result += sideNodeOrdinals[0];
    }
  }
  double stopTime = stk::cpu_time();
  std::cout << "result = " << result << ", time = " << stopTime - startTime << "s" << std::endl;
}

TEST(stk_topology, DISABLED_hex_8_is_equivalent_perf)
{
  stk::topology hex8 = stk::topology::HEX_8;

  const size_t numIterations = 100000000;
  size_t testResult = 0;
  const unsigned numNodes = hex8.num_nodes();
  const unsigned numPermutations = hex8.num_permutations();

  double startTime = stk::cpu_time();
  std::vector<unsigned> baseNodeArray(numNodes);
  std::vector<unsigned> permutedNodeArray(numNodes);
  hex8.permutation_node_ordinals(0, baseNodeArray.data());
  hex8.permutation_node_ordinals(numPermutations-1, permutedNodeArray.data());

  for (size_t i = 0; i < numIterations; ++i) {
    stk::EquivalentPermutation result = hex8.is_equivalent(baseNodeArray.data(), permutedNodeArray.data());
    testResult += result.permutation_number;
  }
  double stopTime = stk::cpu_time();
  std::cout << "result = " << testResult << ", time = " << stopTime - startTime << "s" << std::endl;
}

std::vector<std::vector<uint8_t>> get_gold_edge_node_ordinals_hex8() {
  return std::vector<std::vector<uint8_t>> {
    {0, 1},
    {1, 2},
    {2, 3},
    {3, 0},
    {4, 5},
    {5, 6},
    {6, 7},
    {7, 4},
    {0, 4},
    {1, 5},
    {2, 6},
    {3, 7}
  };
}

std::vector<std::vector<uint8_t>> get_gold_face_node_ordinals_hex8() {
  return std::vector<std::vector<uint8_t>> {
    {0, 1, 5, 4},
    {1, 2, 6, 5},
    {2, 3, 7, 6},
    {0, 4, 7, 3},
    {0, 3, 2, 1},
    {4, 5, 6, 7}
  };
}

std::vector<std::vector<uint8_t>> get_gold_permutation_node_ordinals_hex8() {
  return std::vector<std::vector<uint8_t>> {
    {0, 1, 2, 3, 4, 5, 6, 7},
    {0, 1, 5, 4, 3, 2, 6, 7},
    {0, 4, 7, 3, 1, 5, 6, 2},
    {1, 2, 3, 0, 5, 6, 7, 4},
    {1, 2, 6, 5, 0, 3, 7, 4},
    {1, 5, 4, 0, 2, 6, 7, 3},
    {2, 3, 0, 1, 6, 7, 4, 5},
    {2, 3, 7, 6, 1, 0, 4, 5},
    {2, 6, 5, 1, 3, 7, 4, 0},
    {3, 0, 1, 2, 7, 4, 5, 6},
    {3, 0, 4, 7, 2, 1, 5, 6},
    {3, 7, 6, 2, 0, 4, 5, 1},
    {4, 0, 1, 5, 7, 3, 2, 6},
    {4, 7, 3, 0, 5, 6, 2, 1},
    {4, 7, 6, 5, 0, 3, 2, 1},
    {5, 1, 2, 6, 4, 0, 3, 7},
    {5, 4, 0, 1, 6, 7, 3, 2},
    {5, 4, 7, 6, 1, 0, 3, 2},
    {6, 2, 3, 7, 5, 1, 0, 4},
    {6, 5, 1, 2, 7, 4, 0, 3},
    {6, 5, 4, 7, 2, 1, 0, 3},
    {7, 3, 0, 4, 6, 2, 1, 5},
    {7, 6, 2, 3, 4, 5, 1, 0},
    {7, 6, 5, 4, 3, 2, 1, 0}
  };
}

TEST(stk_topology, hex_8)
{
  stk::topology t = stk::topology::HEX_8;

  EXPECT_TRUE(t.is_valid());
  EXPECT_TRUE(t.has_homogeneous_faces());
  EXPECT_FALSE(t.is_shell());

  EXPECT_EQ(t.rank(),stk::topology::ELEMENT_RANK);
  EXPECT_EQ(t.side_rank(),stk::topology::FACE_RANK);
  EXPECT_EQ(t.num_sides(),6u);

  EXPECT_EQ(t.num_nodes(),8u);
  EXPECT_EQ(t.num_vertices(),8u);
  EXPECT_EQ(t.num_edges(),12u);
  EXPECT_EQ(t.num_faces(),6u);

  EXPECT_FALSE(t.defined_on_spatial_dimension(1));
  EXPECT_FALSE(t.defined_on_spatial_dimension(2));
  EXPECT_TRUE(t.defined_on_spatial_dimension(3));

  EXPECT_EQ(t.base(),stk::topology::HEX_8);

  EXPECT_EQ(t.face_topology(0), stk::topology::QUAD_4);
  EXPECT_EQ(t.face_topology(1), stk::topology::QUAD_4);
  EXPECT_EQ(t.face_topology(2), stk::topology::QUAD_4);
  EXPECT_EQ(t.face_topology(3), stk::topology::QUAD_4);
  EXPECT_EQ(t.face_topology(4), stk::topology::QUAD_4);
  EXPECT_EQ(t.face_topology(5), stk::topology::QUAD_4);

  check_edge_node_ordinals(t, get_gold_edge_node_ordinals_hex8());
  check_edge_nodes(t, get_gold_edge_node_ordinals_hex8());

  check_side_node_ordinals(t, get_gold_face_node_ordinals_hex8());
  check_face_node_ordinals(t, get_gold_face_node_ordinals_hex8());
  check_side_nodes(t, get_gold_face_node_ordinals_hex8());
  check_face_nodes(t, get_gold_face_node_ordinals_hex8());

  check_permutation_node_ordinals(t, get_gold_permutation_node_ordinals_hex8());
  check_permutation_nodes(t, get_gold_permutation_node_ordinals_hex8());

  check_equivalent(t, get_gold_permutation_node_ordinals_hex8());
  check_lexicographical_smallest_permutation(t, get_gold_permutation_node_ordinals_hex8());
}

void check_hex_8_on_device()
{
  OrdinalType goldEdgeNodeOrdinals = fillGoldOrdinals(get_gold_edge_node_ordinals_hex8());
  OrdinalType goldFaceNodeOrdinals = fillGoldOrdinals(get_gold_face_node_ordinals_hex8());
  OrdinalType goldPermutationNodeOrdinals = fillGoldOrdinals(get_gold_permutation_node_ordinals_hex8());

  stk::topology t = stk::topology::HEX_8;
  constexpr unsigned numNodes = stk::topology_detail::topology_data<stk::topology::HEX_8>::num_nodes;
  EXPECT_EQ(8u, numNodes);

  Kokkos::parallel_for(stk::ngp::DeviceRangePolicy(0, 1), KOKKOS_LAMBDA(const int i)
  {
    NGP_EXPECT_TRUE(t.is_valid());
    NGP_EXPECT_TRUE(t.has_homogeneous_faces());
    NGP_EXPECT_FALSE(t.is_shell());

    NGP_EXPECT_EQ(t.rank(),stk::topology::ELEMENT_RANK);
    NGP_EXPECT_EQ(t.side_rank(),stk::topology::FACE_RANK);
    NGP_EXPECT_EQ(t.num_sides(),6u);

    NGP_EXPECT_EQ(t.num_nodes(),8u);
    NGP_EXPECT_EQ(t.num_vertices(),8u);
    NGP_EXPECT_EQ(t.num_edges(),12u);
    NGP_EXPECT_EQ(t.num_faces(),6u);

    NGP_EXPECT_FALSE(t.defined_on_spatial_dimension(1));
    NGP_EXPECT_FALSE(t.defined_on_spatial_dimension(2));
    NGP_EXPECT_TRUE(t.defined_on_spatial_dimension(3));

    NGP_EXPECT_EQ(t.base(),stk::topology::HEX_8);

    NGP_EXPECT_EQ(t.face_topology(0), stk::topology::QUAD_4);
    NGP_EXPECT_EQ(t.face_topology(1), stk::topology::QUAD_4);
    NGP_EXPECT_EQ(t.face_topology(2), stk::topology::QUAD_4);
    NGP_EXPECT_EQ(t.face_topology(3), stk::topology::QUAD_4);
    NGP_EXPECT_EQ(t.face_topology(4), stk::topology::QUAD_4);
    NGP_EXPECT_EQ(t.face_topology(5), stk::topology::QUAD_4);

    check_edge_node_ordinals_ngp<numNodes>(t, goldEdgeNodeOrdinals);
    check_edge_nodes_ngp<numNodes>(t, goldEdgeNodeOrdinals);

    check_side_node_ordinals_ngp<numNodes>(t, goldFaceNodeOrdinals);
    check_face_node_ordinals_ngp<numNodes>(t, goldFaceNodeOrdinals);
    check_side_nodes_ngp<numNodes>(t, goldFaceNodeOrdinals);
    check_face_nodes_ngp<numNodes>(t, goldFaceNodeOrdinals);
  });

  Kokkos::parallel_for(stk::ngp::DeviceRangePolicy(0, 1), KOKKOS_LAMBDA(const int i)
  {
    check_permutation_node_ordinals_ngp<numNodes>(t, goldPermutationNodeOrdinals);
    check_permutation_nodes_ngp<numNodes>(t, goldPermutationNodeOrdinals);

    check_equivalent_ngp<numNodes>(t, goldPermutationNodeOrdinals);
    check_lexicographical_smallest_permutation_ngp<numNodes>(t, goldPermutationNodeOrdinals);
  });
}

NGP_TEST(stk_topology_ngp, hex_8)
{
  check_hex_8_on_device();
}


std::vector<std::vector<uint8_t>> get_gold_edge_node_ordinals_hex20() {
  return std::vector<std::vector<uint8_t>> {
    {0, 1, 8},
    {1, 2, 9},
    {2, 3, 10},
    {3, 0, 11},
    {4, 5, 16},
    {5, 6, 17},
    {6, 7, 18},
    {7, 4, 19},
    {0, 4, 12},
    {1, 5, 13},
    {2, 6, 14},
    {3, 7, 15}
  };
}

std::vector<std::vector<uint8_t>> get_gold_face_node_ordinals_hex20() {
  return std::vector<std::vector<uint8_t>> {
    {0, 1, 5, 4,  8, 13, 16, 12},
    {1, 2, 6, 5,  9, 14, 17, 13},
    {2, 3, 7, 6, 10, 15, 18, 14},
    {0, 4, 7, 3, 12, 19, 15, 11},
    {0, 3, 2, 1, 11, 10,  9,  8},
    {4, 5, 6, 7, 16, 17, 18, 19}
  };
}

std::vector<std::vector<uint8_t>> get_gold_permutation_node_ordinals_hex20() {
  return std::vector<std::vector<uint8_t>> {
    {0, 1, 2, 3, 4, 5, 6, 7,  8,  9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19},
    {0, 1, 5, 4, 3, 2, 6, 7,  8, 13, 16, 12, 11,  9, 17, 19, 10, 14, 18, 15},
    {0, 4, 7, 3, 1, 5, 6, 2, 12, 19, 15, 11,  8, 16, 18, 10, 13, 17, 14,  9},
    {1, 2, 3, 0, 5, 6, 7, 4,  9, 10, 11,  8, 13, 14, 15, 12, 17, 18, 19, 16},
    {1, 2, 6, 5, 0, 3, 7, 4,  9, 14, 17, 13,  8, 10, 18, 16, 11, 15, 19, 12},
    {1, 5, 4, 0, 2, 6, 7, 3, 13, 16, 12,  8,  9, 17, 19, 11, 14, 18, 15, 10},
    {2, 3, 0, 1, 6, 7, 4, 5, 10, 11,  8,  9, 14, 15, 12, 13, 18, 19, 16, 17},
    {2, 3, 7, 6, 1, 0, 4, 5, 10, 15, 18, 14,  9, 11, 19, 17,  8, 12, 16, 13},
    {2, 6, 5, 1, 3, 7, 4, 0, 14, 17, 13,  9, 10, 18, 16,  8, 15, 19, 12, 11},
    {3, 0, 1, 2, 7, 4, 5, 6, 11,  8,  9, 10, 15, 12, 13, 14, 19, 16, 17, 18},
    {3, 0, 4, 7, 2, 1, 5, 6, 11, 12, 19, 15, 10,  8, 16, 18,  9, 13, 17, 14},
    {3, 7, 6, 2, 0, 4, 5, 1, 15, 18, 14, 10, 11, 19, 17,  9, 12, 16, 13,  8},
    {4, 0, 1, 5, 7, 3, 2, 6, 12,  8, 13, 16, 19, 11,  9, 17, 15, 10, 14, 18},
    {4, 7, 3, 0, 5, 6, 2, 1, 19, 15, 11, 12, 16, 18, 10,  8, 17, 14,  9, 13},
    {4, 7, 6, 5, 0, 3, 2, 1, 19, 18, 17, 16, 12, 15, 14, 13, 11, 10,  9,  8},
    {5, 1, 2, 6, 4, 0, 3, 7, 13,  9, 14, 17, 16,  8, 10, 18, 12, 11, 15, 19},
    {5, 4, 0, 1, 6, 7, 3, 2, 16, 12,  8, 13, 17, 19, 11,  9, 18, 15, 10, 14},
    {5, 4, 7, 6, 1, 0, 3, 2, 16, 19, 18, 17, 13, 12, 15, 14,  8, 11, 10,  9},
    {6, 2, 3, 7, 5, 1, 0, 4, 14, 10, 15, 18, 17,  9, 11, 19, 13,  8, 12, 16},
    {6, 5, 1, 2, 7, 4, 0, 3, 17, 13,  9, 14, 18, 16,  8, 10, 19, 12, 11, 15},
    {6, 5, 4, 7, 2, 1, 0, 3, 17, 16, 19, 18, 14, 13, 12, 15,  9,  8, 11, 10},
    {7, 3, 0, 4, 6, 2, 1, 5, 15, 11, 12, 19, 18, 10,  8, 16, 14,  9, 13, 17},
    {7, 6, 2, 3, 4, 5, 1, 0, 18, 14, 10, 15, 19, 17,  9, 11, 16, 13,  8, 12},
    {7, 6, 5, 4, 3, 2, 1, 0, 18, 17, 16, 19, 15, 14, 13, 12, 10,  9,  8, 11}
  };
}

TEST(stk_topology, hex_20)
{
  stk::topology t = stk::topology::HEX_20;

  EXPECT_TRUE(t.is_valid());
  EXPECT_TRUE(t.has_homogeneous_faces());
  EXPECT_FALSE(t.is_shell());

  EXPECT_EQ(t.rank(),stk::topology::ELEMENT_RANK);
  EXPECT_EQ(t.side_rank(),stk::topology::FACE_RANK);
  EXPECT_EQ(t.num_sides(),6u);

  EXPECT_EQ(t.num_nodes(),20u);
  EXPECT_EQ(t.num_vertices(),8u);
  EXPECT_EQ(t.num_edges(),12u);
  EXPECT_EQ(t.num_faces(),6u);

  EXPECT_FALSE(t.defined_on_spatial_dimension(1));
  EXPECT_FALSE(t.defined_on_spatial_dimension(2));
  EXPECT_TRUE(t.defined_on_spatial_dimension(3));

  EXPECT_EQ(t.base(),stk::topology::HEX_8);

  EXPECT_EQ(t.face_topology(0), stk::topology::QUAD_8);
  EXPECT_EQ(t.face_topology(1), stk::topology::QUAD_8);
  EXPECT_EQ(t.face_topology(2), stk::topology::QUAD_8);
  EXPECT_EQ(t.face_topology(3), stk::topology::QUAD_8);
  EXPECT_EQ(t.face_topology(4), stk::topology::QUAD_8);
  EXPECT_EQ(t.face_topology(5), stk::topology::QUAD_8);

  check_edge_node_ordinals(t, get_gold_edge_node_ordinals_hex20());
  check_edge_nodes(t, get_gold_edge_node_ordinals_hex20());

  check_side_node_ordinals(t, get_gold_face_node_ordinals_hex20());
  check_face_node_ordinals(t, get_gold_face_node_ordinals_hex20());
  check_side_nodes(t, get_gold_face_node_ordinals_hex20());
  check_face_nodes(t, get_gold_face_node_ordinals_hex20());

  check_permutation_node_ordinals(t, get_gold_permutation_node_ordinals_hex20());
  check_permutation_nodes(t, get_gold_permutation_node_ordinals_hex20());

  check_equivalent(t, get_gold_permutation_node_ordinals_hex20());
  check_lexicographical_smallest_permutation(t, get_gold_permutation_node_ordinals_hex20());
}

void check_hex_20_on_device()
{
  OrdinalType goldEdgeNodeOrdinals = fillGoldOrdinals(get_gold_edge_node_ordinals_hex20());
  OrdinalType goldFaceNodeOrdinals = fillGoldOrdinals(get_gold_face_node_ordinals_hex20());
  OrdinalType goldPermutationNodeOrdinals = fillGoldOrdinals(get_gold_permutation_node_ordinals_hex20());

  const stk::topology t = stk::topology::HEX_20;
  constexpr unsigned numNodes = stk::topology_detail::topology_data<stk::topology::HEX_20>::num_nodes;
  EXPECT_EQ(20u, numNodes);

  Kokkos::parallel_for(stk::ngp::DeviceRangePolicy(0, 1), KOKKOS_LAMBDA(const int i)
  {
    NGP_EXPECT_TRUE(t.is_valid());
    NGP_EXPECT_TRUE(t.has_homogeneous_faces());
    NGP_EXPECT_FALSE(t.is_shell());

    NGP_EXPECT_EQ(t.rank(),stk::topology::ELEMENT_RANK);
    NGP_EXPECT_EQ(t.side_rank(),stk::topology::FACE_RANK);
    NGP_EXPECT_EQ(t.num_sides(),6u);

    NGP_EXPECT_EQ(t.num_nodes(),20u);
    NGP_EXPECT_EQ(t.num_vertices(),8u);
    NGP_EXPECT_EQ(t.num_edges(),12u);
    NGP_EXPECT_EQ(t.num_faces(),6u);

    NGP_EXPECT_FALSE(t.defined_on_spatial_dimension(1));
    NGP_EXPECT_FALSE(t.defined_on_spatial_dimension(2));
    NGP_EXPECT_TRUE(t.defined_on_spatial_dimension(3));

    NGP_EXPECT_EQ(t.base(),stk::topology::HEX_8);

    NGP_EXPECT_EQ(t.face_topology(0), stk::topology::QUAD_8);
    NGP_EXPECT_EQ(t.face_topology(1), stk::topology::QUAD_8);
    NGP_EXPECT_EQ(t.face_topology(2), stk::topology::QUAD_8);
    NGP_EXPECT_EQ(t.face_topology(3), stk::topology::QUAD_8);
    NGP_EXPECT_EQ(t.face_topology(4), stk::topology::QUAD_8);
    NGP_EXPECT_EQ(t.face_topology(5), stk::topology::QUAD_8);

    check_edge_node_ordinals_ngp<numNodes>(t, goldEdgeNodeOrdinals);
    check_edge_nodes_ngp<numNodes>(t, goldEdgeNodeOrdinals);

    check_side_node_ordinals_ngp<numNodes>(t, goldFaceNodeOrdinals);
    check_face_node_ordinals_ngp<numNodes>(t, goldFaceNodeOrdinals);
    check_side_nodes_ngp<numNodes>(t, goldFaceNodeOrdinals);
    check_face_nodes_ngp<numNodes>(t, goldFaceNodeOrdinals);
  });

  Kokkos::parallel_for(stk::ngp::DeviceRangePolicy(0, 1), KOKKOS_LAMBDA(const int i)
  {
    check_permutation_node_ordinals_ngp<numNodes>(t, goldPermutationNodeOrdinals);
    check_permutation_nodes_ngp<numNodes>(t, goldPermutationNodeOrdinals);

    check_equivalent_ngp<numNodes>(t, goldPermutationNodeOrdinals);
    check_lexicographical_smallest_permutation_ngp<numNodes>(t, goldPermutationNodeOrdinals);
  });
}

NGP_TEST(stk_topology_ngp, hex_20)
{
  check_hex_20_on_device();
}

std::vector<std::vector<uint8_t>> get_gold_edge_node_ordinals_hex27() {
  return std::vector<std::vector<uint8_t>> {
    {0, 1, 8},
    {1, 2, 9},
    {2, 3, 10},
    {3, 0, 11},
    {4, 5, 16},
    {5, 6, 17},
    {6, 7, 18},
    {7, 4, 19},
    {0, 4, 12},
    {1, 5, 13},
    {2, 6, 14},
    {3, 7, 15}
  };
}

std::vector<std::vector<uint8_t>> get_gold_face_node_ordinals_hex27() {
  return std::vector<std::vector<uint8_t>> {
    {0, 1, 5, 4,  8, 13, 16, 12, 25},
    {1, 2, 6, 5,  9, 14, 17, 13, 24},
    {2, 3, 7, 6, 10, 15, 18, 14, 26},
    {0, 4, 7, 3, 12, 19, 15, 11, 23},
    {0, 3, 2, 1, 11, 10,  9,  8, 21},
    {4, 5, 6, 7, 16, 17, 18, 19, 22}
  };
}

std::vector<std::vector<uint8_t>> get_gold_permutation_node_ordinals_hex27() {
  return std::vector<std::vector<uint8_t>> {
    {0, 1, 2, 3, 4, 5, 6, 7,  8,  9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26},
    {0, 1, 5, 4, 3, 2, 6, 7,  8, 13, 16, 12, 11,  9, 17, 19, 10, 14, 18, 15, 20, 25, 26, 23, 24, 21, 22},
    {0, 4, 7, 3, 1, 5, 6, 2, 12, 19, 15, 11,  8, 16, 18, 10, 13, 17, 14,  9, 20, 23, 24, 21, 22, 25, 26},
    {1, 2, 3, 0, 5, 6, 7, 4,  9, 10, 11,  8, 13, 14, 15, 12, 17, 18, 19, 16, 20, 21, 22, 25, 26, 24, 23},
    {1, 2, 6, 5, 0, 3, 7, 4,  9, 14, 17, 13,  8, 10, 18, 16, 11, 15, 19, 12, 20, 24, 23, 25, 26, 21, 22},
    {1, 5, 4, 0, 2, 6, 7, 3, 13, 16, 12,  8,  9, 17, 19, 11, 14, 18, 15, 10, 20, 25, 26, 21, 22, 24, 23},
    {2, 3, 0, 1, 6, 7, 4, 5, 10, 11,  8,  9, 14, 15, 12, 13, 18, 19, 16, 17, 20, 21, 22, 24, 23, 26, 25},
    {2, 3, 7, 6, 1, 0, 4, 5, 10, 15, 18, 14,  9, 11, 19, 17,  8, 12, 16, 13, 20, 26, 25, 24, 23, 21, 22},
    {2, 6, 5, 1, 3, 7, 4, 0, 14, 17, 13,  9, 10, 18, 16,  8, 15, 19, 12, 11, 20, 24, 23, 21, 22, 26, 25},
    {3, 0, 1, 2, 7, 4, 5, 6, 11,  8,  9, 10, 15, 12, 13, 14, 19, 16, 17, 18, 20, 21, 22, 26, 25, 23, 24},
    {3, 0, 4, 7, 2, 1, 5, 6, 11, 12, 19, 15, 10,  8, 16, 18,  9, 13, 17, 14, 20, 23, 24, 26, 25, 21, 22},
    {3, 7, 6, 2, 0, 4, 5, 1, 15, 18, 14, 10, 11, 19, 17,  9, 12, 16, 13,  8, 20, 26, 25, 21, 22, 23, 24},
    {4, 0, 1, 5, 7, 3, 2, 6, 12,  8, 13, 16, 19, 11,  9, 17, 15, 10, 14, 18, 20, 25, 26, 22, 21, 23, 24},
    {4, 7, 3, 0, 5, 6, 2, 1, 19, 15, 11, 12, 16, 18, 10,  8, 17, 14,  9, 13, 20, 23, 24, 25, 26, 22, 21},
    {4, 7, 6, 5, 0, 3, 2, 1, 19, 18, 17, 16, 12, 15, 14, 13, 11, 10,  9,  8, 20, 22, 21, 25, 26, 23, 24},
    {5, 1, 2, 6, 4, 0, 3, 7, 13,  9, 14, 17, 16,  8, 10, 18, 12, 11, 15, 19, 20, 24, 23, 22, 21, 25, 26},
    {5, 4, 0, 1, 6, 7, 3, 2, 16, 12,  8, 13, 17, 19, 11,  9, 18, 15, 10, 14, 20, 25, 26, 24, 23, 22, 21},
    {5, 4, 7, 6, 1, 0, 3, 2, 16, 19, 18, 17, 13, 12, 15, 14,  8, 11, 10,  9, 20, 22, 21, 24, 23, 25, 26},
    {6, 2, 3, 7, 5, 1, 0, 4, 14, 10, 15, 18, 17,  9, 11, 19, 13,  8, 12, 16, 20, 26, 25, 22, 21, 24, 23},
    {6, 5, 1, 2, 7, 4, 0, 3, 17, 13,  9, 14, 18, 16,  8, 10, 19, 12, 11, 15, 20, 24, 23, 26, 25, 22, 21},
    {6, 5, 4, 7, 2, 1, 0, 3, 17, 16, 19, 18, 14, 13, 12, 15,  9,  8, 11, 10, 20, 22, 21, 26, 25, 24, 23},
    {7, 3, 0, 4, 6, 2, 1, 5, 15, 11, 12, 19, 18, 10,  8, 16, 14,  9, 13, 17, 20, 23, 24, 22, 21, 26, 25},
    {7, 6, 2, 3, 4, 5, 1, 0, 18, 14, 10, 15, 19, 17,  9, 11, 16, 13,  8, 12, 20, 26, 25, 23, 24, 22, 21},
    {7, 6, 5, 4, 3, 2, 1, 0, 18, 17, 16, 19, 15, 14, 13, 12, 10,  9,  8, 11, 20, 22, 21, 23, 24, 26, 25}
  };
}

TEST(stk_topology, hex_27)
{
  stk::topology t = stk::topology::HEX_27;

  EXPECT_TRUE(t.is_valid());
  EXPECT_TRUE(t.has_homogeneous_faces());
  EXPECT_FALSE(t.is_shell());

  EXPECT_EQ(t.rank(),stk::topology::ELEMENT_RANK);
  EXPECT_EQ(t.side_rank(),stk::topology::FACE_RANK);
  EXPECT_EQ(t.num_sides(),6u);

  EXPECT_EQ(t.num_nodes(),27u);
  EXPECT_EQ(t.num_vertices(),8u);
  EXPECT_EQ(t.num_edges(),12u);
  EXPECT_EQ(t.num_faces(),6u);

  EXPECT_FALSE(t.defined_on_spatial_dimension(1));
  EXPECT_FALSE(t.defined_on_spatial_dimension(2));
  EXPECT_TRUE(t.defined_on_spatial_dimension(3));

  EXPECT_EQ(t.base(),stk::topology::HEX_8);

  EXPECT_EQ(t.face_topology(0), stk::topology::QUAD_9);
  EXPECT_EQ(t.face_topology(1), stk::topology::QUAD_9);
  EXPECT_EQ(t.face_topology(2), stk::topology::QUAD_9);
  EXPECT_EQ(t.face_topology(3), stk::topology::QUAD_9);
  EXPECT_EQ(t.face_topology(4), stk::topology::QUAD_9);
  EXPECT_EQ(t.face_topology(5), stk::topology::QUAD_9);

  check_edge_node_ordinals(t, get_gold_edge_node_ordinals_hex27());
  check_edge_nodes(t, get_gold_edge_node_ordinals_hex27());

  check_side_node_ordinals(t, get_gold_face_node_ordinals_hex27());
  check_face_node_ordinals(t, get_gold_face_node_ordinals_hex27());
  check_side_nodes(t, get_gold_face_node_ordinals_hex27());
  check_face_nodes(t, get_gold_face_node_ordinals_hex27());

  check_permutation_node_ordinals(t, get_gold_permutation_node_ordinals_hex27());
  check_permutation_nodes(t, get_gold_permutation_node_ordinals_hex27());

  check_equivalent(t, get_gold_permutation_node_ordinals_hex27());
  check_lexicographical_smallest_permutation(t, get_gold_permutation_node_ordinals_hex27());
}

void check_hex_27_on_device()
{
  OrdinalType goldEdgeNodeOrdinals = fillGoldOrdinals(get_gold_edge_node_ordinals_hex27());
  OrdinalType goldFaceNodeOrdinals = fillGoldOrdinals(get_gold_face_node_ordinals_hex27());
  OrdinalType goldPermutationNodeOrdinals = fillGoldOrdinals(get_gold_permutation_node_ordinals_hex27());

  const stk::topology t = stk::topology::HEX_27;
  const unsigned numNodes = stk::topology_detail::topology_data<stk::topology::HEX_27>::num_nodes;
  EXPECT_EQ(27u, numNodes);

  Kokkos::parallel_for(stk::ngp::DeviceRangePolicy(0, 1), KOKKOS_LAMBDA(const int i)
  {
    NGP_EXPECT_TRUE(t.is_valid());
    NGP_EXPECT_TRUE(t.has_homogeneous_faces());
    NGP_EXPECT_FALSE(t.is_shell());

    NGP_EXPECT_EQ(t.rank(),stk::topology::ELEMENT_RANK);
    NGP_EXPECT_EQ(t.side_rank(),stk::topology::FACE_RANK);
    NGP_EXPECT_EQ(t.num_sides(),6u);

    NGP_EXPECT_EQ(t.num_nodes(),27u);
    NGP_EXPECT_EQ(t.num_vertices(),8u);
    NGP_EXPECT_EQ(t.num_edges(),12u);
    NGP_EXPECT_EQ(t.num_faces(),6u);

    NGP_EXPECT_FALSE(t.defined_on_spatial_dimension(1));
    NGP_EXPECT_FALSE(t.defined_on_spatial_dimension(2));
    NGP_EXPECT_TRUE(t.defined_on_spatial_dimension(3));

    NGP_EXPECT_EQ(t.base(),stk::topology::HEX_8);

    NGP_EXPECT_EQ(t.face_topology(0), stk::topology::QUAD_9);
    NGP_EXPECT_EQ(t.face_topology(1), stk::topology::QUAD_9);
    NGP_EXPECT_EQ(t.face_topology(2), stk::topology::QUAD_9);
    NGP_EXPECT_EQ(t.face_topology(3), stk::topology::QUAD_9);
    NGP_EXPECT_EQ(t.face_topology(4), stk::topology::QUAD_9);
    NGP_EXPECT_EQ(t.face_topology(5), stk::topology::QUAD_9);

    check_edge_node_ordinals_ngp<numNodes>(t, goldEdgeNodeOrdinals);
    check_edge_nodes_ngp<numNodes>(t, goldEdgeNodeOrdinals);

    check_side_node_ordinals_ngp<numNodes>(t, goldFaceNodeOrdinals);
    check_face_node_ordinals_ngp<numNodes>(t, goldFaceNodeOrdinals);
    check_side_nodes_ngp<numNodes>(t, goldFaceNodeOrdinals);
    check_face_nodes_ngp<numNodes>(t, goldFaceNodeOrdinals);
  });

  Kokkos::parallel_for(stk::ngp::DeviceRangePolicy(0, 1), KOKKOS_LAMBDA(const int i)
  {
    check_permutation_node_ordinals_ngp<numNodes>(t, goldPermutationNodeOrdinals);
    check_permutation_nodes_ngp<numNodes>(t, goldPermutationNodeOrdinals);

    check_equivalent_ngp<numNodes>(t, goldPermutationNodeOrdinals);
    check_lexicographical_smallest_permutation_ngp<numNodes>(t, goldPermutationNodeOrdinals);
  });
}

NGP_TEST(stk_topology_ngp, hex_27)
{
  check_hex_27_on_device();
}

}

