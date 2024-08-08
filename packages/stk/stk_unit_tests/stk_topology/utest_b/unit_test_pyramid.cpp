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
#include "stk_topology/topology.hpp"  // for topology, topology::TRI_6, topology::PYRAMID_5, top...
#include "topology_test_utils.hpp"    // for INVALID, check_edge_node_ordinals, check_edge_node_...
#include <vector>                     // for vector

namespace {

std::vector<std::vector<uint8_t>> get_gold_edge_node_ordinals_pyramid5() {
  return std::vector<std::vector<uint8_t>> {
    {0, 1},
    {1, 2},
    {2, 3},
    {3, 0},
    {0, 4},
    {1, 4},
    {2, 4},
    {3, 4}
  };
}

std::vector<std::vector<uint8_t>> get_gold_face_node_ordinals_pyramid5() {
  return std::vector<std::vector<uint8_t>> {
    {0, 1, 4},
    {1, 2, 4},
    {2, 3, 4},
    {0, 4, 3},
    {0, 3, 2, 1}
  };
}

std::vector<std::vector<uint8_t>> get_gold_permutation_node_ordinals_pyramid5() {
  return std::vector<std::vector<uint8_t>> {
    {0, 1, 2, 3, 4},
    {1, 2, 3, 0, 4},
    {2, 3, 0, 1, 4},
    {3, 0, 1, 2, 4}
  };
}

TEST(stk_topology, pyramid_5)
{
  stk::topology t = stk::topology::PYRAMID_5;

  EXPECT_TRUE(t.is_valid());
  EXPECT_FALSE(t.has_homogeneous_faces());
  EXPECT_FALSE(t.is_shell());

  EXPECT_EQ(t.rank(),stk::topology::ELEMENT_RANK);
  EXPECT_EQ(t.side_rank(),stk::topology::FACE_RANK);
  EXPECT_EQ(t.num_sides(),5u);

  EXPECT_EQ(t.num_nodes(),5u);
  EXPECT_EQ(t.num_vertices(),5u);
  EXPECT_EQ(t.num_edges(),8u);
  EXPECT_EQ(t.num_faces(),5u);

  EXPECT_FALSE(t.defined_on_spatial_dimension(1));
  EXPECT_FALSE(t.defined_on_spatial_dimension(2));
  EXPECT_TRUE(t.defined_on_spatial_dimension(3));

  EXPECT_EQ(t.base(),stk::topology::PYRAMID_5);

  EXPECT_EQ(t.face_topology(0), stk::topology::TRI_3);
  EXPECT_EQ(t.face_topology(1), stk::topology::TRI_3);
  EXPECT_EQ(t.face_topology(2), stk::topology::TRI_3);
  EXPECT_EQ(t.face_topology(3), stk::topology::TRI_3);
  EXPECT_EQ(t.face_topology(4), stk::topology::QUAD_4);

  check_edge_node_ordinals(t, get_gold_edge_node_ordinals_pyramid5());
  check_edge_nodes(t, get_gold_edge_node_ordinals_pyramid5());

  check_side_node_ordinals(t, get_gold_face_node_ordinals_pyramid5());
  check_face_node_ordinals(t, get_gold_face_node_ordinals_pyramid5());
  check_side_nodes(t, get_gold_face_node_ordinals_pyramid5());
  check_face_nodes(t, get_gold_face_node_ordinals_pyramid5());

  check_permutation_node_ordinals(t, get_gold_permutation_node_ordinals_pyramid5());
  check_permutation_nodes(t, get_gold_permutation_node_ordinals_pyramid5());

  check_equivalent(t, get_gold_permutation_node_ordinals_pyramid5());
  check_lexicographical_smallest_permutation(t, get_gold_permutation_node_ordinals_pyramid5());
}

void check_pyramid_5_on_device()
{
  OrdinalType goldEdgeNodeOrdinals = fillGoldOrdinals(get_gold_edge_node_ordinals_pyramid5());
  OrdinalType goldFaceNodeOrdinals = fillGoldOrdinals(get_gold_face_node_ordinals_pyramid5());
  OrdinalType goldPermutationNodeOrdinals = fillGoldOrdinals(get_gold_permutation_node_ordinals_pyramid5());

  stk::topology t = stk::topology::PYRAMID_5;
  constexpr unsigned numNodes = stk::topology_detail::topology_data<stk::topology::PYRAMID_5>::num_nodes;
  EXPECT_EQ(5u, numNodes);

  Kokkos::parallel_for(stk::ngp::DeviceRangePolicy(0, 1), KOKKOS_LAMBDA(const int i)
  {
    NGP_EXPECT_TRUE(t.is_valid());
    NGP_EXPECT_FALSE(t.has_homogeneous_faces());
    NGP_EXPECT_FALSE(t.is_shell());

    NGP_EXPECT_EQ(t.rank(),stk::topology::ELEMENT_RANK);
    NGP_EXPECT_EQ(t.side_rank(),stk::topology::FACE_RANK);
    NGP_EXPECT_EQ(t.num_sides(),5u);

    NGP_EXPECT_EQ(t.num_nodes(),5u);
    NGP_EXPECT_EQ(t.num_vertices(),5u);
    NGP_EXPECT_EQ(t.num_edges(),8u);
    NGP_EXPECT_EQ(t.num_faces(),5u);

    NGP_EXPECT_FALSE(t.defined_on_spatial_dimension(1));
    NGP_EXPECT_FALSE(t.defined_on_spatial_dimension(2));
    NGP_EXPECT_TRUE(t.defined_on_spatial_dimension(3));

    NGP_EXPECT_EQ(t.base(),stk::topology::PYRAMID_5);

    NGP_EXPECT_EQ(t.face_topology(0), stk::topology::TRI_3);
    NGP_EXPECT_EQ(t.face_topology(1), stk::topology::TRI_3);
    NGP_EXPECT_EQ(t.face_topology(2), stk::topology::TRI_3);
    NGP_EXPECT_EQ(t.face_topology(3), stk::topology::TRI_3);
    NGP_EXPECT_EQ(t.face_topology(4), stk::topology::QUAD_4);

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

NGP_TEST(stk_topology_ngp, pyramid_5)
{
  check_pyramid_5_on_device();
}


std::vector<std::vector<uint8_t>> get_gold_edge_node_ordinals_pyramid13() {
  return std::vector<std::vector<uint8_t>> {
    {0, 1, 5},
    {1, 2, 6},
    {2, 3, 7},
    {3, 0, 8},
    {0, 4, 9},
    {1, 4, 10},
    {2, 4, 11},
    {3, 4, 12}
  };
}

std::vector<std::vector<uint8_t>> get_gold_face_node_ordinals_pyramid13() {
  return std::vector<std::vector<uint8_t>> {
    {0, 1, 4, 5, 10,  9},
    {1, 2, 4, 6, 11, 10},
    {2, 3, 4, 7, 12, 11},
    {3, 0, 4, 8,  9, 12},
    {0, 3, 2, 1,  8,  7, 6, 5}
  };
}

std::vector<std::vector<uint8_t>> get_gold_permutation_node_ordinals_pyramid13() {
  return std::vector<std::vector<uint8_t>> {
    {0, 1, 2, 3, 4, 5, 6, 7, 8,  9, 10, 11, 12},
    {1, 2, 3, 0, 4, 6, 7, 8, 5, 10, 11, 12,  9},
    {2, 3, 0, 1, 4, 7, 8, 5, 6, 11, 12,  9, 10},
    {3, 0, 1, 2, 4, 8, 5, 6, 7, 12,  9, 10, 11}
  };
}

TEST(stk_topology, pyramid_13)
{
  stk::topology t = stk::topology::PYRAMID_13;

  EXPECT_TRUE(t.is_valid());
  EXPECT_FALSE(t.has_homogeneous_faces());
  EXPECT_FALSE(t.is_shell());

  EXPECT_EQ(t.rank(),stk::topology::ELEMENT_RANK);
  EXPECT_EQ(t.side_rank(),stk::topology::FACE_RANK);
  EXPECT_EQ(t.num_sides(),5u);

  EXPECT_EQ(t.num_nodes(),13u);
  EXPECT_EQ(t.num_vertices(),5u);
  EXPECT_EQ(t.num_edges(),8u);
  EXPECT_EQ(t.num_faces(),5u);

  EXPECT_FALSE(t.defined_on_spatial_dimension(1));
  EXPECT_FALSE(t.defined_on_spatial_dimension(2));
  EXPECT_TRUE(t.defined_on_spatial_dimension(3));

  EXPECT_EQ(t.base(),stk::topology::PYRAMID_5);

  EXPECT_EQ(t.face_topology(0), stk::topology::TRI_6);
  EXPECT_EQ(t.face_topology(1), stk::topology::TRI_6);
  EXPECT_EQ(t.face_topology(2), stk::topology::TRI_6);
  EXPECT_EQ(t.face_topology(3), stk::topology::TRI_6);
  EXPECT_EQ(t.face_topology(4), stk::topology::QUAD_8);

  check_edge_node_ordinals(t, get_gold_edge_node_ordinals_pyramid13());
  check_edge_nodes(t, get_gold_edge_node_ordinals_pyramid13());

  check_side_node_ordinals(t, get_gold_face_node_ordinals_pyramid13());
  check_face_node_ordinals(t, get_gold_face_node_ordinals_pyramid13());
  check_side_nodes(t, get_gold_face_node_ordinals_pyramid13());
  check_face_nodes(t, get_gold_face_node_ordinals_pyramid13());

  check_permutation_node_ordinals(t, get_gold_permutation_node_ordinals_pyramid13());
  check_permutation_nodes(t, get_gold_permutation_node_ordinals_pyramid13());

  check_equivalent(t, get_gold_permutation_node_ordinals_pyramid13());
  check_lexicographical_smallest_permutation(t, get_gold_permutation_node_ordinals_pyramid13());
}

void check_pyramid_13_on_device()
{
  OrdinalType goldEdgeNodeOrdinals = fillGoldOrdinals(get_gold_edge_node_ordinals_pyramid13());
  OrdinalType goldFaceNodeOrdinals = fillGoldOrdinals(get_gold_face_node_ordinals_pyramid13());
  OrdinalType goldPermutationNodeOrdinals = fillGoldOrdinals(get_gold_permutation_node_ordinals_pyramid13());

  stk::topology t = stk::topology::PYRAMID_13;
  constexpr unsigned numNodes = stk::topology_detail::topology_data<stk::topology::PYRAMID_13>::num_nodes;
  EXPECT_EQ(13u, numNodes);

  Kokkos::parallel_for(stk::ngp::DeviceRangePolicy(0, 1), KOKKOS_LAMBDA(const int i)
  {
    NGP_EXPECT_TRUE(t.is_valid());
    NGP_EXPECT_FALSE(t.has_homogeneous_faces());
    NGP_EXPECT_FALSE(t.is_shell());

    NGP_EXPECT_EQ(t.rank(),stk::topology::ELEMENT_RANK);
    NGP_EXPECT_EQ(t.side_rank(),stk::topology::FACE_RANK);
    NGP_EXPECT_EQ(t.num_sides(),5u);

    NGP_EXPECT_EQ(t.num_nodes(),13u);
    NGP_EXPECT_EQ(t.num_vertices(),5u);
    NGP_EXPECT_EQ(t.num_edges(),8u);
    NGP_EXPECT_EQ(t.num_faces(),5u);

    NGP_EXPECT_FALSE(t.defined_on_spatial_dimension(1));
    NGP_EXPECT_FALSE(t.defined_on_spatial_dimension(2));
    NGP_EXPECT_TRUE(t.defined_on_spatial_dimension(3));

    NGP_EXPECT_EQ(t.base(),stk::topology::PYRAMID_5);

    NGP_EXPECT_EQ(t.face_topology(0), stk::topology::TRI_6);
    NGP_EXPECT_EQ(t.face_topology(1), stk::topology::TRI_6);
    NGP_EXPECT_EQ(t.face_topology(2), stk::topology::TRI_6);
    NGP_EXPECT_EQ(t.face_topology(3), stk::topology::TRI_6);
    NGP_EXPECT_EQ(t.face_topology(4), stk::topology::QUAD_8);

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

NGP_TEST(stk_topology_ngp, pyramid_13)
{
  check_pyramid_13_on_device();
}


std::vector<std::vector<uint8_t>> get_gold_edge_node_ordinals_pyramid14() {
  return std::vector<std::vector<uint8_t>> {
    {0, 1, 5},
    {1, 2, 6},
    {2, 3, 7},
    {3, 0, 8},
    {0, 4, 9},
    {1, 4, 10},
    {2, 4, 11},
    {3, 4, 12}
  };
}

std::vector<std::vector<uint8_t>> get_gold_face_node_ordinals_pyramid14() {
  return std::vector<std::vector<uint8_t>> {
    {0, 1, 4, 5, 10,  9},
    {1, 2, 4, 6, 11, 10},
    {2, 3, 4, 7, 12, 11},
    {3, 0, 4, 8,  9, 12},
    {0, 3, 2, 1,  8,  7, 6, 5, 13}
  };
}

std::vector<std::vector<uint8_t>> get_gold_permutation_node_ordinals_pyramid14() {
  return std::vector<std::vector<uint8_t>> {
    {0, 1, 2, 3, 4, 5, 6, 7, 8,  9, 10, 11, 12, 13},
    {1, 2, 3, 0, 4, 6, 7, 8, 5, 10, 11, 12,  9, 13},
    {2, 3, 0, 1, 4, 7, 8, 5, 6, 11, 12,  9, 10, 13},
    {3, 0, 1, 2, 4, 8, 5, 6, 7, 12,  9, 10, 11, 13}
  };
}

TEST(stk_topology, pyramid_14)
{
  stk::topology t = stk::topology::PYRAMID_14;

  EXPECT_TRUE(t.is_valid());
  EXPECT_FALSE(t.has_homogeneous_faces());
  EXPECT_FALSE(t.is_shell());

  EXPECT_EQ(t.rank(),stk::topology::ELEMENT_RANK);
  EXPECT_EQ(t.side_rank(),stk::topology::FACE_RANK);
  EXPECT_EQ(t.num_sides(),5u);

  EXPECT_EQ(t.num_nodes(),14u);
  EXPECT_EQ(t.num_vertices(),5u);
  EXPECT_EQ(t.num_edges(),8u);
  EXPECT_EQ(t.num_faces(),5u);

  EXPECT_FALSE(t.defined_on_spatial_dimension(1));
  EXPECT_FALSE(t.defined_on_spatial_dimension(2));
  EXPECT_TRUE(t.defined_on_spatial_dimension(3));

  EXPECT_EQ(t.base(),stk::topology::PYRAMID_5);

  EXPECT_EQ(t.face_topology(0), stk::topology::TRI_6);
  EXPECT_EQ(t.face_topology(1), stk::topology::TRI_6);
  EXPECT_EQ(t.face_topology(2), stk::topology::TRI_6);
  EXPECT_EQ(t.face_topology(3), stk::topology::TRI_6);
  EXPECT_EQ(t.face_topology(4), stk::topology::QUAD_9);

  check_edge_node_ordinals(t, get_gold_edge_node_ordinals_pyramid14());
  check_edge_nodes(t, get_gold_edge_node_ordinals_pyramid14());

  check_side_node_ordinals(t, get_gold_face_node_ordinals_pyramid14());
  check_face_node_ordinals(t, get_gold_face_node_ordinals_pyramid14());
  check_side_nodes(t, get_gold_face_node_ordinals_pyramid14());
  check_face_nodes(t, get_gold_face_node_ordinals_pyramid14());

  check_permutation_node_ordinals(t, get_gold_permutation_node_ordinals_pyramid14());
  check_permutation_nodes(t, get_gold_permutation_node_ordinals_pyramid14());

  check_equivalent(t, get_gold_permutation_node_ordinals_pyramid14());
  check_lexicographical_smallest_permutation(t, get_gold_permutation_node_ordinals_pyramid14());
}

void check_pyramid_14_on_device()
{
  OrdinalType goldEdgeNodeOrdinals = fillGoldOrdinals(get_gold_edge_node_ordinals_pyramid14());
  OrdinalType goldFaceNodeOrdinals = fillGoldOrdinals(get_gold_face_node_ordinals_pyramid14());
  OrdinalType goldPermutationNodeOrdinals = fillGoldOrdinals(get_gold_permutation_node_ordinals_pyramid14());

  stk::topology t = stk::topology::PYRAMID_14;
  constexpr unsigned numNodes = stk::topology_detail::topology_data<stk::topology::PYRAMID_14>::num_nodes;
  EXPECT_EQ(14u, numNodes);

  Kokkos::parallel_for(stk::ngp::DeviceRangePolicy(0, 1), KOKKOS_LAMBDA(const int i)
  {
    NGP_EXPECT_TRUE(t.is_valid());
    NGP_EXPECT_FALSE(t.has_homogeneous_faces());
    NGP_EXPECT_FALSE(t.is_shell());

    NGP_EXPECT_EQ(t.rank(),stk::topology::ELEMENT_RANK);
    NGP_EXPECT_EQ(t.side_rank(),stk::topology::FACE_RANK);
    NGP_EXPECT_EQ(t.num_sides(),5u);

    NGP_EXPECT_EQ(t.num_nodes(),14u);
    NGP_EXPECT_EQ(t.num_vertices(),5u);
    NGP_EXPECT_EQ(t.num_edges(),8u);
    NGP_EXPECT_EQ(t.num_faces(),5u);

    NGP_EXPECT_FALSE(t.defined_on_spatial_dimension(1));
    NGP_EXPECT_FALSE(t.defined_on_spatial_dimension(2));
    NGP_EXPECT_TRUE(t.defined_on_spatial_dimension(3));

    NGP_EXPECT_EQ(t.base(),stk::topology::PYRAMID_5);

    NGP_EXPECT_EQ(t.face_topology(0), stk::topology::TRI_6);
    NGP_EXPECT_EQ(t.face_topology(1), stk::topology::TRI_6);
    NGP_EXPECT_EQ(t.face_topology(2), stk::topology::TRI_6);
    NGP_EXPECT_EQ(t.face_topology(3), stk::topology::TRI_6);
    NGP_EXPECT_EQ(t.face_topology(4), stk::topology::QUAD_9);

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

NGP_TEST(stk_topology_ngp, pyramid_14)
{
  check_pyramid_14_on_device();
}

}
