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
#include "topology_test_utils.hpp"           // for check_edge_node_ordinals, check_edge_node_or...
#include <cstddef>                           // for size_t
#include <vector>                            // for vector

namespace {

std::vector<std::vector<uint8_t>> get_gold_edge_node_ordinals_tet4() {
  return std::vector<std::vector<uint8_t>> {
    {0, 1},
    {1, 2},
    {2, 0},
    {0, 3},
    {1, 3},
    {2, 3}
  };
}

std::vector<std::vector<uint8_t>> get_gold_face_node_ordinals_tet4() {
  return std::vector<std::vector<uint8_t>> {
    {0, 1, 3},
    {1, 2, 3},
    {0, 3, 2},
    {0, 2, 1}
  };
}

std::vector<std::vector<uint8_t>> get_gold_permutation_node_ordinals_tet4() {
  return std::vector<std::vector<uint8_t>> {
    {0, 1, 2, 3},
    {1, 2, 0, 3},
    {2, 0, 1, 3},
    {0, 3, 1, 2},
    {3, 1, 0, 2},
    {1, 0, 3, 2},
    {0, 2, 3, 1},
    {2, 3, 0, 1},
    {3, 0, 2, 1},
    {1, 3, 2, 0},
    {3, 2, 1, 0},
    {2, 1, 3, 0}
  };
}

TEST(stk_topology, tet_4)
{
  stk::topology t = stk::topology::TET_4;

  EXPECT_TRUE(t.is_valid());
  EXPECT_TRUE(t.has_homogeneous_faces());
  EXPECT_FALSE(t.is_shell());

  EXPECT_EQ(t.rank(), stk::topology::ELEMENT_RANK);
  EXPECT_EQ(t.side_rank(), stk::topology::FACE_RANK);
  EXPECT_EQ(t.num_sides(), 4u);

  EXPECT_EQ(t.num_nodes(), 4u);
  EXPECT_EQ(t.num_vertices(), 4u);
  EXPECT_EQ(t.num_edges(), 6u);
  EXPECT_EQ(t.num_faces(), 4u);

  EXPECT_FALSE(t.defined_on_spatial_dimension(1));
  EXPECT_FALSE(t.defined_on_spatial_dimension(2));
  EXPECT_TRUE(t.defined_on_spatial_dimension(3));

  EXPECT_EQ(t.base(), stk::topology::TET_4);

  EXPECT_EQ(t.face_topology(0), stk::topology::TRI_3);
  EXPECT_EQ(t.face_topology(1), stk::topology::TRI_3);
  EXPECT_EQ(t.face_topology(2), stk::topology::TRI_3);
  EXPECT_EQ(t.face_topology(3), stk::topology::TRI_3);

  check_edge_node_ordinals(t, get_gold_edge_node_ordinals_tet4());
  check_edge_nodes(t, get_gold_edge_node_ordinals_tet4());

  check_side_node_ordinals(t, get_gold_face_node_ordinals_tet4());
  check_face_node_ordinals(t, get_gold_face_node_ordinals_tet4());
  check_side_nodes(t, get_gold_face_node_ordinals_tet4());
  check_face_nodes(t, get_gold_face_node_ordinals_tet4());

  check_permutation_node_ordinals(t, get_gold_permutation_node_ordinals_tet4());
  check_permutation_nodes(t, get_gold_permutation_node_ordinals_tet4());

  check_equivalent(t, get_gold_permutation_node_ordinals_tet4());
  check_lexicographical_smallest_permutation(t, get_gold_permutation_node_ordinals_tet4());
}

void check_tet_4_on_device()
{
  OrdinalType goldEdgeNodeOrdinals = fillGoldOrdinals(get_gold_edge_node_ordinals_tet4());
  OrdinalType goldFaceNodeOrdinals = fillGoldOrdinals(get_gold_face_node_ordinals_tet4());
  OrdinalType goldPermutationNodeOrdinals = fillGoldOrdinals(get_gold_permutation_node_ordinals_tet4());

  stk::topology t = stk::topology::TET_4;
  constexpr unsigned numNodes = stk::topology_detail::topology_data<stk::topology::TET_4>::num_nodes;
  EXPECT_EQ(4u, numNodes);

  Kokkos::parallel_for(stk::ngp::DeviceRangePolicy(0, 1), KOKKOS_LAMBDA(const int i)
  {
    NGP_EXPECT_TRUE(t.is_valid());
    NGP_EXPECT_TRUE(t.has_homogeneous_faces());
    NGP_EXPECT_FALSE(t.is_shell());

    NGP_EXPECT_EQ(t.rank(), stk::topology::ELEMENT_RANK);
    NGP_EXPECT_EQ(t.side_rank(), stk::topology::FACE_RANK);
    NGP_EXPECT_EQ(t.num_sides(), 4u);

    NGP_EXPECT_EQ(t.num_nodes(), 4u);
    NGP_EXPECT_EQ(t.num_vertices(), 4u);
    NGP_EXPECT_EQ(t.num_edges(), 6u);
    NGP_EXPECT_EQ(t.num_faces(), 4u);

    NGP_EXPECT_FALSE(t.defined_on_spatial_dimension(1));
    NGP_EXPECT_FALSE(t.defined_on_spatial_dimension(2));
    NGP_EXPECT_TRUE(t.defined_on_spatial_dimension(3));

    NGP_EXPECT_EQ(t.base(), stk::topology::TET_4);

    NGP_EXPECT_EQ(t.face_topology(0), stk::topology::TRI_3);
    NGP_EXPECT_EQ(t.face_topology(1), stk::topology::TRI_3);
    NGP_EXPECT_EQ(t.face_topology(2), stk::topology::TRI_3);
    NGP_EXPECT_EQ(t.face_topology(3), stk::topology::TRI_3);

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

NGP_TEST(stk_topology_ngp, tet_4)
{
  check_tet_4_on_device();
}


std::vector<std::vector<uint8_t>> get_gold_edge_node_ordinals_tet8() {
  return std::vector<std::vector<uint8_t>> {
    {0, 1},
    {1, 2},
    {2, 0},
    {0, 3},
    {1, 3},
    {2, 3}
  };
}

std::vector<std::vector<uint8_t>> get_gold_face_node_ordinals_tet8() {
  return std::vector<std::vector<uint8_t>> {
    {0, 1, 3,  4},
    {1, 2, 3,  5},
    {0, 3, 2,  7},
    {0, 2, 1,  6}
  };
}

std::vector<std::vector<uint8_t>> get_gold_permutation_node_ordinals_tet8() {
  return std::vector<std::vector<uint8_t>> {
    {0, 1, 2, 3,  4, 5, 6, 7}
  };
}

TEST(stk_topology, tet_8)
{
  stk::topology t = stk::topology::TET_8;

  EXPECT_TRUE(t.is_valid());
  EXPECT_TRUE(t.has_homogeneous_faces());
  EXPECT_FALSE(t.is_shell());

  EXPECT_EQ(t.rank(), stk::topology::ELEMENT_RANK);
  EXPECT_EQ(t.side_rank(), stk::topology::FACE_RANK);
  EXPECT_EQ(t.num_sides(), 4u);

  EXPECT_EQ(t.num_nodes(), 8u);
  EXPECT_EQ(t.num_vertices(), 4u);
  EXPECT_EQ(t.num_edges(), 6u);
  EXPECT_EQ(t.num_faces(), 4u);

  EXPECT_FALSE(t.defined_on_spatial_dimension(1));
  EXPECT_FALSE(t.defined_on_spatial_dimension(2));
  EXPECT_TRUE(t.defined_on_spatial_dimension(3));

  EXPECT_EQ(t.base(), stk::topology::TET_4);

  EXPECT_EQ(t.face_topology(0), stk::topology::TRI_4);
  EXPECT_EQ(t.face_topology(1), stk::topology::TRI_4);
  EXPECT_EQ(t.face_topology(2), stk::topology::TRI_4);
  EXPECT_EQ(t.face_topology(3), stk::topology::TRI_4);

  check_edge_node_ordinals(t, get_gold_edge_node_ordinals_tet8());
  check_edge_nodes(t, get_gold_edge_node_ordinals_tet8());

  check_side_node_ordinals(t, get_gold_face_node_ordinals_tet8());
  check_face_node_ordinals(t, get_gold_face_node_ordinals_tet8());
  check_side_nodes(t, get_gold_face_node_ordinals_tet8());
  check_face_nodes(t, get_gold_face_node_ordinals_tet8());

  check_permutation_node_ordinals(t, get_gold_permutation_node_ordinals_tet8());
  check_permutation_nodes(t, get_gold_permutation_node_ordinals_tet8());

  check_equivalent(t, get_gold_permutation_node_ordinals_tet8());
  check_lexicographical_smallest_permutation(t, get_gold_permutation_node_ordinals_tet8());
}

void check_tet_8_on_device()
{
  OrdinalType goldEdgeNodeOrdinals = fillGoldOrdinals(get_gold_edge_node_ordinals_tet8());
  OrdinalType goldFaceNodeOrdinals = fillGoldOrdinals(get_gold_face_node_ordinals_tet8());
  OrdinalType goldPermutationNodeOrdinals = fillGoldOrdinals(get_gold_permutation_node_ordinals_tet8());

  stk::topology t = stk::topology::TET_8;
  constexpr unsigned numNodes = stk::topology_detail::topology_data<stk::topology::TET_8>::num_nodes;
  EXPECT_EQ(8u, numNodes);

  Kokkos::parallel_for(stk::ngp::DeviceRangePolicy(0, 1), KOKKOS_LAMBDA(const int i)
  {
    NGP_EXPECT_TRUE(t.is_valid());
    NGP_EXPECT_TRUE(t.has_homogeneous_faces());
    NGP_EXPECT_FALSE(t.is_shell());

    NGP_EXPECT_EQ(t.rank(), stk::topology::ELEMENT_RANK);
    NGP_EXPECT_EQ(t.side_rank(), stk::topology::FACE_RANK);
    NGP_EXPECT_EQ(t.num_sides(), 4u);

    NGP_EXPECT_EQ(t.num_nodes(), 8u);
    NGP_EXPECT_EQ(t.num_vertices(), 4u);
    NGP_EXPECT_EQ(t.num_edges(), 6u);
    NGP_EXPECT_EQ(t.num_faces(), 4u);

    NGP_EXPECT_FALSE(t.defined_on_spatial_dimension(1));
    NGP_EXPECT_FALSE(t.defined_on_spatial_dimension(2));
    NGP_EXPECT_TRUE(t.defined_on_spatial_dimension(3));

    NGP_EXPECT_EQ(t.base(), stk::topology::TET_4);

    NGP_EXPECT_EQ(t.face_topology(0), stk::topology::TRI_4);
    NGP_EXPECT_EQ(t.face_topology(1), stk::topology::TRI_4);
    NGP_EXPECT_EQ(t.face_topology(2), stk::topology::TRI_4);
    NGP_EXPECT_EQ(t.face_topology(3), stk::topology::TRI_4);

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

NGP_TEST(stk_topology_ngp, tet_8)
{
  check_tet_8_on_device();
}


std::vector<std::vector<uint8_t>> get_gold_edge_node_ordinals_tet10() {
  return std::vector<std::vector<uint8_t>> {
    {0, 1,  4},
    {1, 2,  5},
    {2, 0,  6},
    {0, 3,  7},
    {1, 3,  8},
    {2, 3,  9}
  };
}

std::vector<std::vector<uint8_t>> get_gold_face_node_ordinals_tet10() {
  return std::vector<std::vector<uint8_t>> {
    {0, 1, 3,  4, 8, 7},
    {1, 2, 3,  5, 9, 8},
    {0, 3, 2,  7, 9, 6},
    {0, 2, 1,  6, 5, 4}
  };
}

std::vector<std::vector<uint8_t>> get_gold_permutation_node_ordinals_tet10() {
  return std::vector<std::vector<uint8_t>> {
    {0, 1, 2, 3,  4, 5, 6, 7, 8, 9},
    {1, 2, 0, 3,  5, 6, 4, 8, 9, 7},
    {2, 0, 1, 3,  6, 4, 5, 9, 7, 8},
    {0, 3, 1, 2,  7, 8, 4, 6, 9, 5},
    {3, 1, 0, 2,  8, 4, 7, 9, 5, 6},
    {1, 0, 3, 2,  4, 7, 8, 5, 6, 9},
    {0, 2, 3, 1,  6, 9, 7, 4, 5, 8},
    {2, 3, 0, 1,  9, 7, 6, 5, 8, 4},
    {3, 0, 2, 1,  7, 6, 9, 8, 4, 5},
    {1, 3, 2, 0,  8, 9, 5, 4, 7, 6},
    {3, 2, 1, 0,  9, 5, 8, 7, 6, 4},
    {2, 1, 3, 0,  5, 8, 9, 6, 4, 7}
  };
}

TEST(stk_topology, tet_10)
{
  stk::topology t = stk::topology::TET_10;

  EXPECT_TRUE(t.is_valid());
  EXPECT_TRUE(t.has_homogeneous_faces());
  EXPECT_FALSE(t.is_shell());

  EXPECT_EQ(t.rank(), stk::topology::ELEMENT_RANK);
  EXPECT_EQ(t.side_rank(), stk::topology::FACE_RANK);
  EXPECT_EQ(t.num_sides(), 4u);

  EXPECT_EQ(t.num_nodes(), 10u);
  EXPECT_EQ(t.num_vertices(), 4u);
  EXPECT_EQ(t.num_edges(), 6u);
  EXPECT_EQ(t.num_faces(), 4u);

  EXPECT_FALSE(t.defined_on_spatial_dimension(1));
  EXPECT_FALSE(t.defined_on_spatial_dimension(2));
  EXPECT_TRUE(t.defined_on_spatial_dimension(3));

  EXPECT_EQ(t.base(), stk::topology::TET_4);

  EXPECT_EQ(t.face_topology(0), stk::topology::TRI_6);
  EXPECT_EQ(t.face_topology(1), stk::topology::TRI_6);
  EXPECT_EQ(t.face_topology(2), stk::topology::TRI_6);
  EXPECT_EQ(t.face_topology(3), stk::topology::TRI_6);

  check_edge_node_ordinals(t, get_gold_edge_node_ordinals_tet10());
  check_edge_nodes(t, get_gold_edge_node_ordinals_tet10());

  check_side_node_ordinals(t, get_gold_face_node_ordinals_tet10());
  check_face_node_ordinals(t, get_gold_face_node_ordinals_tet10());
  check_side_nodes(t, get_gold_face_node_ordinals_tet10());
  check_face_nodes(t, get_gold_face_node_ordinals_tet10());

  check_permutation_node_ordinals(t, get_gold_permutation_node_ordinals_tet10());
  check_permutation_nodes(t, get_gold_permutation_node_ordinals_tet10());

  check_equivalent(t, get_gold_permutation_node_ordinals_tet10());
  check_lexicographical_smallest_permutation(t, get_gold_permutation_node_ordinals_tet10());
}

void check_tet_10_on_device()
{
  OrdinalType goldEdgeNodeOrdinals = fillGoldOrdinals(get_gold_edge_node_ordinals_tet10());
  OrdinalType goldFaceNodeOrdinals = fillGoldOrdinals(get_gold_face_node_ordinals_tet10());
  OrdinalType goldPermutationNodeOrdinals = fillGoldOrdinals(get_gold_permutation_node_ordinals_tet10());

  stk::topology t = stk::topology::TET_10;
  constexpr unsigned numNodes = stk::topology_detail::topology_data<stk::topology::TET_10>::num_nodes;
  EXPECT_EQ(10u, numNodes);

  Kokkos::parallel_for(stk::ngp::DeviceRangePolicy(0, 1), KOKKOS_LAMBDA(const int i)
  {
    NGP_EXPECT_TRUE(t.is_valid());
    NGP_EXPECT_TRUE(t.has_homogeneous_faces());
    NGP_EXPECT_FALSE(t.is_shell());

    NGP_EXPECT_EQ(t.rank(), stk::topology::ELEMENT_RANK);
    NGP_EXPECT_EQ(t.side_rank(), stk::topology::FACE_RANK);
    NGP_EXPECT_EQ(t.num_sides(), 4u);

    NGP_EXPECT_EQ(t.num_nodes(), 10u);
    NGP_EXPECT_EQ(t.num_vertices(), 4u);
    NGP_EXPECT_EQ(t.num_edges(), 6u);
    NGP_EXPECT_EQ(t.num_faces(), 4u);

    NGP_EXPECT_FALSE(t.defined_on_spatial_dimension(1));
    NGP_EXPECT_FALSE(t.defined_on_spatial_dimension(2));
    NGP_EXPECT_TRUE(t.defined_on_spatial_dimension(3));

    NGP_EXPECT_EQ(t.base(), stk::topology::TET_4);

    NGP_EXPECT_EQ(t.face_topology(0), stk::topology::TRI_6);
    NGP_EXPECT_EQ(t.face_topology(1), stk::topology::TRI_6);
    NGP_EXPECT_EQ(t.face_topology(2), stk::topology::TRI_6);
    NGP_EXPECT_EQ(t.face_topology(3), stk::topology::TRI_6);

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

NGP_TEST(stk_topology_ngp, tet_10)
{
  check_tet_10_on_device();
}


std::vector<std::vector<uint8_t>> get_gold_edge_node_ordinals_tet11() {
  return std::vector<std::vector<uint8_t>> {
    {0, 1,  4},
    {1, 2,  5},
    {2, 0,  6},
    {0, 3,  7},
    {1, 3,  8},
    {2, 3,  9}
  };
}

std::vector<std::vector<uint8_t>> get_gold_face_node_ordinals_tet11() {
  return std::vector<std::vector<uint8_t>> {
    {0, 1, 3,  4, 8, 7},
    {1, 2, 3,  5, 9, 8},
    {0, 3, 2,  7, 9, 6},
    {0, 2, 1,  6, 5, 4}
  };
}

std::vector<std::vector<uint8_t>> get_gold_permutation_node_ordinals_tet11() {
  return std::vector<std::vector<uint8_t>> {
    {0, 1, 2, 3,  4, 5, 6, 7, 8, 9,  10},
    {1, 2, 0, 3,  5, 6, 4, 8, 9, 7,  10},
    {2, 0, 1, 3,  6, 4, 5, 9, 7, 8,  10},
    {0, 3, 1, 2,  7, 8, 4, 6, 9, 5,  10},
    {3, 1, 0, 2,  8, 4, 7, 9, 5, 6,  10},
    {1, 0, 3, 2,  4, 7, 8, 5, 6, 9,  10},
    {0, 2, 3, 1,  6, 9, 7, 4, 5, 8,  10},
    {2, 3, 0, 1,  9, 7, 6, 5, 8, 4,  10},
    {3, 0, 2, 1,  7, 6, 9, 8, 4, 5,  10},
    {1, 3, 2, 0,  8, 9, 5, 4, 7, 6,  10},
    {3, 2, 1, 0,  9, 5, 8, 7, 6, 4,  10},
    {2, 1, 3, 0,  5, 8, 9, 6, 4, 7,  10}
  };
}

TEST(stk_topology, tet_11)
{
  stk::topology t = stk::topology::TET_11;

  EXPECT_TRUE(t.is_valid());
  EXPECT_TRUE(t.has_homogeneous_faces());
  EXPECT_FALSE(t.is_shell());

  EXPECT_EQ(t.rank(), stk::topology::ELEMENT_RANK);
  EXPECT_EQ(t.side_rank(), stk::topology::FACE_RANK);
  EXPECT_EQ(t.num_sides(), 4u);

  EXPECT_EQ(t.num_nodes(), 11u);
  EXPECT_EQ(t.num_vertices(), 4u);
  EXPECT_EQ(t.num_edges(), 6u);
  EXPECT_EQ(t.num_faces(), 4u);

  EXPECT_FALSE(t.defined_on_spatial_dimension(1));
  EXPECT_FALSE(t.defined_on_spatial_dimension(2));
  EXPECT_TRUE(t.defined_on_spatial_dimension(3));

  EXPECT_EQ(t.base(), stk::topology::TET_4);

  EXPECT_EQ(t.face_topology(0), stk::topology::TRI_6);
  EXPECT_EQ(t.face_topology(1), stk::topology::TRI_6);
  EXPECT_EQ(t.face_topology(2), stk::topology::TRI_6);
  EXPECT_EQ(t.face_topology(3), stk::topology::TRI_6);

  check_edge_node_ordinals(t, get_gold_edge_node_ordinals_tet11());
  check_edge_nodes(t, get_gold_edge_node_ordinals_tet11());

  check_side_node_ordinals(t, get_gold_face_node_ordinals_tet11());
  check_face_node_ordinals(t, get_gold_face_node_ordinals_tet11());
  check_side_nodes(t, get_gold_face_node_ordinals_tet11());
  check_face_nodes(t, get_gold_face_node_ordinals_tet11());

  check_permutation_node_ordinals(t, get_gold_permutation_node_ordinals_tet11());
  check_permutation_nodes(t, get_gold_permutation_node_ordinals_tet11());

  check_equivalent(t, get_gold_permutation_node_ordinals_tet11());
  check_lexicographical_smallest_permutation(t, get_gold_permutation_node_ordinals_tet11());
}

void check_tet_11_on_device()
{
  OrdinalType goldEdgeNodeOrdinals = fillGoldOrdinals(get_gold_edge_node_ordinals_tet11());
  OrdinalType goldFaceNodeOrdinals = fillGoldOrdinals(get_gold_face_node_ordinals_tet11());
  OrdinalType goldPermutationNodeOrdinals = fillGoldOrdinals(get_gold_permutation_node_ordinals_tet11());

  stk::topology t = stk::topology::TET_11;
  constexpr unsigned numNodes = stk::topology_detail::topology_data<stk::topology::TET_11>::num_nodes;
  EXPECT_EQ(11u, numNodes);

  Kokkos::parallel_for(stk::ngp::DeviceRangePolicy(0, 1), KOKKOS_LAMBDA(const int i)
  {
    NGP_EXPECT_TRUE(t.is_valid());
    NGP_EXPECT_TRUE(t.has_homogeneous_faces());
    NGP_EXPECT_FALSE(t.is_shell());

    NGP_EXPECT_EQ(t.rank(), stk::topology::ELEMENT_RANK);
    NGP_EXPECT_EQ(t.side_rank(), stk::topology::FACE_RANK);
    NGP_EXPECT_EQ(t.num_sides(), 4u);

    NGP_EXPECT_EQ(t.num_nodes(), 11u);
    NGP_EXPECT_EQ(t.num_vertices(), 4u);
    NGP_EXPECT_EQ(t.num_edges(), 6u);
    NGP_EXPECT_EQ(t.num_faces(), 4u);

    NGP_EXPECT_FALSE(t.defined_on_spatial_dimension(1));
    NGP_EXPECT_FALSE(t.defined_on_spatial_dimension(2));
    NGP_EXPECT_TRUE(t.defined_on_spatial_dimension(3));

    NGP_EXPECT_EQ(t.base(), stk::topology::TET_4);

    NGP_EXPECT_EQ(t.face_topology(0), stk::topology::TRI_6);
    NGP_EXPECT_EQ(t.face_topology(1), stk::topology::TRI_6);
    NGP_EXPECT_EQ(t.face_topology(2), stk::topology::TRI_6);
    NGP_EXPECT_EQ(t.face_topology(3), stk::topology::TRI_6);

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

NGP_TEST(stk_topology_ngp, tet_11)
{
  check_tet_11_on_device();
}

}
