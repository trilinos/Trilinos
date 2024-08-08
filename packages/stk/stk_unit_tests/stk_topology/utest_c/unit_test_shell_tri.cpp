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

std::vector<std::vector<uint8_t>> get_gold_edge_node_ordinals_shell_tri3() {
  return std::vector<std::vector<uint8_t>> {
    {0, 1},
    {1, 2},
    {2, 0}
  };
}

std::vector<std::vector<uint8_t>> get_gold_face_node_ordinals_shell_tri3() {
  return std::vector<std::vector<uint8_t>> {
    {0, 1, 2},
    {0, 2, 1}
  };
}

std::vector<std::vector<uint8_t>> get_gold_side_node_ordinals_shell_tri3() {
  return std::vector<std::vector<uint8_t>> {
    {0, 1, 2},
    {0, 2, 1},
    {0, 1},
    {1, 2},
    {2, 0}
  };
}

TEST(stk_topology, shell_tri_3)
{
  stk::topology t = stk::topology::SHELL_TRI_3;

  EXPECT_TRUE(t.is_valid());
  EXPECT_TRUE(t.has_homogeneous_faces());
  EXPECT_TRUE(t.is_shell());

  EXPECT_EQ(t.rank(),stk::topology::ELEMENT_RANK);
  EXPECT_EQ(t.side_rank(),stk::topology::FACE_RANK);

  EXPECT_EQ(t.num_nodes(),3u);
  EXPECT_EQ(t.num_vertices(),3u);
  EXPECT_EQ(t.num_edges(),3u);
  EXPECT_EQ(t.num_faces(),2u);
  EXPECT_EQ(t.num_sides(),5u);

  EXPECT_FALSE(t.defined_on_spatial_dimension(1));
  EXPECT_FALSE(t.defined_on_spatial_dimension(2));
  EXPECT_TRUE(t.defined_on_spatial_dimension(3));

  EXPECT_EQ(t.base(),stk::topology::SHELL_TRI_3);

  EXPECT_EQ(t.edge_topology(0), stk::topology::LINE_2);
  EXPECT_EQ(t.edge_topology(1), stk::topology::LINE_2);
  EXPECT_EQ(t.edge_topology(2), stk::topology::LINE_2);

  EXPECT_EQ(t.face_topology(0), stk::topology::TRI_3);
  EXPECT_EQ(t.face_topology(1), stk::topology::TRI_3);

  EXPECT_EQ(t.side_topology(0), stk::topology::TRI_3);
  EXPECT_EQ(t.side_topology(1), stk::topology::TRI_3);
  EXPECT_EQ(t.side_topology(2), stk::topology::SHELL_SIDE_BEAM_2);
  EXPECT_EQ(t.side_topology(3), stk::topology::SHELL_SIDE_BEAM_2);
  EXPECT_EQ(t.side_topology(4), stk::topology::SHELL_SIDE_BEAM_2);

  check_edge_node_ordinals(t, get_gold_edge_node_ordinals_shell_tri3());
  check_edge_nodes(t, get_gold_edge_node_ordinals_shell_tri3());

  check_face_node_ordinals(t, get_gold_face_node_ordinals_shell_tri3());
  check_face_nodes(t, get_gold_face_node_ordinals_shell_tri3());

  check_side_node_ordinals(t, get_gold_side_node_ordinals_shell_tri3());
  check_side_nodes(t, get_gold_side_node_ordinals_shell_tri3());
}

void check_shell_tri_3_on_device()
{
  OrdinalType goldEdgeNodeOrdinals = fillGoldOrdinals(get_gold_edge_node_ordinals_shell_tri3());
  OrdinalType goldFaceNodeOrdinals = fillGoldOrdinals(get_gold_face_node_ordinals_shell_tri3());
  OrdinalType goldSideNodeOrdinals = fillGoldOrdinals(get_gold_side_node_ordinals_shell_tri3());

  stk::topology t = stk::topology::SHELL_TRI_3;
  constexpr unsigned numNodes = stk::topology_detail::topology_data<stk::topology::SHELL_TRI_3>::num_nodes;
  EXPECT_EQ(3u, numNodes);

  Kokkos::parallel_for(stk::ngp::DeviceRangePolicy(0, 1), KOKKOS_LAMBDA(const int i)
  {
    NGP_EXPECT_TRUE(t.is_valid());
    NGP_EXPECT_TRUE(t.has_homogeneous_faces());
    NGP_EXPECT_TRUE(t.is_shell());

    NGP_EXPECT_EQ(t.rank(),stk::topology::ELEMENT_RANK);
    NGP_EXPECT_EQ(t.side_rank(),stk::topology::FACE_RANK);

    NGP_EXPECT_EQ(t.num_nodes(),3u);
    NGP_EXPECT_EQ(t.num_vertices(),3u);
    NGP_EXPECT_EQ(t.num_edges(),3u);
    NGP_EXPECT_EQ(t.num_faces(),2u);
    NGP_EXPECT_EQ(t.num_sides(),5u);

    NGP_EXPECT_FALSE(t.defined_on_spatial_dimension(1));
    NGP_EXPECT_FALSE(t.defined_on_spatial_dimension(2));
    NGP_EXPECT_TRUE(t.defined_on_spatial_dimension(3));

    NGP_EXPECT_EQ(t.base(),stk::topology::SHELL_TRI_3);

    NGP_EXPECT_EQ(t.edge_topology(0), stk::topology::LINE_2);
    NGP_EXPECT_EQ(t.edge_topology(1), stk::topology::LINE_2);
    NGP_EXPECT_EQ(t.edge_topology(2), stk::topology::LINE_2);
  
    NGP_EXPECT_EQ(t.face_topology(0), stk::topology::TRI_3);
    NGP_EXPECT_EQ(t.face_topology(1), stk::topology::TRI_3);
  
    NGP_EXPECT_EQ(t.side_topology(0), stk::topology::TRI_3);
    NGP_EXPECT_EQ(t.side_topology(1), stk::topology::TRI_3);
    NGP_EXPECT_EQ(t.side_topology(2), stk::topology::SHELL_SIDE_BEAM_2);
    NGP_EXPECT_EQ(t.side_topology(3), stk::topology::SHELL_SIDE_BEAM_2);
    NGP_EXPECT_EQ(t.side_topology(4), stk::topology::SHELL_SIDE_BEAM_2);

    check_edge_node_ordinals_ngp<numNodes>(t, goldEdgeNodeOrdinals);
    check_edge_nodes_ngp<numNodes>(t, goldEdgeNodeOrdinals);
  });

  Kokkos::parallel_for(stk::ngp::DeviceRangePolicy(0, 1), KOKKOS_LAMBDA(const int i)
  {
    check_face_node_ordinals_ngp<numNodes>(t, goldFaceNodeOrdinals);
    check_face_nodes_ngp<numNodes>(t, goldFaceNodeOrdinals);

    check_side_node_ordinals_ngp<numNodes>(t, goldSideNodeOrdinals);
    check_side_nodes_ngp<numNodes>(t, goldSideNodeOrdinals);
  });
}

NGP_TEST(stk_topology_ngp, shell_tri_3)
{
  check_shell_tri_3_on_device();
}

std::vector<std::vector<uint8_t>> get_gold_edge_node_ordinals_shell_tri4() {
  return std::vector<std::vector<uint8_t>> {
    {0, 1},
    {1, 2},
    {2, 0}
  };
}

std::vector<std::vector<uint8_t>> get_gold_face_node_ordinals_shell_tri4() {
  return std::vector<std::vector<uint8_t>> {
    {0, 1, 2, 3},
    {0, 2, 1, 3}
  };
}

std::vector<std::vector<uint8_t>> get_gold_side_node_ordinals_shell_tri4() {
  return std::vector<std::vector<uint8_t>> {
    {0, 1, 2, 3},
    {0, 2, 1, 3},
    {0, 1},
    {1, 2},
    {2, 0}
  };
}

TEST(stk_topology, shell_tri_4)
{
  stk::topology t = stk::topology::SHELL_TRI_4;

  EXPECT_TRUE(t.is_valid());
  EXPECT_TRUE(t.has_homogeneous_faces());
  EXPECT_TRUE(t.is_shell());

  EXPECT_EQ(t.rank(),stk::topology::ELEMENT_RANK);
  EXPECT_EQ(t.side_rank(),stk::topology::FACE_RANK);

  EXPECT_EQ(t.num_nodes(),4u);
  EXPECT_EQ(t.num_vertices(),3u);
  EXPECT_EQ(t.num_faces(),2u);
  EXPECT_EQ(t.num_edges(),3u);
  EXPECT_EQ(t.num_sides(),5u);

  EXPECT_FALSE(t.defined_on_spatial_dimension(1));
  EXPECT_FALSE(t.defined_on_spatial_dimension(2));
  EXPECT_TRUE(t.defined_on_spatial_dimension(3));

  EXPECT_EQ(t.base(),stk::topology::SHELL_TRI_3);

  EXPECT_EQ(t.edge_topology(0), stk::topology::LINE_2);
  EXPECT_EQ(t.edge_topology(1), stk::topology::LINE_2);
  EXPECT_EQ(t.edge_topology(2), stk::topology::LINE_2);

  EXPECT_EQ(t.face_topology(0), stk::topology::TRI_4);
  EXPECT_EQ(t.face_topology(1), stk::topology::TRI_4);

  EXPECT_EQ(t.side_topology(0), stk::topology::TRI_4);
  EXPECT_EQ(t.side_topology(1), stk::topology::TRI_4);
  EXPECT_EQ(t.side_topology(2), stk::topology::SHELL_SIDE_BEAM_2);
  EXPECT_EQ(t.side_topology(3), stk::topology::SHELL_SIDE_BEAM_2);
  EXPECT_EQ(t.side_topology(4), stk::topology::SHELL_SIDE_BEAM_2);

  check_edge_node_ordinals(t, get_gold_edge_node_ordinals_shell_tri4());
  check_edge_nodes(t, get_gold_edge_node_ordinals_shell_tri4());

  check_face_node_ordinals(t, get_gold_face_node_ordinals_shell_tri4());
  check_face_nodes(t, get_gold_face_node_ordinals_shell_tri4());

  check_side_node_ordinals(t, get_gold_side_node_ordinals_shell_tri4());
  check_side_nodes(t, get_gold_side_node_ordinals_shell_tri4());
}

void check_shell_tri_4_on_device()
{
  OrdinalType goldEdgeNodeOrdinals = fillGoldOrdinals(get_gold_edge_node_ordinals_shell_tri4());
  OrdinalType goldFaceNodeOrdinals = fillGoldOrdinals(get_gold_face_node_ordinals_shell_tri4());
  OrdinalType goldSideNodeOrdinals = fillGoldOrdinals(get_gold_side_node_ordinals_shell_tri4());

  stk::topology t = stk::topology::SHELL_TRI_4;
  constexpr unsigned numNodes = stk::topology_detail::topology_data<stk::topology::SHELL_TRI_4>::num_nodes;
  EXPECT_EQ(4u, numNodes);

  Kokkos::parallel_for(stk::ngp::DeviceRangePolicy(0, 1), KOKKOS_LAMBDA(const int i)
  {
    NGP_EXPECT_TRUE(t.is_valid());
    NGP_EXPECT_TRUE(t.has_homogeneous_faces());
    NGP_EXPECT_TRUE(t.is_shell());

    NGP_EXPECT_EQ(t.rank(),stk::topology::ELEMENT_RANK);
    NGP_EXPECT_EQ(t.side_rank(),stk::topology::FACE_RANK);

    NGP_EXPECT_EQ(t.num_nodes(),4u);
    NGP_EXPECT_EQ(t.num_vertices(),3u);
    NGP_EXPECT_EQ(t.num_edges(),3u);
    NGP_EXPECT_EQ(t.num_faces(),2u);
    NGP_EXPECT_EQ(t.num_sides(),5u);

    NGP_EXPECT_FALSE(t.defined_on_spatial_dimension(1));
    NGP_EXPECT_FALSE(t.defined_on_spatial_dimension(2));
    NGP_EXPECT_TRUE(t.defined_on_spatial_dimension(3));

    NGP_EXPECT_EQ(t.base(),stk::topology::SHELL_TRI_3);

    NGP_EXPECT_EQ(t.edge_topology(0), stk::topology::LINE_2);
    NGP_EXPECT_EQ(t.edge_topology(1), stk::topology::LINE_2);
    NGP_EXPECT_EQ(t.edge_topology(2), stk::topology::LINE_2);
  
    NGP_EXPECT_EQ(t.face_topology(0), stk::topology::TRI_4);
    NGP_EXPECT_EQ(t.face_topology(1), stk::topology::TRI_4);
  
    NGP_EXPECT_EQ(t.side_topology(0), stk::topology::TRI_4);
    NGP_EXPECT_EQ(t.side_topology(1), stk::topology::TRI_4);
    NGP_EXPECT_EQ(t.side_topology(2), stk::topology::SHELL_SIDE_BEAM_2);
    NGP_EXPECT_EQ(t.side_topology(3), stk::topology::SHELL_SIDE_BEAM_2);
    NGP_EXPECT_EQ(t.side_topology(4), stk::topology::SHELL_SIDE_BEAM_2);

    check_edge_node_ordinals_ngp<numNodes>(t, goldEdgeNodeOrdinals);
    check_edge_nodes_ngp<numNodes>(t, goldEdgeNodeOrdinals);
  });

  Kokkos::parallel_for(stk::ngp::DeviceRangePolicy(0, 1), KOKKOS_LAMBDA(const int i)
  {
    check_face_node_ordinals_ngp<numNodes>(t, goldFaceNodeOrdinals);
    check_face_nodes_ngp<numNodes>(t, goldFaceNodeOrdinals);

    check_side_node_ordinals_ngp<numNodes>(t, goldSideNodeOrdinals);
    check_side_nodes_ngp<numNodes>(t, goldSideNodeOrdinals);
  });
}

NGP_TEST(stk_topology_ngp, shell_tri_4)
{
  check_shell_tri_4_on_device();
}

std::vector<std::vector<uint8_t>> get_gold_edge_node_ordinals_shell_tri6() {
  return std::vector<std::vector<uint8_t>> {
    {0, 1, 3},
    {1, 2, 4},
    {2, 0, 5}
  };
}

std::vector<std::vector<uint8_t>> get_gold_face_node_ordinals_shell_tri6() {
  return std::vector<std::vector<uint8_t>> {
    {0, 1, 2,  3, 4, 5},
    {0, 2, 1,  5, 4, 3}
  };
}

std::vector<std::vector<uint8_t>> get_gold_side_node_ordinals_shell_tri6() {
  return std::vector<std::vector<uint8_t>> {
    {0, 1, 2,  3, 4, 5},
    {0, 2, 1,  5, 4, 3},
    {0, 1, 3},
    {1, 2, 4},
    {2, 0, 5}
  };
}

TEST(stk_topology, shell_tri_6)
{
  stk::topology t = stk::topology::SHELL_TRI_6;

  EXPECT_TRUE(t.is_valid());
  EXPECT_TRUE(t.has_homogeneous_faces());
  EXPECT_TRUE(t.is_shell());

  EXPECT_EQ(t.rank(),stk::topology::ELEMENT_RANK);
  EXPECT_EQ(t.side_rank(),stk::topology::FACE_RANK);

  EXPECT_EQ(t.num_nodes(),6u);
  EXPECT_EQ(t.num_vertices(),3u);
  EXPECT_EQ(t.num_faces(),2u);
  EXPECT_EQ(t.num_edges(),3u);
  EXPECT_EQ(t.num_sides(),5u);

  EXPECT_FALSE(t.defined_on_spatial_dimension(1));
  EXPECT_FALSE(t.defined_on_spatial_dimension(2));
  EXPECT_TRUE(t.defined_on_spatial_dimension(3));

  EXPECT_EQ(t.base(),stk::topology::SHELL_TRI_3);

  EXPECT_EQ(t.edge_topology(0), stk::topology::LINE_3);
  EXPECT_EQ(t.edge_topology(1), stk::topology::LINE_3);
  EXPECT_EQ(t.edge_topology(2), stk::topology::LINE_3);

  EXPECT_EQ(t.face_topology(0), stk::topology::TRI_6);
  EXPECT_EQ(t.face_topology(1), stk::topology::TRI_6);

  EXPECT_EQ(t.side_topology(0), stk::topology::TRI_6);
  EXPECT_EQ(t.side_topology(1), stk::topology::TRI_6);
  EXPECT_EQ(t.side_topology(2), stk::topology::SHELL_SIDE_BEAM_3);
  EXPECT_EQ(t.side_topology(3), stk::topology::SHELL_SIDE_BEAM_3);
  EXPECT_EQ(t.side_topology(4), stk::topology::SHELL_SIDE_BEAM_3);

  check_edge_node_ordinals(t, get_gold_edge_node_ordinals_shell_tri6());
  check_edge_nodes(t, get_gold_edge_node_ordinals_shell_tri6());

  check_face_node_ordinals(t, get_gold_face_node_ordinals_shell_tri6());
  check_face_nodes(t, get_gold_face_node_ordinals_shell_tri6());

  check_side_node_ordinals(t, get_gold_side_node_ordinals_shell_tri6());
  check_side_nodes(t, get_gold_side_node_ordinals_shell_tri6());
}

void check_shell_tri_6_on_device()
{
  OrdinalType goldEdgeNodeOrdinals = fillGoldOrdinals(get_gold_edge_node_ordinals_shell_tri6());
  OrdinalType goldFaceNodeOrdinals = fillGoldOrdinals(get_gold_face_node_ordinals_shell_tri6());
  OrdinalType goldSideNodeOrdinals = fillGoldOrdinals(get_gold_side_node_ordinals_shell_tri6());

  stk::topology t = stk::topology::SHELL_TRI_6;
  constexpr unsigned numNodes = stk::topology_detail::topology_data<stk::topology::SHELL_TRI_6>::num_nodes;
  EXPECT_EQ(6u, numNodes);

  Kokkos::parallel_for(stk::ngp::DeviceRangePolicy(0, 1), KOKKOS_LAMBDA(const int i)
  {
    NGP_EXPECT_TRUE(t.is_valid());
    NGP_EXPECT_TRUE(t.has_homogeneous_faces());
    NGP_EXPECT_TRUE(t.is_shell());

    NGP_EXPECT_EQ(t.rank(),stk::topology::ELEMENT_RANK);
    NGP_EXPECT_EQ(t.side_rank(),stk::topology::FACE_RANK);

    NGP_EXPECT_EQ(t.num_nodes(),6u);
    NGP_EXPECT_EQ(t.num_vertices(),3u);
    NGP_EXPECT_EQ(t.num_edges(),3u);
    NGP_EXPECT_EQ(t.num_faces(),2u);
    NGP_EXPECT_EQ(t.num_sides(),5u);

    NGP_EXPECT_FALSE(t.defined_on_spatial_dimension(1));
    NGP_EXPECT_FALSE(t.defined_on_spatial_dimension(2));
    NGP_EXPECT_TRUE(t.defined_on_spatial_dimension(3));

    NGP_EXPECT_EQ(t.base(),stk::topology::SHELL_TRI_3);

    NGP_EXPECT_EQ(t.edge_topology(0), stk::topology::LINE_3);
    NGP_EXPECT_EQ(t.edge_topology(1), stk::topology::LINE_3);
    NGP_EXPECT_EQ(t.edge_topology(2), stk::topology::LINE_3);
  
    NGP_EXPECT_EQ(t.face_topology(0), stk::topology::TRI_6);
    NGP_EXPECT_EQ(t.face_topology(1), stk::topology::TRI_6);

    NGP_EXPECT_EQ(t.shell_side_topology(0), stk::topology::SHELL_SIDE_BEAM_3);
    NGP_EXPECT_EQ(t.shell_side_topology(1), stk::topology::SHELL_SIDE_BEAM_3);
    NGP_EXPECT_EQ(t.shell_side_topology(2), stk::topology::SHELL_SIDE_BEAM_3);
  
    NGP_EXPECT_EQ(t.side_topology(0), stk::topology::TRI_6);
    NGP_EXPECT_EQ(t.side_topology(1), stk::topology::TRI_6);
    NGP_EXPECT_EQ(t.side_topology(2), stk::topology::SHELL_SIDE_BEAM_3);
    NGP_EXPECT_EQ(t.side_topology(3), stk::topology::SHELL_SIDE_BEAM_3);
    NGP_EXPECT_EQ(t.side_topology(4), stk::topology::SHELL_SIDE_BEAM_3);

    check_edge_node_ordinals_ngp<numNodes>(t, goldEdgeNodeOrdinals);
    check_edge_nodes_ngp<numNodes>(t, goldEdgeNodeOrdinals);
  });

  Kokkos::parallel_for(stk::ngp::DeviceRangePolicy(0, 1), KOKKOS_LAMBDA(const int i)
  {
    check_face_node_ordinals_ngp<numNodes>(t, goldFaceNodeOrdinals);
    check_face_nodes_ngp<numNodes>(t, goldFaceNodeOrdinals);

    check_side_node_ordinals_ngp<numNodes>(t, goldSideNodeOrdinals);
    check_side_nodes_ngp<numNodes>(t, goldSideNodeOrdinals);
  });
}

NGP_TEST(stk_topology_ngp, shell_tri_6)
{
  check_shell_tri_6_on_device();
}

}
