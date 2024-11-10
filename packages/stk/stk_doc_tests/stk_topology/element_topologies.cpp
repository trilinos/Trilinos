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

#include "gtest/gtest.h"              // for AssertionResult, Message, TestPartResult, EXPECT_TRUE
#include "stk_topology/topology.hpp"  // for topology, topology::FACE_RANK, topology::EDGE_RANK
#include <vector>                     // for vector

namespace {

//ZeroDimElement
TEST(stk_topology_understanding, zero_dim_element)
{
  stk::topology sphere = stk::topology::PARTICLE;

  EXPECT_TRUE(sphere.is_valid());
  EXPECT_FALSE(sphere.has_homogeneous_faces());
  EXPECT_FALSE(sphere.is_shell());

  EXPECT_TRUE(sphere.rank() != stk::topology::NODE_RANK);
  EXPECT_TRUE(sphere.rank() != stk::topology::EDGE_RANK);
  EXPECT_TRUE(sphere.rank() != stk::topology::FACE_RANK);
  EXPECT_TRUE(sphere.rank() != stk::topology::CONSTRAINT_RANK);
  EXPECT_TRUE(sphere.rank() == stk::topology::ELEMENT_RANK);

  EXPECT_EQ(sphere.side_rank(), stk::topology::INVALID_RANK);

  EXPECT_EQ(sphere.dimension(),1u);
  EXPECT_EQ(sphere.num_nodes(),1u);
  EXPECT_EQ(sphere.num_vertices(),1u);
  EXPECT_EQ(sphere.num_edges(),0u);
  EXPECT_EQ(sphere.num_faces(),0u);
  EXPECT_EQ(sphere.num_sides(),0u);
  EXPECT_EQ(sphere.num_permutations(),1u);
  EXPECT_EQ(sphere.num_positive_permutations(),1u);

  EXPECT_FALSE(sphere.defined_on_spatial_dimension(0));

  EXPECT_TRUE(sphere.defined_on_spatial_dimension(1));
  EXPECT_TRUE(sphere.defined_on_spatial_dimension(2));
  EXPECT_TRUE(sphere.defined_on_spatial_dimension(3));

  EXPECT_EQ(sphere.base(),stk::topology::PARTICLE);
}

//ConstraintRank
TEST(stk_topology_understanding, constraint_ranks)
{
  // that's it folks, there is no topology that provides a constraint rank
}

void checkForValidOffsets(const unsigned numNodes, const unsigned *offsets, const unsigned *expectedNodeOffsets)
{
  for (unsigned i=0;i<numNodes;i++)
  {
    EXPECT_EQ(expectedNodeOffsets[i], offsets[i]);
  }
}

void checkPermutedNodeIds(const unsigned numNodes, const unsigned *offsets, const unsigned *nodeIds, const unsigned *goldIds)
{
  for (unsigned i=0;i<numNodes;i++)
  {
    EXPECT_EQ(goldIds[offsets[i]], nodeIds[i]);
  }
}

void checkNodeOrderingAndOffsetsForFaces(const stk::topology &element, const unsigned *elementNodes, const unsigned *expectedNodeOffsets)
{
  ASSERT_TRUE(element.has_homogeneous_faces());
  unsigned numNodesPerFace = element.face_topology(0).num_nodes();

  unsigned *offsets = new unsigned[numNodesPerFace];
  unsigned *nodeIds = new unsigned[numNodesPerFace];

  for (unsigned index=0; index<element.num_faces();index++)
  {
    element.face_node_ordinals(index, offsets);
    element.face_nodes(elementNodes, index, nodeIds);

    checkForValidOffsets(numNodesPerFace, offsets, expectedNodeOffsets + numNodesPerFace * index);
    checkPermutedNodeIds(numNodesPerFace, offsets, nodeIds, elementNodes);
  }

  delete [] nodeIds; nodeIds = 0;
  delete [] offsets; offsets = 0;
}

void checkNodeOrderingAndOffsetsForEdges(const stk::topology &element, const unsigned *elementNodes, const unsigned *expectedNodeOffsets)
{
  std::vector<unsigned> offsets;
  std::vector<unsigned> nodeIds;
  for (unsigned index=0; index<element.num_edges();index++) {
    unsigned numNodesPerEdge = element.edge_topology(index).num_nodes();
    offsets.resize(numNodesPerEdge);
    nodeIds.resize(numNodesPerEdge);

    element.edge_node_ordinals(index, offsets.data());
    element.edge_nodes(elementNodes, index, nodeIds.data());

    checkForValidOffsets(numNodesPerEdge, offsets.data(), expectedNodeOffsets + numNodesPerEdge * index);
    checkPermutedNodeIds(numNodesPerEdge, offsets.data(), nodeIds.data(), elementNodes);
  }
}

void checkNodeOrderingAndOffsetsForPermutations(const stk::topology &element, const unsigned *elementNodes, const unsigned *expectedNodeOffsets)
{
  unsigned numNodes = element.num_nodes();
  unsigned *offsets = new unsigned[numNodes];
  unsigned *nodeIds = new unsigned[numNodes];

  for (unsigned index=0; index<element.num_permutations();index++)
  {
    element.permutation_node_ordinals(index, offsets);
    element.permutation_nodes(elementNodes, index, nodeIds);

    checkForValidOffsets(numNodes, offsets, expectedNodeOffsets + numNodes * index);
    checkPermutedNodeIds(numNodes, offsets, nodeIds, elementNodes);
  }

  delete [] nodeIds; nodeIds = 0;
  delete [] offsets; offsets = 0;
}

//OneDimElement
TEST(stk_topology_understanding, one_dim_higher_order_element)
{
  stk::topology secondOrderBeam = stk::topology::BEAM_3;

  EXPECT_TRUE(secondOrderBeam.is_valid());
  EXPECT_FALSE(secondOrderBeam.has_homogeneous_faces());
  EXPECT_FALSE(secondOrderBeam.is_shell());

  EXPECT_TRUE(secondOrderBeam.rank() != stk::topology::NODE_RANK);
  EXPECT_TRUE(secondOrderBeam.rank() != stk::topology::EDGE_RANK);
  EXPECT_TRUE(secondOrderBeam.rank() != stk::topology::FACE_RANK);
  EXPECT_TRUE(secondOrderBeam.rank() != stk::topology::CONSTRAINT_RANK);
  EXPECT_TRUE(secondOrderBeam.rank() == stk::topology::ELEMENT_RANK);

  EXPECT_TRUE(secondOrderBeam.side_rank() == stk::topology::EDGE_RANK);

  EXPECT_EQ(2u, secondOrderBeam.dimension());
  EXPECT_EQ(3u, secondOrderBeam.num_nodes());
  EXPECT_EQ(2u, secondOrderBeam.num_vertices());
  EXPECT_EQ(1u, secondOrderBeam.num_edges());

  EXPECT_EQ(0u, secondOrderBeam.num_faces());
  EXPECT_EQ(1u, secondOrderBeam.num_positive_permutations());
  EXPECT_EQ(2u, secondOrderBeam.num_permutations());

  EXPECT_FALSE(secondOrderBeam.defined_on_spatial_dimension(0));
  EXPECT_FALSE(secondOrderBeam.defined_on_spatial_dimension(1));

  EXPECT_TRUE(secondOrderBeam.defined_on_spatial_dimension(2));
  EXPECT_TRUE(secondOrderBeam.defined_on_spatial_dimension(3));

  EXPECT_TRUE(secondOrderBeam.base() == stk::topology::BEAM_2);

  unsigned beamNodes[3] = { 10, 20, 14 }; // 10 *-------*-------* 20
  //            14
  {
    unsigned expectedNodeOffsets[3] = { 0, 1, 2 };
    //unit-test checking utility:
    checkNodeOrderingAndOffsetsForEdges(secondOrderBeam, beamNodes, expectedNodeOffsets);
  }

  {
    unsigned expectedNodeOffsets[6] = {
      0, 1, 2,
      1, 0, 2
    };

    //unit-test checking utility:
    checkNodeOrderingAndOffsetsForPermutations(secondOrderBeam, beamNodes, expectedNodeOffsets);
  }
}

//TwoDimElement
TEST(stk_topology_understanding, two_dim_higher_order_element)
{
  stk::topology secondOrderTriShell = stk::topology::SHELL_TRIANGLE_6;
  EXPECT_TRUE(secondOrderTriShell == stk::topology::SHELL_TRI_6);

  EXPECT_TRUE(secondOrderTriShell.is_valid());
  EXPECT_TRUE(secondOrderTriShell.has_homogeneous_faces());
  EXPECT_TRUE(secondOrderTriShell.is_shell());

  EXPECT_TRUE(secondOrderTriShell.rank() != stk::topology::NODE_RANK);
  EXPECT_TRUE(secondOrderTriShell.rank() != stk::topology::EDGE_RANK);
  EXPECT_TRUE(secondOrderTriShell.rank() != stk::topology::FACE_RANK);
  EXPECT_TRUE(secondOrderTriShell.rank() != stk::topology::CONSTRAINT_RANK);
  EXPECT_TRUE(secondOrderTriShell.rank() == stk::topology::ELEMENT_RANK);

  EXPECT_TRUE(secondOrderTriShell.side_rank() == stk::topology::FACE_RANK);

  EXPECT_EQ(3u, secondOrderTriShell.dimension());
  EXPECT_EQ(6u, secondOrderTriShell.num_nodes());
  EXPECT_EQ(3u, secondOrderTriShell.num_vertices());
  EXPECT_EQ(3u, secondOrderTriShell.num_edges());
  EXPECT_EQ(2u, secondOrderTriShell.num_faces());

  // permutations are the number of ways the number of vertices can be permuted
  EXPECT_EQ(6u, secondOrderTriShell.num_permutations());
  // positive permutations are ones that the normal is maintained
  EXPECT_EQ(3u, secondOrderTriShell.num_positive_permutations());

  EXPECT_FALSE(secondOrderTriShell.defined_on_spatial_dimension(0));
  EXPECT_FALSE(secondOrderTriShell.defined_on_spatial_dimension(1));
  EXPECT_FALSE(secondOrderTriShell.defined_on_spatial_dimension(2));

  EXPECT_TRUE(secondOrderTriShell.defined_on_spatial_dimension(3));

  EXPECT_TRUE(secondOrderTriShell.base() == stk::topology::SHELL_TRI_3);
  EXPECT_TRUE(secondOrderTriShell.base() == stk::topology::SHELL_TRIANGLE_3);

  unsigned shellNodes[6] = { 10, 11, 12, 100, 101, 102 }; // first 3 are vertex nodes (picture?)

  {
    unsigned goldValuesEdgeOffsets[9] = {
      0, 1, 3,
      1, 2, 4,
      2, 0, 5
    };

    //unit-test checking utility:
    checkNodeOrderingAndOffsetsForEdges(secondOrderTriShell, shellNodes, goldValuesEdgeOffsets);
  }

  {
    unsigned goldValuesFaceNodeOffsets[12] = {
      0, 1, 2, 3, 4, 5,
      0, 2, 1, 5, 4, 3
    };

    //unit-test checking utility:
    checkNodeOrderingAndOffsetsForFaces(secondOrderTriShell, shellNodes, goldValuesFaceNodeOffsets);
  }

  {
    unsigned goldValueOffsetsPerm[36] = {
      0, 1, 2, 3, 4, 5,
      2, 0, 1, 5, 3, 4,
      1, 2, 0, 4, 5, 3,
      0, 2, 1, 5, 4, 3,
      2, 1, 0, 4, 3, 5,
      1, 0, 2, 3, 5, 4
    };

    //unit-test checking utility:
    checkNodeOrderingAndOffsetsForPermutations(secondOrderTriShell, shellNodes, goldValueOffsetsPerm);
  }
}


//ThreeDimElement
TEST(stk_topology_understanding, three_dim_linear_element)
{
  stk::topology hex8 = stk::topology::HEX_8;
  EXPECT_TRUE(hex8 == stk::topology::HEXAHEDRON_8);

  EXPECT_TRUE(hex8.is_valid());
  EXPECT_TRUE(hex8.has_homogeneous_faces());
  EXPECT_FALSE(hex8.is_shell());

  EXPECT_TRUE(hex8.rank() != stk::topology::NODE_RANK);
  EXPECT_TRUE(hex8.rank() != stk::topology::EDGE_RANK);
  EXPECT_TRUE(hex8.rank() != stk::topology::FACE_RANK);
  EXPECT_TRUE(hex8.rank() != stk::topology::CONSTRAINT_RANK);
  EXPECT_TRUE(hex8.rank() == stk::topology::ELEMENT_RANK);

  EXPECT_TRUE(hex8.side_rank() == stk::topology::FACE_RANK);

  EXPECT_EQ(3u, hex8.dimension());
  EXPECT_EQ(8u, hex8.num_nodes());
  EXPECT_EQ(8u, hex8.num_vertices());
  EXPECT_EQ(12u, hex8.num_edges());
  EXPECT_EQ(6u, hex8.num_faces());

  if (stk::topology::topology_type<stk::topology::HEX_8>::num_permutations > 1) {
    // permutations are the number of ways the number of vertices can be permuted
    EXPECT_EQ(24u, hex8.num_permutations());
    // positive permutations are ones that the normal is maintained
    EXPECT_EQ(24u, hex8.num_positive_permutations());
  }

  EXPECT_FALSE(hex8.defined_on_spatial_dimension(0));
  EXPECT_FALSE(hex8.defined_on_spatial_dimension(1));
  EXPECT_FALSE(hex8.defined_on_spatial_dimension(2));

  EXPECT_TRUE(hex8.defined_on_spatial_dimension(3));

  EXPECT_TRUE(hex8.base() == stk::topology::HEX_8);

  unsigned hex8Nodes[8] = { 0, 1, 2, 3, 4, 5, 6, 7 };

  {
    for(unsigned i = 0; i < hex8.num_edges(); i++) {
      EXPECT_EQ(hex8.edge_topology(i), stk::topology::LINE_2);
      ASSERT_EQ(hex8.edge_topology(i).num_nodes(), 2u);
    }

    unsigned goldValuesEdgeOffsets[24] = {
      0, 1,
      1, 2,
      2, 3,
      3, 0,
      4, 5,
      5, 6,
      6, 7,
      7, 4,
      0, 4,
      1, 5,
      2, 6,
      3, 7
    };

    //unit-test checking utility:
    checkNodeOrderingAndOffsetsForEdges(hex8, hex8Nodes, goldValuesEdgeOffsets);
  }

  {
    stk::topology goldFaceTopology = stk::topology::QUAD_4;
    unsigned goldNumNodesPerFace = 4;
    for (unsigned faceIndex=0;faceIndex<hex8.num_faces();faceIndex++)
    {
      EXPECT_EQ(goldFaceTopology, hex8.face_topology(faceIndex));
      ASSERT_EQ(goldNumNodesPerFace, hex8.face_topology(faceIndex).num_nodes());
    }

    unsigned goldValuesFaceOffsets[24] = {
      0, 1, 5, 4,
      1, 2, 6, 5,
      2, 3, 7, 6,
      0, 4, 7, 3,
      0, 3, 2, 1,
      4, 5, 6, 7
    };

    //unit-test checking utility:
    checkNodeOrderingAndOffsetsForFaces(hex8, hex8Nodes, goldValuesFaceOffsets);
  }

}
//EndThreeDim

}

