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

#include <gtest/gtest.h>
#include <stk_ngp_test/ngp_test.hpp>
#include <stk_topology/topology.hpp>
#include "stk_util/util/ReportHandler.hpp"
#include <limits>
#define INVALID UINT_MAX

inline void check_side_node_ordinals(stk::topology topology, std::vector<std::vector<unsigned>> & gold_side_node_ordinals)
{
  for (unsigned side = 0; side < topology.num_sides(); ++side) {
    stk::topology sideTopo = topology.side_topology(side);
    // If this is a topology where the "sides" are nodes, num_nodes() on the side topology
    // will come back as zero when there is actually just one node on the "side".
    unsigned numSideNodes = (sideTopo.num_nodes() > 0) ? sideTopo.num_nodes() : 1;
    std::vector<unsigned> side_node_ordinals(numSideNodes);
    topology.side_node_ordinals(side, side_node_ordinals.data());
    EXPECT_EQ(gold_side_node_ordinals[side], side_node_ordinals);
  }
}

inline void check_edge_node_ordinals(stk::topology topology, std::vector<std::vector<unsigned>> & gold_edge_node_ordinals)
{
  std::vector<unsigned> edge_node_ordinals;
  for (unsigned edge = 0; edge < topology.num_edges(); ++edge) {
    stk::topology edgeTopo = topology.edge_topology(edge);
    unsigned numEdgeNodes = edgeTopo.num_nodes();
    edge_node_ordinals.resize(numEdgeNodes);
    topology.edge_node_ordinals(edge, edge_node_ordinals.data());
    EXPECT_EQ(gold_edge_node_ordinals[edge], edge_node_ordinals);
  }
}

inline void check_face_node_ordinals(stk::topology topology, std::vector<std::vector<unsigned>> & gold_face_node_ordinals)
{
  for (unsigned face = 0; face < topology.num_faces(); ++face) {
    stk::topology faceTopo = topology.face_topology(face);
    unsigned numFaceNodes = faceTopo.num_nodes();
    std::vector<unsigned> face_node_ordinals(numFaceNodes);
    topology.face_node_ordinals(face, face_node_ordinals.data());
    EXPECT_EQ(gold_face_node_ordinals[face], face_node_ordinals);
  }
}

inline void check_permutation_node_ordinals(stk::topology topology, std::vector<std::vector<unsigned>> & gold_permutation_node_ordinals)
{
  const unsigned numPermutations = topology.num_permutations();
  for (unsigned permutationOrdinal = 0; permutationOrdinal < numPermutations; ++permutationOrdinal) {
    std::vector<unsigned> permutation_node_ordinals(topology.num_nodes());
    topology.permutation_node_ordinals(permutationOrdinal, permutation_node_ordinals.data());
    EXPECT_EQ(gold_permutation_node_ordinals[permutationOrdinal], permutation_node_ordinals);
  }
}

inline void fill_permutation(unsigned permutation, const std::vector<std::vector<unsigned>> & permutationNodeOrdinals, std::vector<unsigned> & permutedOrdinals)
{
  const unsigned numNodes = permutationNodeOrdinals[permutation].size();
  permutedOrdinals.resize(numNodes);
  for (unsigned node = 0; node < numNodes; ++node) {
    permutedOrdinals[node] = permutationNodeOrdinals[permutation][node];
  }
}

inline void check_equivalent(stk::topology topology, std::vector<std::vector<unsigned>> & gold_permutation_node_ordinals)
{
  const unsigned numPermutations = topology.num_permutations();
  std::vector<unsigned> baseNodeArray(topology.num_nodes());
  fill_permutation(0, gold_permutation_node_ordinals, baseNodeArray);

  for (unsigned permutationOrdinal = 0; permutationOrdinal < numPermutations; ++permutationOrdinal) {
    std::vector<unsigned> permutedNodeArray(topology.num_nodes());
    fill_permutation(permutationOrdinal, gold_permutation_node_ordinals, permutedNodeArray);
    stk::EquivalentPermutation result = topology.is_equivalent(baseNodeArray.data(), permutedNodeArray.data());
    EXPECT_TRUE(result.is_equivalent);
    EXPECT_EQ(permutationOrdinal, result.permutation_number);
  }

  // Throw in a few duds that don't match, just for good measure
  for (unsigned node = 0; node < topology.num_nodes(); ++node) {
    std::vector<unsigned> permutedNodeArray(baseNodeArray);
    permutedNodeArray[node] = topology.num_nodes();  // Tweak one node to not match
    stk::EquivalentPermutation result = topology.is_equivalent(baseNodeArray.data(), permutedNodeArray.data());
    EXPECT_FALSE(result.is_equivalent);
    EXPECT_EQ(0u, result.permutation_number);
  }
}

inline void check_lexicographical_smallest_permutation(stk::topology topology, std::vector<std::vector<unsigned>> & gold_permutation_node_ordinals)
{
  const unsigned numPermutations = topology.num_permutations();
  const unsigned numNodes = (topology.num_nodes() > 0) ? topology.num_nodes() : 1;
  std::vector<unsigned> nodeArray(numNodes);

  const unsigned firstPermutation = 0;
  const unsigned lastPermutation = (numPermutations > 0) ? (numPermutations - 1) : 0;

  for (unsigned permWithSmallestSorting = firstPermutation; permWithSmallestSorting <= lastPermutation; ++permWithSmallestSorting) {
    // Fill the nodes with IDs so that the target permutation will sort first
    for (unsigned node = 0; node < numNodes; ++node) {
      unsigned nodeIndex = gold_permutation_node_ordinals[permWithSmallestSorting][node];
      nodeArray[nodeIndex] = node + 1;
    }

    // Check if we retrieve this target permutation as the lexicographically-smallest
    unsigned result = topology.lexicographical_smallest_permutation(nodeArray.data());
    EXPECT_EQ(permWithSmallestSorting, result);

    // Check if we retrieve a positive permutation as the lexicographically-smallest when we set
    // the flag.  Can't predict which one we will get, so just check if it's positive.
    const bool onlyPositivePermutations = true;
    result = topology.lexicographical_smallest_permutation(nodeArray.data(), onlyPositivePermutations);
    if (permWithSmallestSorting < topology.num_positive_permutations() || (topology.num_positive_permutations() == 0)) {
      EXPECT_EQ(permWithSmallestSorting, result);  // Match exactly
    }
    else {
      EXPECT_LT(result, topology.num_positive_permutations());  // Confirm positive
    }
  }
}

template <unsigned MAX_NODES>
STK_INLINE_FUNCTION
void check_side_node_ordinals_ngp(stk::topology topology, unsigned gold_side_node_ordinals[][MAX_NODES])
{
  for (unsigned side = 0; side < topology.num_sides(); ++side) {
    stk::topology sideTopo = topology.side_topology(side);
    // If this is a topology where the "sides" are nodes, num_nodes() on the side topology
    // will come back as zero when there is actually just one node on the "side".
    unsigned numSideNodes = (sideTopo.num_nodes() > 0) ? sideTopo.num_nodes() : 1;
    unsigned side_node_ordinals[MAX_NODES];
    topology.side_node_ordinals(side, side_node_ordinals);
    for (unsigned side_node = 0; side_node < numSideNodes; ++side_node) {
      NGP_EXPECT_EQ(gold_side_node_ordinals[side][side_node], side_node_ordinals[side_node]);
    }
  }
}

template <unsigned MAX_NODES>
STK_INLINE_FUNCTION
void check_edge_node_ordinals_ngp(stk::topology topology, unsigned gold_edge_node_ordinals[][MAX_NODES])
{
  unsigned edge_node_ordinals[MAX_NODES];
  for (unsigned edge = 0; edge < topology.num_edges(); ++edge) {
    stk::topology edgeTopo = topology.edge_topology(edge);
    unsigned numEdgeNodes = edgeTopo.num_nodes();
    topology.edge_node_ordinals(edge, edge_node_ordinals);
    for (unsigned edge_node = 0; edge_node < numEdgeNodes; ++edge_node) {
      NGP_EXPECT_EQ(gold_edge_node_ordinals[edge][edge_node], edge_node_ordinals[edge_node]);
    }
  }
}

template <unsigned MAX_NODES>
STK_INLINE_FUNCTION
void check_face_node_ordinals_ngp(stk::topology topology, unsigned gold_face_node_ordinals[][MAX_NODES])
{
  for (unsigned face = 0; face < topology.num_faces(); ++face) {
    stk::topology faceTopo = topology.face_topology(face);
    unsigned numFaceNodes = faceTopo.num_nodes();
    unsigned face_node_ordinals[MAX_NODES];
    topology.face_node_ordinals(face, face_node_ordinals);
    for (unsigned face_node = 0; face_node < numFaceNodes; ++face_node) {
      NGP_EXPECT_EQ(gold_face_node_ordinals[face][face_node], face_node_ordinals[face_node]);
    }
  }
}

template <unsigned NUM_NODES>
STK_INLINE_FUNCTION
void check_permutation_node_ordinals_ngp(stk::topology topology, unsigned gold_permutation_node_ordinals[][NUM_NODES])
{
  const unsigned numPermutations = topology.num_permutations();
  for (unsigned permutationOrdinal = 0; permutationOrdinal < numPermutations; ++permutationOrdinal) {
    unsigned permutation_node_ordinals[NUM_NODES];
    topology.permutation_node_ordinals(permutationOrdinal, permutation_node_ordinals);
    for (unsigned node = 0; node < NUM_NODES; ++node) {
      NGP_EXPECT_EQ(gold_permutation_node_ordinals[permutationOrdinal][node], permutation_node_ordinals[node]);
    }
  }
}

template <unsigned NUM_NODES>
STK_INLINE_FUNCTION
void fill_permutation_ngp(unsigned permutation, const unsigned permutationNodeOrdinals[][NUM_NODES], unsigned permutedOrdinals[NUM_NODES])
{
  for (unsigned node = 0; node < NUM_NODES; ++node) {
    permutedOrdinals[node] = permutationNodeOrdinals[permutation][node];
  }
}

template <unsigned NUM_NODES>
STK_INLINE_FUNCTION
void check_equivalent_ngp(stk::topology topology, unsigned gold_permutation_node_ordinals[][NUM_NODES])
{
  const unsigned numPermutations = topology.num_permutations();
  unsigned baseNodeArray[NUM_NODES];
  fill_permutation_ngp(0, gold_permutation_node_ordinals, baseNodeArray);

  for (unsigned permutationOrdinal = 0; permutationOrdinal < numPermutations; ++permutationOrdinal) {
    unsigned permutedNodeArray[NUM_NODES];
    fill_permutation_ngp(permutationOrdinal, gold_permutation_node_ordinals, permutedNodeArray);
    stk::EquivalentPermutation result = topology.is_equivalent(baseNodeArray, permutedNodeArray);
    NGP_EXPECT_TRUE(result.is_equivalent);
    NGP_EXPECT_EQ(permutationOrdinal, result.permutation_number);
  }

  // Throw in a few duds that don't match, just for good measure
  for (unsigned node = 0; node < NUM_NODES; ++node) {
    unsigned permutedNodeArray[NUM_NODES];
    for (unsigned initNode = 0; initNode < NUM_NODES; ++initNode) {
      permutedNodeArray[initNode] = baseNodeArray[initNode];
    }
    permutedNodeArray[node] = NUM_NODES;  // Tweak one node to not match
    stk::EquivalentPermutation result = topology.is_equivalent(baseNodeArray, permutedNodeArray);
    NGP_EXPECT_FALSE(result.is_equivalent);
    NGP_EXPECT_EQ(0u, result.permutation_number);
  }
}

inline void check_side_nodes(stk::topology topology, std::vector<std::vector<unsigned>> & gold_side_node_ordinals)
{
  std::vector<unsigned> allElemNodes(topology.num_nodes());
  for (unsigned nodeOrdinal = 0; nodeOrdinal < topology.num_nodes(); ++nodeOrdinal) {
    allElemNodes[nodeOrdinal] = nodeOrdinal + 100;
  }

  for (unsigned side = 0; side < topology.num_sides(); ++side) {
    stk::topology sideTopo = topology.side_topology(side);
    unsigned numSideNodes = (sideTopo.num_nodes() > 0) ? sideTopo.num_nodes() : 1;
    std::vector<unsigned> sideNodes(numSideNodes);

    topology.side_nodes(allElemNodes.data(), side, sideNodes.data());
    for (unsigned sideNodeOrdinal = 0; sideNodeOrdinal < numSideNodes; ++sideNodeOrdinal) {
      EXPECT_EQ(gold_side_node_ordinals[side][sideNodeOrdinal] + 100, sideNodes[sideNodeOrdinal]);
    }
  }
}

inline void check_edge_nodes(stk::topology topology, std::vector<std::vector<unsigned>> & gold_edge_node_ordinals)
{
  std::vector<unsigned> allElemNodes(topology.num_nodes());
  for (unsigned nodeOrdinal = 0; nodeOrdinal < topology.num_nodes(); ++nodeOrdinal) {
    allElemNodes[nodeOrdinal] = nodeOrdinal + 100;
  }

  std::vector<unsigned> edgeNodes;
  for (unsigned edge = 0; edge < topology.num_edges(); ++edge) {
    stk::topology edgeTopo = topology.edge_topology(edge);
    unsigned numEdgeNodes = edgeTopo.num_nodes();
    edgeNodes.resize(numEdgeNodes);
    topology.edge_nodes(allElemNodes.data(), edge, edgeNodes.data());
    for (unsigned edgeNodeOrdinal = 0; edgeNodeOrdinal < numEdgeNodes; ++edgeNodeOrdinal) {
      EXPECT_EQ(gold_edge_node_ordinals[edge][edgeNodeOrdinal] + 100, edgeNodes[edgeNodeOrdinal]);
    }
  }
}

inline void check_face_nodes(stk::topology topology, std::vector<std::vector<unsigned>> & gold_face_node_ordinals)
{
  std::vector<unsigned> allElemNodes(topology.num_nodes());
  for (unsigned nodeOrdinal = 0; nodeOrdinal < topology.num_nodes(); ++nodeOrdinal) {
    allElemNodes[nodeOrdinal] = nodeOrdinal + 100;
  }

  for (unsigned face = 0; face < topology.num_faces(); ++face) {
    stk::topology faceTopo = topology.face_topology(face);
    unsigned numFaceNodes = faceTopo.num_nodes();
    std::vector<unsigned> faceNodes(numFaceNodes);

    topology.face_nodes(allElemNodes.data(), face, faceNodes.data());
    for (unsigned faceNodeOrdinal = 0; faceNodeOrdinal < numFaceNodes; ++faceNodeOrdinal) {
      EXPECT_EQ(gold_face_node_ordinals[face][faceNodeOrdinal] + 100, faceNodes[faceNodeOrdinal]);
    }
  }
}

inline void check_permutation_nodes(stk::topology topology, std::vector<std::vector<unsigned>> & gold_permutation_node_ordinals)
{
  std::vector<unsigned> allNodes(topology.num_nodes());
  for (unsigned nodeOrdinal = 0; nodeOrdinal < topology.num_nodes(); ++nodeOrdinal) {
    allNodes[nodeOrdinal] = nodeOrdinal + 100;
  }

  const unsigned numPermutations = topology.num_permutations();
  for (unsigned permutationOrdinal = 0; permutationOrdinal < numPermutations; ++permutationOrdinal) {
    std::vector<unsigned> permutationNodes(topology.num_nodes());
    topology.permutation_nodes(allNodes.data(), permutationOrdinal, permutationNodes.data());
    for (unsigned nodeOrdinal = 0; nodeOrdinal < topology.num_nodes(); ++nodeOrdinal) {
      EXPECT_EQ(gold_permutation_node_ordinals[permutationOrdinal][nodeOrdinal] + 100, permutationNodes[nodeOrdinal]);
    }
  }
}

constexpr unsigned MAX_NODES_PER_ELEM = 100;

template <unsigned MAX_NODES>
STK_INLINE_FUNCTION
void check_side_nodes_ngp(stk::topology topology, unsigned gold_side_node_ordinals[][MAX_NODES])
{
  unsigned allElemNodes[MAX_NODES_PER_ELEM];
  NGP_EXPECT_TRUE(topology.num_nodes() < MAX_NODES_PER_ELEM);

  for (unsigned nodeOrdinal = 0; nodeOrdinal < topology.num_nodes(); ++nodeOrdinal) {
    allElemNodes[nodeOrdinal] = nodeOrdinal + 100;
  }

  for (unsigned side = 0; side < topology.num_sides(); ++side) {
    stk::topology sideTopo = topology.side_topology(side);
    unsigned numSideNodes = (sideTopo.num_nodes() > 0) ? sideTopo.num_nodes() : 1;
    unsigned sideNodes[MAX_NODES];

    topology.side_nodes(allElemNodes, side, sideNodes);
    for (unsigned sideNodeOrdinal = 0; sideNodeOrdinal < numSideNodes; ++sideNodeOrdinal) {
      NGP_EXPECT_EQ(gold_side_node_ordinals[side][sideNodeOrdinal] + 100, sideNodes[sideNodeOrdinal]);
    }
  }
}

template <unsigned MAX_NODES>
STK_INLINE_FUNCTION
void check_edge_nodes_ngp(stk::topology topology, unsigned gold_edge_node_ordinals[][MAX_NODES])
{
  unsigned allElemNodes[MAX_NODES_PER_ELEM];
  NGP_EXPECT_TRUE(topology.num_nodes() < MAX_NODES_PER_ELEM);

  for (unsigned nodeOrdinal = 0; nodeOrdinal < topology.num_nodes(); ++nodeOrdinal) {
    allElemNodes[nodeOrdinal] = nodeOrdinal + 100;
  }

  unsigned edgeNodes[MAX_NODES];
  for (unsigned edge = 0; edge < topology.num_edges(); ++edge) {
    stk::topology edgeTopo = topology.edge_topology(edge);
    unsigned numEdgeNodes = edgeTopo.num_nodes();
    topology.edge_nodes(allElemNodes, edge, edgeNodes);
    for (unsigned edgeNodeOrdinal = 0; edgeNodeOrdinal < numEdgeNodes; ++edgeNodeOrdinal) {
      NGP_EXPECT_EQ(gold_edge_node_ordinals[edge][edgeNodeOrdinal] + 100, edgeNodes[edgeNodeOrdinal]);
    }
  }
}

template <unsigned MAX_NODES>
STK_INLINE_FUNCTION
void check_face_nodes_ngp(stk::topology topology, unsigned gold_face_node_ordinals[][MAX_NODES])
{
  unsigned allElemNodes[MAX_NODES_PER_ELEM];
  NGP_EXPECT_TRUE(topology.num_nodes() < MAX_NODES_PER_ELEM);

  for (unsigned nodeOrdinal = 0; nodeOrdinal < topology.num_nodes(); ++nodeOrdinal) {
    allElemNodes[nodeOrdinal] = nodeOrdinal + 100;
  }

  for (unsigned face = 0; face < topology.num_faces(); ++face) {
    stk::topology faceTopo = topology.face_topology(face);
    unsigned numFaceNodes = faceTopo.num_nodes();
    unsigned faceNodes[MAX_NODES];

    topology.face_nodes(allElemNodes, face, faceNodes);
    for (unsigned faceNodeOrdinal = 0; faceNodeOrdinal < numFaceNodes; ++faceNodeOrdinal) {
      NGP_EXPECT_EQ(gold_face_node_ordinals[face][faceNodeOrdinal] + 100, faceNodes[faceNodeOrdinal]);
    }
  }
}

template <unsigned NUM_NODES>
STK_INLINE_FUNCTION
void check_permutation_nodes_ngp(stk::topology topology, unsigned gold_permutation_node_ordinals[][NUM_NODES])
{
  unsigned allNodes[NUM_NODES];
  for (unsigned nodeOrdinal = 0; nodeOrdinal < NUM_NODES; ++nodeOrdinal) {
    allNodes[nodeOrdinal] = nodeOrdinal + 100;
  }

  const unsigned numPermutations = topology.num_permutations();
  for (unsigned permutationOrdinal = 0; permutationOrdinal < numPermutations; ++permutationOrdinal) {
    unsigned permutationNodes[NUM_NODES];
    topology.permutation_nodes(allNodes, permutationOrdinal, permutationNodes);
    for (unsigned nodeOrdinal = 0; nodeOrdinal < NUM_NODES; ++nodeOrdinal) {
      NGP_EXPECT_EQ(gold_permutation_node_ordinals[permutationOrdinal][nodeOrdinal] + 100, permutationNodes[nodeOrdinal]);
    }
  }
}

template <unsigned NUM_NODES>
STK_INLINE_FUNCTION
void check_lexicographical_smallest_permutation_ngp(stk::topology topology, unsigned gold_permutation_node_ordinals[][NUM_NODES])
{
  if (topology.num_nodes() == 0) {
    return;
  }

  const unsigned numPermutations = topology.num_permutations();
  unsigned nodeArray[NUM_NODES];

  const unsigned firstPermutation = 0;
  const unsigned lastPermutation = (numPermutations > 0) ? (numPermutations - 1) : 0;

  for (unsigned permWithSmallestSorting = firstPermutation; permWithSmallestSorting <= lastPermutation; ++permWithSmallestSorting) {
    // Fill the nodes with IDs so that the target permutation will sort first
    for (unsigned node = 0; node < NUM_NODES; ++node) {
      unsigned nodeIndex = gold_permutation_node_ordinals[permWithSmallestSorting][node];
      nodeArray[nodeIndex] = node + 1;
    }

    // Check if we retrieve this target permutation as the lexicographically-smallest
    unsigned result = topology.lexicographical_smallest_permutation(nodeArray);
    NGP_EXPECT_EQ(permWithSmallestSorting, result);

    // Check if we retrieve a positive permutation as the lexicographically-smallest when we set
    // the flag.  Can't predict which one we will get, so just check if it's positive.
    const bool onlyPositivePermutations = true;
    result = topology.lexicographical_smallest_permutation(nodeArray, onlyPositivePermutations);
    if (permWithSmallestSorting < topology.num_positive_permutations() || (topology.num_positive_permutations() == 0)) {
      NGP_EXPECT_EQ(permWithSmallestSorting, result);  // Match exactly
    }
    else {
      NGP_EXPECT_LT(result, topology.num_positive_permutations());  // Confirm positive
    }
  }
}
