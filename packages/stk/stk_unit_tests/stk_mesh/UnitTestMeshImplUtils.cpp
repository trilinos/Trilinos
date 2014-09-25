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

#include <stddef.h>                     // for size_t
#include <stk_mesh/base/BulkData.hpp>   // for BulkData
#include <stk_mesh/base/GetEntities.hpp>       // for comm_mesh_counts, count_entities
#include <stk_mesh/base/Comm.hpp>       // for comm_mesh_counts
#include <stk_mesh/base/MetaData.hpp>   // for MetaData
#include <gtest/gtest.h>
#include <vector>                       // for vector, vector<>::iterator
#include "stk_mesh/base/Bucket.hpp"     // for Bucket
#include "stk_mesh/base/Types.hpp"      // for BucketVector, EntityRank
#include "stk_topology/topology.hpp"    // for topology, etc
#include "unit_tests/Setup2Block2HexMesh.hpp"
#include "stk_mesh/baseImpl/MeshImplUtils.hpp"

TEST ( UnitTestMeshImplUtils, find_elements_these_nodes_have_in_common )
{
  stk::ParallelMachine communicator = MPI_COMM_WORLD;
  int numProcs = stk::parallel_machine_size(communicator);
  if (numProcs > 2) {
    return;
  }

  const unsigned spatialDim = 3;
  stk::mesh::MetaData meta(spatialDim);
  stk::mesh::BulkData bulk(meta, communicator);

  setup2Block2HexMesh(bulk);

  const unsigned numNodesPerEdge = 2;

  //edge 2-6 is connected to elements 1 and 2
  stk::mesh::EntityId edge_2_6_nodeIds[] = {2, 6};
  stk::mesh::Entity edge_2_6_nodes[numNodesPerEdge];
  edge_2_6_nodes[0] = bulk.get_entity(stk::topology::NODE_RANK, edge_2_6_nodeIds[0]);
  edge_2_6_nodes[1] = bulk.get_entity(stk::topology::NODE_RANK, edge_2_6_nodeIds[1]);

  //edge 5-6 is only connected to element 1
  stk::mesh::EntityId edge_5_6_nodeIds[] = {5, 6};
  stk::mesh::Entity edge_5_6_nodes[numNodesPerEdge];
  edge_5_6_nodes[0] = bulk.get_entity(stk::topology::NODE_RANK, edge_5_6_nodeIds[0]);
  edge_5_6_nodes[1] = bulk.get_entity(stk::topology::NODE_RANK, edge_5_6_nodeIds[1]);

  std::vector<stk::mesh::Entity> elements;

  stk::mesh::impl::find_elements_these_nodes_have_in_common(bulk, numNodesPerEdge, edge_2_6_nodes, elements);

  size_t expected_num_elements = 2;
  EXPECT_EQ(expected_num_elements, elements.size());

  stk::mesh::impl::find_elements_these_nodes_have_in_common(bulk, numNodesPerEdge, edge_5_6_nodes, elements);

  expected_num_elements = 1;
  EXPECT_EQ(expected_num_elements, elements.size());
}

TEST ( UnitTestMeshImplUtils, find_locally_owned_elements_these_nodes_have_in_common )
{
  stk::ParallelMachine communicator = MPI_COMM_WORLD;
  int numProcs = stk::parallel_machine_size(communicator);
  if (numProcs > 2) {
    return;
  }

  const unsigned spatialDim = 3;
  stk::mesh::MetaData meta(spatialDim);
  stk::mesh::BulkData bulk(meta, communicator);

  setup2Block2HexMesh(bulk);

  const unsigned numNodesPerEdge = 2;

  //edge 2-6 is connected to elements 1 and 2
  stk::mesh::EntityId edge_2_6_nodeIds[] = {2, 6};
  stk::mesh::Entity edge_2_6_nodes[numNodesPerEdge];
  edge_2_6_nodes[0] = bulk.get_entity(stk::topology::NODE_RANK, edge_2_6_nodeIds[0]);
  edge_2_6_nodes[1] = bulk.get_entity(stk::topology::NODE_RANK, edge_2_6_nodeIds[1]);

  //edge 5-6 is only connected to element 1
  stk::mesh::EntityId edge_5_6_nodeIds[] = {5, 6};
  stk::mesh::Entity edge_5_6_nodes[numNodesPerEdge];
  edge_5_6_nodes[0] = bulk.get_entity(stk::topology::NODE_RANK, edge_5_6_nodeIds[0]);
  edge_5_6_nodes[1] = bulk.get_entity(stk::topology::NODE_RANK, edge_5_6_nodeIds[1]);

  std::vector<stk::mesh::Entity> elements;

  stk::mesh::impl::find_locally_owned_elements_these_nodes_have_in_common(bulk, numNodesPerEdge, edge_2_6_nodes, elements);

  size_t expected_num_elements = 1;
  if (numProcs == 1) {
    expected_num_elements = 2;
  }
  EXPECT_EQ(expected_num_elements, elements.size());

  stk::mesh::impl::find_locally_owned_elements_these_nodes_have_in_common(bulk, numNodesPerEdge, edge_5_6_nodes, elements);

  expected_num_elements = 0;
  if (bulk.parallel_rank() == 0) {
    //edge_5_6 is connected to element 1, which is locally-owned on proc 0
    expected_num_elements = 1;
  }
  EXPECT_EQ(expected_num_elements, elements.size());
}

TEST ( UnitTestMeshImplUtils, find_element_edge_ordinal_and_equivalent_nodes )
{
  stk::ParallelMachine communicator = MPI_COMM_WORLD;
  int numProcs = stk::parallel_machine_size(communicator);
  if (numProcs > 2) {
    return;
  }

  const unsigned spatialDim = 3;
  stk::mesh::MetaData meta(spatialDim);
  stk::mesh::BulkData bulk(meta, communicator);

  setup2Block2HexMesh(bulk);

  const unsigned numNodesPerEdge = 2;

  //we know that edge 2-6 is edge-ordinal 9 on element 1 and
  //edge-ordinal 8 on element 2
  stk::mesh::EntityId edge_2_6_nodeIds[] = {2, 6};
  stk::mesh::Entity edge_2_6_nodes[numNodesPerEdge];
  edge_2_6_nodes[0] = bulk.get_entity(stk::topology::NODE_RANK, edge_2_6_nodeIds[0]);
  edge_2_6_nodes[1] = bulk.get_entity(stk::topology::NODE_RANK, edge_2_6_nodeIds[1]);

  stk::mesh::EntityId elem1Id = 1;
  stk::mesh::EntityId elem2Id = 2;
  stk::mesh::Entity elem1 = bulk.get_entity(stk::topology::ELEM_RANK, elem1Id);
  stk::mesh::Entity elem2 = bulk.get_entity(stk::topology::ELEM_RANK, elem2Id);

  stk::mesh::Entity elemEdgeNodes[numNodesPerEdge];
  unsigned elemEdgeOrdinal = 999;

  bool found_it = stk::mesh::impl::find_element_edge_ordinal_and_equivalent_nodes(bulk, elem1, numNodesPerEdge, edge_2_6_nodes,
                      elemEdgeOrdinal, elemEdgeNodes);

  EXPECT_EQ(true, found_it);
  unsigned expectedElemEdgeOrdinal = 9;
  EXPECT_EQ(expectedElemEdgeOrdinal, elemEdgeOrdinal);
  EXPECT_EQ(edge_2_6_nodes[0], elemEdgeNodes[0]);
  EXPECT_EQ(edge_2_6_nodes[1], elemEdgeNodes[1]);

  found_it = stk::mesh::impl::find_element_edge_ordinal_and_equivalent_nodes(bulk, elem2, numNodesPerEdge, edge_2_6_nodes,
                      elemEdgeOrdinal, elemEdgeNodes);

  EXPECT_EQ(true, found_it);
  expectedElemEdgeOrdinal = 8;
  EXPECT_EQ(expectedElemEdgeOrdinal, elemEdgeOrdinal);
  EXPECT_EQ(edge_2_6_nodes[0], elemEdgeNodes[0]);
  EXPECT_EQ(edge_2_6_nodes[1], elemEdgeNodes[1]);
}

