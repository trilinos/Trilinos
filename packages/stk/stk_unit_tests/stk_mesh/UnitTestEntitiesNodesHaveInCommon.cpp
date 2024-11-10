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

#include "Setup2Block2HexMesh.hpp"
#include "stk_io/StkMeshIoBroker.hpp"
#include "stk_mesh/base/Bucket.hpp"     // for Bucket
#include "stk_mesh/base/Types.hpp"      // for BucketVector, EntityRank
#include "stk_mesh/baseImpl/MeshImplUtils.hpp"
#include "stk_topology/topology.hpp"    // for topology, etc
#include <gtest/gtest.h>
#include <stddef.h>                     // for size_t
#include <stk_mesh/base/BulkData.hpp>   // for BulkData
#include <stk_mesh/base/Comm.hpp>       // for comm_mesh_counts
#include <stk_mesh/base/CreateFaces.hpp>
#include <stk_mesh/base/GetEntities.hpp> // for count_entities
#include <stk_mesh/base/MetaData.hpp>   // for MetaData
#include <stk_unit_test_utils/BulkDataTester.hpp>
#include <stk_unit_test_utils/MeshFixture.hpp>
#include <stk_unit_test_utils/getOption.h>
#include <stk_unit_test_utils/BuildMesh.hpp>
#include <vector>                       // for vector, vector<>::iterator

namespace {
using stk::unit_test_util::build_mesh;

TEST ( MeshImplUtils, find_elements_these_nodes_have_in_common )
{
  stk::ParallelMachine communicator = MPI_COMM_WORLD;
  int numProcs = stk::parallel_machine_size(communicator);
  if (numProcs > 2) {
    return;
  }

  const unsigned spatialDim = 3;
  std::shared_ptr<stk::mesh::BulkData> bulkPtr = build_mesh(spatialDim, communicator);
  stk::mesh::BulkData& bulk = *bulkPtr;

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

  stk::mesh::impl::find_entities_these_nodes_have_in_common(bulk, stk::topology::ELEMENT_RANK, numNodesPerEdge, edge_2_6_nodes, elements);

  size_t expected_num_elements = 2;
  EXPECT_EQ(expected_num_elements, elements.size());

  stk::mesh::impl::find_entities_these_nodes_have_in_common(bulk, stk::topology::ELEMENT_RANK, numNodesPerEdge, edge_5_6_nodes, elements);

  expected_num_elements = 1;
  EXPECT_EQ(expected_num_elements, elements.size());
}

TEST ( MeshImplUtils, find_faces_these_nodes_have_in_common )
{
  stk::ParallelMachine communicator = MPI_COMM_WORLD;
  int numProcs = stk::parallel_machine_size(communicator);
  if (numProcs > 2) {
    return;
  }

  const unsigned spatialDim = 3;
  std::shared_ptr<stk::mesh::BulkData> bulkPtr = build_mesh(spatialDim, communicator);
  stk::mesh::BulkData& bulk = *bulkPtr;

  setup2Block2HexMesh(bulk);
  stk::mesh::create_faces(bulk);

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

  std::vector<stk::mesh::Entity> faces;

  stk::mesh::impl::find_entities_these_nodes_have_in_common(bulk, stk::topology::FACE_RANK, numNodesPerEdge, &edge_2_6_nodes[0],faces);

  size_t expected_num_faces = 3;
  EXPECT_EQ(expected_num_faces, faces.size());

  stk::mesh::impl::find_entities_these_nodes_have_in_common(bulk, stk::topology::FACE_RANK, numNodesPerEdge, &edge_5_6_nodes[0],faces);

  expected_num_faces = 2;
  EXPECT_EQ(expected_num_faces, faces.size());
}

TEST ( MeshImplUtils, find_locally_owned_elements_these_nodes_have_in_common )
{
  stk::ParallelMachine communicator = MPI_COMM_WORLD;
  int numProcs = stk::parallel_machine_size(communicator);
  if (numProcs > 2) {
    return;
  }

  const unsigned spatialDim = 3;
  std::shared_ptr<stk::mesh::BulkData> bulkPtr = build_mesh(spatialDim, communicator);
  stk::mesh::BulkData& bulk = *bulkPtr;

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

TEST(MeshImplUtils, do_these_nodes_have_any_shell_elements_in_common_no_nodes)
{
  const int spatialDim = 3;
  std::shared_ptr<stk::mesh::BulkData> bulkPtr = build_mesh(spatialDim, MPI_COMM_WORLD);
  stk::mesh::BulkData& mesh = *bulkPtr;
  EXPECT_FALSE(stk::mesh::impl::do_these_nodes_have_any_shell_elements_in_common(mesh, 0, NULL));
}

TEST(MeshImplUtils, do_these_nodes_have_any_shell_elements_in_common_hexshell)
{
  const int spatialDim = 3;
  std::shared_ptr<stk::mesh::BulkData> bulkPtr = build_mesh(spatialDim, MPI_COMM_WORLD);
  stk::mesh::BulkData& mesh = *bulkPtr;
  if (mesh.parallel_size() == 1) {
    stk::io::StkMeshIoBroker exodus_file_reader(MPI_COMM_WORLD);
    exodus_file_reader.set_bulk_data(mesh);
    exodus_file_reader.add_mesh_database("generated:1x1x1|shell:X", stk::io::READ_MESH);
    exodus_file_reader.create_input_mesh();
    exodus_file_reader.populate_bulk_data();
    {
      stk::mesh::Entity shell = (*mesh.buckets(stk::topology::ELEMENT_RANK)[1])[0];
      ASSERT_TRUE(mesh.bucket(shell).topology().is_shell());
      stk::mesh::EntityVector nodes(4);
      for (unsigned i=0 ; i<4 ; ++i) {
        nodes[i] = mesh.begin_nodes(shell)[i];
      }
      EXPECT_TRUE(stk::mesh::impl::do_these_nodes_have_any_shell_elements_in_common(mesh, 4, nodes.data()));
    }
    {
      stk::mesh::Entity hex = (*mesh.buckets(stk::topology::ELEMENT_RANK)[0])[0];
      ASSERT_FALSE(mesh.bucket(hex).topology().is_shell());
      stk::mesh::EntityVector nodes(8);
      for (unsigned i=0 ; i<8 ; ++i) {
        nodes[i] = mesh.begin_nodes(hex)[i];
      }
      EXPECT_FALSE(stk::mesh::impl::do_these_nodes_have_any_shell_elements_in_common(mesh, 8, nodes.data()));
    }

  }
}

TEST(MeshImplUtils, do_these_nodes_have_any_shell_elements_in_common_hexshellwrap)
{
  const int spatialDim = 3;
  std::shared_ptr<stk::mesh::BulkData> bulkPtr = build_mesh(spatialDim, MPI_COMM_WORLD);
  stk::mesh::BulkData& mesh = *bulkPtr;
  if (mesh.parallel_size() == 1) {
    stk::io::StkMeshIoBroker exodus_file_reader(MPI_COMM_WORLD);
    exodus_file_reader.set_bulk_data(mesh);
    exodus_file_reader.add_mesh_database("generated:1x1x1|shell:xXyY", stk::io::READ_MESH);
    exodus_file_reader.create_input_mesh();
    exodus_file_reader.populate_bulk_data();
    stk::mesh::Entity hex = (*mesh.buckets(stk::topology::ELEMENT_RANK)[0])[0];
    ASSERT_FALSE(mesh.bucket(hex).topology().is_shell());
    stk::mesh::EntityVector nodes(4);
    for (unsigned i=0 ; i<4 ; ++i) {
      nodes[i] = mesh.begin_nodes(hex)[i];
    }
    EXPECT_FALSE(stk::mesh::impl::do_these_nodes_have_any_shell_elements_in_common(mesh, 4, nodes.data()));
  }
}

template<typename T>
std::string vec_to_string(const std::vector<T> &vec)
{
  std::string s = "[";
  for(const T &item : vec)
    s += " " + std::to_string(item);
  return s + " ]";
}

class EntitiesNodesHaveInCommon : public stk::unit_test_util::MeshFixture
{
protected:
  void expect_nodes_have_elems_in_common(stk::mesh::EntityId elemId,
                                         const stk::mesh::EntityIdVector &nodeIds,
                                         const stk::mesh::EntityIdVector &goldElemIds)
  {
    stk::mesh::EntityVector nodes(nodeIds.size());
    for(unsigned i=0; i<nodes.size(); i++)
      nodes[i] = get_bulk().get_entity(stk::topology::NODE_RANK, nodeIds[i]);

    stk::mesh::EntityVector elementsInCommon;
    stk::mesh::impl::find_entities_with_larger_ids_these_nodes_have_in_common_and_locally_owned(elemId,
                                                                                                get_bulk(),
                                                                                                stk::topology::ELEM_RANK,
                                                                                                nodes.size(),
                                                                                                nodes.data(),
                                                                                                elementsInCommon);
    stk::mesh::EntityIdVector elemIds = entities_to_ids(elementsInCommon);
    ASSERT_EQ(goldElemIds.size(), elemIds.size())
        << "expected " << vec_to_string(goldElemIds) << " got " << vec_to_string(elemIds);
    std::sort(elemIds.begin(), elemIds.end());
    for(size_t i=0; i<goldElemIds.size(); i++)
      EXPECT_EQ(goldElemIds[i], elemIds[i]);
  }

  stk::mesh::EntityIdVector entities_to_ids(const stk::mesh::EntityVector &entities)
  {
    stk::mesh::EntityIdVector ids(entities.size());
    for(size_t i=0; i<entities.size(); i++)
      ids[i] = get_bulk().identifier(entities[i]);
    return ids;
  }
};

TEST_F(EntitiesNodesHaveInCommon, 2elems)
{
  if(stk::parallel_machine_size(get_comm()) == 1)
  {
    setup_mesh("generated:1x1x2", stk::mesh::BulkData::NO_AUTO_AURA);

    expect_nodes_have_elems_in_common(0, {1, 2, 3, 4}, {1});
    expect_nodes_have_elems_in_common(1, {1, 2, 3, 4}, {});
    expect_nodes_have_elems_in_common(0, {5, 6, 7, 8}, {1, 2});
    expect_nodes_have_elems_in_common(1, {5, 6, 7, 8}, {2});
    expect_nodes_have_elems_in_common(0, {9, 10, 11, 12}, {2});

    expect_nodes_have_elems_in_common(0, {5, 9}, {2});
  }
}

TEST_F(EntitiesNodesHaveInCommon, 8elems)
{
  if(stk::parallel_machine_size(get_comm()) == 1)
  {
    setup_mesh("generated:2x2x2", stk::mesh::BulkData::NO_AUTO_AURA);

    expect_nodes_have_elems_in_common(0, {14}, {1, 2, 3, 4, 5, 6, 7, 8});
    expect_nodes_have_elems_in_common(0, {14, 13}, {1, 3, 5, 7});
    expect_nodes_have_elems_in_common(0, {14, 13, 10}, {1, 5});
    expect_nodes_have_elems_in_common(0, {14, 13, 10, 1}, {1});
    expect_nodes_have_elems_in_common(0, {14, 13, 10, 11}, {1, 5});

    expect_nodes_have_elems_in_common(1, {14}, {2, 3, 4, 5, 6, 7, 8});
    expect_nodes_have_elems_in_common(1, {14, 13}, {3, 5, 7});
    expect_nodes_have_elems_in_common(1, {14, 13, 10}, {5});
    expect_nodes_have_elems_in_common(1, {14, 13, 10, 11}, {5});
  }
}

}
