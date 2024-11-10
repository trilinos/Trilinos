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

#include <gtest/gtest.h>
#include <stk_topology/topology.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/GetEntities.hpp>
#include <stk_mesh/baseImpl/MeshImplUtils.hpp>
#include <stk_mesh/baseImpl/ConnectEdgesImpl.hpp>
#include <stk_mesh/base/Types.hpp>
#include <stk_mesh/base/Comm.hpp>
#include <stk_mesh/base/FEMHelpers.hpp>
#include <stk_mesh/base/SkinBoundary.hpp>
#include <stk_unit_test_utils/GetMeshSpec.hpp>
#include <stk_unit_test_utils/TextMesh.hpp>
#include <stk_unit_test_utils/BuildMesh.hpp>
#include <stk_io/FillMesh.hpp>
#include <stk_io/IossBridge.hpp>
#include <stk_io/StkMeshIoBroker.hpp>
#include <stk_util/parallel/Parallel.hpp>

using stk::unit_test_util::build_mesh;

bool is_valid_element_edge(const stk::mesh::BulkData& bulk, stk::mesh::Entity elem, stk::mesh::Entity edge)
{
  EXPECT_EQ(stk::topology::ELEM_RANK, bulk.entity_rank(elem));
  EXPECT_EQ(stk::topology::EDGE_RANK, bulk.entity_rank(edge));

  stk::mesh::EntityVector edgeNodes(bulk.begin_nodes(edge), bulk.end_nodes(edge));

  stk::mesh::OrdinalAndPermutation ordinalAndPerm = stk::mesh::get_ordinal_and_permutation(bulk, elem,
                                                                                           stk::topology::EDGE_RANK, edgeNodes);

  return (ordinalAndPerm.first != stk::mesh::INVALID_CONNECTIVITY_ORDINAL);
}

TEST(EdgeConnectorTest, CheckValidEdge)
{
  if(stk::parallel_machine_size(MPI_COMM_WORLD) > 1) { return; }

  std::shared_ptr<stk::mesh::BulkData> bulkPtr = build_mesh(3, MPI_COMM_WORLD);
  stk::mesh::BulkData& bulk = *bulkPtr;
  stk::io::fill_mesh("generated:1x1x1", bulk);;

  stk::mesh::Entity elem = bulk.get_entity(stk::topology::ELEM_RANK, 1u);
  bulk.modification_begin();
  stk::mesh::Entity edge = stk::mesh::declare_element_edge(bulk, 1u, elem, 0u);
  bulk.modification_end();
  const stk::mesh::Entity* nodes = bulk.begin_nodes(edge);
  unsigned numNodes = bulk.num_nodes(edge);

  EXPECT_EQ(2u, numNodes);
  EXPECT_EQ(1u, bulk.identifier(nodes[0]));
  EXPECT_EQ(2u, bulk.identifier(nodes[1]));

  EXPECT_TRUE(is_valid_element_edge(bulk, elem, edge));
}

TEST(EdgeConnectorTest, CheckInvalidEdge)
{
  if(stk::parallel_machine_size(MPI_COMM_WORLD) > 1) { return; }

  std::shared_ptr<stk::mesh::BulkData> bulkPtr = build_mesh(3, MPI_COMM_WORLD);
  stk::mesh::BulkData& bulk = *bulkPtr;
  stk::io::fill_mesh("generated:1x1x1", bulk);;

  stk::mesh::Entity elem = bulk.get_entity(stk::topology::ELEM_RANK, 1u);
  stk::mesh::Entity node1 = bulk.get_entity(stk::topology::NODE_RANK, 1u);
  stk::mesh::Entity node7 = bulk.get_entity(stk::topology::NODE_RANK, 7u);

  stk::mesh::Part& edgePart = bulk.mesh_meta_data().get_topology_root_part(stk::topology::LINE_2);
  stk::mesh::PartVector edgeParts = { &edgePart };

  bulk.modification_begin();
  stk::mesh::Entity edge = bulk.declare_edge(1u, edgeParts);
  bulk.declare_relation(edge, node1, 0u);
  bulk.declare_relation(edge, node7, 1u);
  bulk.modification_end();

  const stk::mesh::Entity* nodes = bulk.begin_nodes(edge);
  unsigned numNodes = bulk.num_nodes(edge);

  EXPECT_EQ(2u, numNodes);
  EXPECT_EQ(1u, bulk.identifier(nodes[0]));
  EXPECT_EQ(7u, bulk.identifier(nodes[1]));

  EXPECT_FALSE(is_valid_element_edge(bulk, elem, edge));
}

TEST(EdgeConnectorTest, CheckConnectElemToValidEdge)
{
  if(stk::parallel_machine_size(MPI_COMM_WORLD) > 1) { return; }

  std::shared_ptr<stk::mesh::BulkData> bulkPtr = build_mesh(3, MPI_COMM_WORLD);
  stk::mesh::BulkData& bulk = *bulkPtr;
  stk::io::fill_mesh("generated:1x1x1", bulk);;

  stk::mesh::Entity elem = bulk.get_entity(stk::topology::ELEM_RANK, 1u);
  stk::mesh::Entity node1 = bulk.get_entity(stk::topology::NODE_RANK, 1u);
  stk::mesh::Entity node2 = bulk.get_entity(stk::topology::NODE_RANK, 2u);
  stk::mesh::Part& edgePart = bulk.mesh_meta_data().get_topology_root_part(stk::topology::LINE_2);
  stk::mesh::PartVector edgeParts = { &edgePart };

  bulk.modification_begin();
  stk::mesh::Entity edge = bulk.declare_edge(1u, edgeParts);
  bulk.declare_relation(edge, node1, 0u);
  bulk.declare_relation(edge, node2, 1u);
  stk::mesh::impl::connect_edge_to_elements(bulk, edge);
  bulk.modification_end();

  EXPECT_EQ(1u, bulk.num_elements(edge));
  const stk::mesh::Entity* edgeElement = bulk.begin_elements(edge);
  EXPECT_EQ(edgeElement[0], elem);

  EXPECT_EQ(1u, bulk.num_edges(elem));
  const stk::mesh::Entity* elementEdge = bulk.begin_edges(elem);
  EXPECT_EQ(elementEdge[0], edge);
}

TEST(EdgeConnectorTest, CheckConnectElemToInvalidEdge)
{
  if(stk::parallel_machine_size(MPI_COMM_WORLD) > 1) { return; }

  std::shared_ptr<stk::mesh::BulkData> bulkPtr = build_mesh(3, MPI_COMM_WORLD);
  stk::mesh::BulkData& bulk = *bulkPtr;
  stk::io::fill_mesh("generated:1x1x1", bulk);;

  stk::mesh::Entity elem = bulk.get_entity(stk::topology::ELEM_RANK, 1u);
  stk::mesh::Entity node1 = bulk.get_entity(stk::topology::NODE_RANK, 1u);
  stk::mesh::Entity node7 = bulk.get_entity(stk::topology::NODE_RANK, 7u);
  stk::mesh::Part& edgePart = bulk.mesh_meta_data().get_topology_root_part(stk::topology::LINE_2);
  stk::mesh::PartVector edgeParts = { &edgePart };

  bulk.modification_begin();
  stk::mesh::Entity edge = bulk.declare_edge(1u, edgeParts);
  bulk.declare_relation(edge, node1, 0u);
  bulk.declare_relation(edge, node7, 1u);
  EXPECT_FALSE(stk::mesh::impl::connect_edge_to_elements(bulk, edge));
  bulk.modification_end();

  EXPECT_EQ(0u, bulk.num_elements(edge));
  EXPECT_EQ(0u, bulk.num_edges(elem));
}

TEST(StkEdgeIo, ParallelWriteMesh)
{
  if(stk::parallel_machine_size(MPI_COMM_WORLD) != 2) { return; }

  std::string filename = "output.exo";
  {
    std::string meshDesc = stk::unit_test_util::get_many_block_mesh_desc(2, 2);
    std::vector<double> coords = stk::unit_test_util::get_many_block_coordinates(2);
    std::shared_ptr<stk::mesh::BulkData> bulkPtr = build_mesh(3, MPI_COMM_WORLD);
    stk::mesh::BulkData& bulk = *bulkPtr;
    stk::mesh::MetaData& meta = bulkPtr->mesh_meta_data();
    stk::mesh::Part* part = &meta.declare_part_with_topology("edgeBlock", stk::topology::LINE_2);
    stk::io::put_edge_block_io_part_attribute(*part);
    stk::unit_test_util::setup_text_mesh(bulk, stk::unit_test_util::get_full_text_mesh_desc(meshDesc, coords));
    stk::mesh::create_edges(bulk, meta.universal_part(), part);
    const stk::mesh::BucketVector& buckets = bulk.buckets(stk::topology::EDGE_RANK);

    unsigned numLocalEdges = stk::mesh::count_selected_entities(meta.locally_owned_part(), buckets);
    unsigned expectedVal = (bulk.parallel_rank() == 0) ? 12u : 8u;
    EXPECT_EQ(expectedVal, numLocalEdges);

    stk::io::StkMeshIoBroker stkIo;
    stkIo.set_bulk_data(bulk);
    size_t outputFileIndex = stkIo.create_output_mesh(filename, stk::io::WRITE_RESULTS);
    stkIo.write_output_mesh(outputFileIndex);
  }

  {
    std::shared_ptr<stk::mesh::BulkData> bulkPtr = build_mesh(MPI_COMM_WORLD);
    stk::mesh::MetaData& meta = bulkPtr->mesh_meta_data();
    stk::mesh::BulkData& bulk = *bulkPtr;
    stk::io::fill_mesh(filename, bulk);
    const stk::mesh::BucketVector& buckets = bulk.buckets(stk::topology::EDGE_RANK);

    std::vector<size_t> entityCounts;
    stk::mesh::comm_mesh_counts(bulk, entityCounts);
    EXPECT_EQ(20u, entityCounts[stk::topology::EDGE_RANK]);
    EXPECT_EQ(2u, entityCounts[stk::topology::ELEM_RANK]);

    unsigned numLocalEdges = stk::mesh::count_selected_entities(meta.locally_owned_part(), buckets);
    unsigned expectedVal = (bulk.parallel_rank() == 0) ? 12u : 8u;
    EXPECT_EQ(expectedVal, numLocalEdges);
  }

  std::string pllFileName = filename+"."+std::to_string(stk::parallel_machine_size(MPI_COMM_WORLD))+"."+std::to_string(stk::parallel_machine_rank(MPI_COMM_WORLD));
  unlink(pllFileName.c_str());
}

TEST(StkEdgeIo, ParallelWriteMeshWithFace)
{
  if(stk::parallel_machine_size(MPI_COMM_WORLD) != 2) { return; }

  std::string filename = "output.exo";
  {
    std::string meshDesc = stk::unit_test_util::get_many_block_mesh_desc(2, 2);
    std::vector<double> coords = stk::unit_test_util::get_many_block_coordinates(2);
    std::shared_ptr<stk::mesh::BulkData> bulkPtr = build_mesh(3, MPI_COMM_WORLD);
    stk::mesh::MetaData& meta = bulkPtr->mesh_meta_data();
    stk::mesh::BulkData& bulk = *bulkPtr;
    stk::mesh::Part* edgePart = &meta.declare_part_with_topology("edgeBlock", stk::topology::LINE_2);
    stk::mesh::Part* facePart = &meta.declare_part_with_topology("faceBlock", stk::topology::QUAD_4);
    stk::io::put_edge_block_io_part_attribute(*edgePart);
    stk::io::put_io_part_attribute(*facePart);
    stk::unit_test_util::setup_text_mesh(bulk, stk::unit_test_util::get_full_text_mesh_desc(meshDesc, coords));

    stk::mesh::PartVector faceParts = {facePart};
    stk::mesh::create_interior_block_boundary_sides(bulk, meta.universal_part(), faceParts);
    stk::mesh::create_edges(bulk, meta.universal_part(), edgePart);
    const stk::mesh::BucketVector& buckets = bulk.buckets(stk::topology::EDGE_RANK);

    unsigned numLocalEdges = stk::mesh::count_selected_entities(meta.locally_owned_part(), buckets);
    unsigned expectedVal = (bulk.parallel_rank() == 0) ? 12u : 8u;
    EXPECT_EQ(expectedVal, numLocalEdges);

    stk::mesh::EntityVector faces;
    stk::mesh::get_selected_entities(*facePart, bulk.buckets(stk::topology::FACE_RANK), faces);
    EXPECT_EQ(1u, faces.size());
    EXPECT_EQ(4u, bulk.num_edges(faces[0]));

    stk::io::StkMeshIoBroker stkIo;
    stkIo.set_bulk_data(bulk);
    size_t outputFileIndex = stkIo.create_output_mesh(filename, stk::io::WRITE_RESULTS);
    stkIo.write_output_mesh(outputFileIndex);
  }

  {
    std::shared_ptr<stk::mesh::BulkData> bulkPtr = build_mesh(3, MPI_COMM_WORLD);
    stk::mesh::MetaData& meta = bulkPtr->mesh_meta_data();
    stk::mesh::BulkData& bulk = *bulkPtr;
    stk::io::fill_mesh(filename, bulk);
    const stk::mesh::BucketVector& buckets = bulk.buckets(stk::topology::EDGE_RANK);

    std::vector<size_t> entityCounts;
    stk::mesh::comm_mesh_counts(bulk, entityCounts);
    EXPECT_EQ(20u, entityCounts[stk::topology::EDGE_RANK]);
    EXPECT_EQ(2u, entityCounts[stk::topology::ELEM_RANK]);

    unsigned numLocalEdges = stk::mesh::count_selected_entities(meta.locally_owned_part(), buckets);
    unsigned expectedVal = (bulk.parallel_rank() == 0) ? 12u : 8u;
    EXPECT_EQ(expectedVal, numLocalEdges);

    stk::mesh::EntityVector faces;
    stk::mesh::get_selected_entities(meta.universal_part(), bulk.buckets(stk::topology::FACE_RANK), faces);
    EXPECT_EQ(1u, faces.size());
    EXPECT_EQ(4u, bulk.num_edges(faces[0]));
  }

  std::string pllFileName = filename+"."+std::to_string(stk::parallel_machine_size(MPI_COMM_WORLD))+"."+std::to_string(stk::parallel_machine_rank(MPI_COMM_WORLD));
  unlink(pllFileName.c_str());
}
