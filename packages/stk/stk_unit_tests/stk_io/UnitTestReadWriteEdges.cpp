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
#include <stk_mesh/base/Field.hpp>
#include <stk_mesh/base/FEMHelpers.hpp>
#include <stk_mesh/base/SkinBoundary.hpp>
#include <stk_unit_test_utils/GetMeshSpec.hpp>
#include <stk_unit_test_utils/TextMesh.hpp>
#include <stk_unit_test_utils/MeshFixture.hpp>
#include <stk_io/FillMesh.hpp>
#include <stk_io/IossBridge.hpp>
#include <stk_io/StkMeshIoBroker.hpp>
#include <stk_util/parallel/Parallel.hpp>
#include "UnitTestReadWriteEdges.hpp"

void StkEdgeIoTest::setup_edge_mesh(unsigned numBlocks)
{
  std::string meshDesc = stk::unit_test_util::get_many_block_mesh_desc(numBlocks, stk::parallel_machine_size(MPI_COMM_WORLD));
  std::vector<double> coords = stk::unit_test_util::get_many_block_coordinates(numBlocks);
  stk::mesh::Part* edgePart = &get_meta().declare_part_with_topology(edgePartName, stk::topology::LINE_2);
  stk::io::put_io_part_attribute(*edgePart);
  stk::unit_test_util::setup_text_mesh(get_bulk(), meshDesc, coords);

  stk::mesh::create_edges(get_bulk(), get_meta().universal_part(), edgePart);
}

void StkEdgeIoTest::setup_mesh_with_edges(unsigned numBlocks)
{
  setup_empty_mesh(stk::mesh::BulkData::NO_AUTO_AURA);
  setup_edge_mesh(numBlocks);
}

void StkEdgeIoTest::setup_mesh_with_edges_and_faces(unsigned numBlocks)
{
  setup_empty_mesh(stk::mesh::BulkData::NO_AUTO_AURA);
  std::string meshDesc = stk::unit_test_util::get_many_block_mesh_desc(numBlocks, stk::parallel_machine_size(MPI_COMM_WORLD));
  std::vector<double> coords = stk::unit_test_util::get_many_block_coordinates(numBlocks);
  stk::mesh::Part* edgePart = &get_meta().declare_part_with_topology(edgePartName, stk::topology::LINE_2);
  stk::mesh::Part* facePart = &get_meta().declare_part_with_topology(facePartName, stk::topology::QUAD_4);
  stk::io::put_io_part_attribute(*facePart);
  stk::io::put_io_part_attribute(*edgePart);
  stk::unit_test_util::setup_text_mesh(get_bulk(), meshDesc, coords);

  stk::mesh::PartVector faceParts = {facePart};
  stk::mesh::create_interior_block_boundary_sides(get_bulk(), get_meta().universal_part(), faceParts);
  stk::mesh::create_edges(get_bulk(), get_meta().universal_part(), edgePart);
}

void StkEdgeIoTest::test_connectivity_to_element(const stk::mesh::BulkData& bulk, stk::mesh::EntityRank entityRank)
{
  stk::mesh::EntityVector entities;
  const stk::mesh::MetaData& meta = bulk.mesh_meta_data();
  const stk::mesh::BucketVector& buckets = bulk.buckets(entityRank);

  if(bulk.parallel_size() > 1) {

    stk::mesh::Selector ownedAndNotShared = meta.locally_owned_part() & !meta.globally_shared_part();
    stk::mesh::get_selected_entities(ownedAndNotShared, buckets, entities);

    for(stk::mesh::Entity entity: entities) {
      EXPECT_EQ(1u, bulk.num_elements(entity));
    }

    stk::mesh::Selector shared = meta.globally_shared_part();
    stk::mesh::get_selected_entities(shared, buckets, entities);

    unsigned expectedNumElems = bulk.is_automatic_aura_on() ? 2u : 1u;
    for(stk::mesh::Entity entity : entities) {
      EXPECT_EQ(expectedNumElems, bulk.num_elements(entity));
    }
  }
}

void StkEdgeIoTest::test_entity_count(const stk::mesh::BulkData& bulk, stk::mesh::EntityRank entityRank,
                                      unsigned expectedNumLocalEntities, unsigned expectedNumEntities)
{
  stk::mesh::EntityVector entities;
  const stk::mesh::MetaData& meta = bulk.mesh_meta_data();
  const stk::mesh::BucketVector& buckets = bulk.buckets(entityRank);

  unsigned numLocalEntities = stk::mesh::count_selected_entities(meta.locally_owned_part(), buckets);
  EXPECT_EQ(expectedNumLocalEntities, numLocalEntities);

  stk::mesh::Selector selector = meta.locally_owned_part() | meta.globally_shared_part();
  unsigned numEntities = stk::mesh::count_selected_entities(selector, buckets);
  EXPECT_EQ(expectedNumEntities, numEntities);
}

void StkEdgeIoTest::test_edges(const stk::mesh::BulkData& bulk)
{
  test_entity_count(bulk, stk::topology::EDGE_RANK, expectedValues.numLocalEdgesPerProc[bulk.parallel_rank()],
                    expectedValues.numEdgesPerProc[bulk.parallel_rank()]);
  test_connectivity_to_element(bulk, stk::topology::EDGE_RANK);
}

void StkEdgeIoTest::test_faces(const stk::mesh::BulkData& bulk)
{
  test_entity_count(bulk, stk::topology::FACE_RANK, expectedValues.numLocalFacesPerProc[bulk.parallel_rank()],
                    expectedValues.numFacesPerProc[bulk.parallel_rank()]);
  test_connectivity_to_element(bulk, stk::topology::FACE_RANK);
}

void StkEdgeIoTest::output_mesh()
{
  stk::io::StkMeshIoBroker stkIo;
  stkIo.set_bulk_data(get_bulk());
  stkIo.enable_edge_io();
  size_t outputFileIndex = stkIo.create_output_mesh(fileName, stk::io::WRITE_RESULTS);

  stkIo.write_output_mesh(outputFileIndex);
}

void StkEdgeIoTest::test_output_mesh()
{
  stk::mesh::MetaData meta;
  stk::mesh::BulkData bulk(meta, MPI_COMM_WORLD);
  load_output_mesh(bulk);
  test_output_mesh(bulk);
}

void StkEdgeIoTest::load_output_mesh(stk::mesh::BulkData& bulk)
{
  stk::io::fill_mesh(fileName, bulk);
}

void StkEdgeIoTest::test_output_mesh(stk::mesh::BulkData& bulk)
{
  std::vector<size_t> entityCounts;
  stk::mesh::comm_mesh_counts(bulk, entityCounts);
  EXPECT_EQ(expectedValues.globalEdgeCount, entityCounts[stk::topology::EDGE_RANK]);
  EXPECT_EQ(expectedValues.globalElemCount, entityCounts[stk::topology::ELEM_RANK]);

  test_edges(bulk);
  test_faces(bulk);
}

void StkEdgeIoTest::set_expected_values(ExpectedValues& expectedValues_)
{
  expectedValues = expectedValues_;
}

TEST_F(StkEdgeIoTest, SerialWriteMesh)
{
  if(stk::parallel_machine_size(MPI_COMM_WORLD) != 1) { return; }
  ExpectedValues expectedValues;
  expectedValues.numEdgesPerProc = std::vector<unsigned>{12};
  expectedValues.numLocalEdgesPerProc = std::vector<unsigned>{12};
  expectedValues.numFacesPerProc = std::vector<unsigned>{0};
  expectedValues.numLocalFacesPerProc = std::vector<unsigned>{0};
  expectedValues.numConnectedEdges = 0;
  expectedValues.globalEdgeCount = 12;
  expectedValues.globalElemCount = 1;

  setup_mesh_with_edges(1);
  set_expected_values(expectedValues);
  test_edges(get_bulk());
  output_mesh();

  test_output_mesh();
}

TEST_F(StkEdgeIoTest, ParallelWriteMesh)
{
  if(stk::parallel_machine_size(MPI_COMM_WORLD) != 2) { return; }
  ExpectedValues expectedValues;
  expectedValues.numEdgesPerProc = std::vector<unsigned>{12, 12};
  expectedValues.numLocalEdgesPerProc = std::vector<unsigned>{12, 8};
  expectedValues.numFacesPerProc = std::vector<unsigned>{0, 0};
  expectedValues.numLocalFacesPerProc = std::vector<unsigned>{0, 0};
  expectedValues.numConnectedEdges = 0;
  expectedValues.globalEdgeCount = 20;
  expectedValues.globalElemCount = 2;

  setup_mesh_with_edges(2);
  set_expected_values(expectedValues);
  test_edges(get_bulk());
  output_mesh();

  test_output_mesh();
}

TEST_F(StkEdgeIoTest, SerialWriteMeshWithFace)
{
  if(stk::parallel_machine_size(MPI_COMM_WORLD) != 1) { return; }
  ExpectedValues expectedValues;
  expectedValues.numEdgesPerProc = std::vector<unsigned>{20};
  expectedValues.numLocalEdgesPerProc = std::vector<unsigned>{20};
  expectedValues.numFacesPerProc = std::vector<unsigned>{1};
  expectedValues.numLocalFacesPerProc = std::vector<unsigned>{1};
  expectedValues.numConnectedEdges = 4;
  expectedValues.globalEdgeCount = 20;
  expectedValues.globalElemCount = 2;

  setup_mesh_with_edges_and_faces(2);
  set_expected_values(expectedValues);
  test_edges(get_bulk());
  test_faces(get_bulk());
  output_mesh();

  test_output_mesh();
}

TEST_F(StkEdgeIoTest, ParallelWriteMeshWithFace)
{
  if(stk::parallel_machine_size(MPI_COMM_WORLD) != 2) { return; }
  ExpectedValues expectedValues;
  expectedValues.numEdgesPerProc = std::vector<unsigned>{12, 12};
  expectedValues.numLocalEdgesPerProc = std::vector<unsigned>{12, 8};
  expectedValues.numFacesPerProc = std::vector<unsigned>{1, 1};
  expectedValues.numLocalFacesPerProc = std::vector<unsigned>{1, 0};
  expectedValues.numConnectedEdges = 4;
  expectedValues.globalEdgeCount = 20;
  expectedValues.globalElemCount = 2;

  setup_mesh_with_edges_and_faces(2);
  set_expected_values(expectedValues);
  test_edges(get_bulk());
  test_faces(get_bulk());
  output_mesh();

  test_output_mesh();
}