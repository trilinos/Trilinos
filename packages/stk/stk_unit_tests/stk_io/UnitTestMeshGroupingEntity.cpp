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

#include <stddef.h>                     // for size_t
#include <stdlib.h>                     // for rand, srand, RAND_MAX
#include <stk_io/IossBridge.hpp>        // for is_part_io_part
#include <stk_io/StkMeshIoBroker.hpp>   // for StkMeshIoBroker
#include <stk_mesh/base/BulkData.hpp>   // for BulkData
#include <stk_mesh/base/MetaData.hpp>   // for MetaData
#include <gtest/gtest.h>
#include <string>                       // for string
#include <vector>                       // for vector, etc
#include "gtest/gtest.h"                // for AssertHelper, ASSERT_TRUE
#include "stk_io/DatabasePurpose.hpp"   // for DatabasePurpose::READ_MESH, etc
#include "stk_mesh/base/Bucket.hpp"     // for Bucket
#include "stk_mesh/base/Entity.hpp"     // for Entity
#include "stk_mesh/base/Part.hpp"       // for Part
#include "stk_mesh/base/Selector.hpp"   // for Selector, operator&
#include "stk_mesh/base/Types.hpp"      // for PartVector, BucketVector, etc
#include "stk_topology/topology.hpp"    // for topology, etc
#include "stk_util/parallel/Parallel.hpp"  // for ParallelMachine
#include "stk_mesh/baseImpl/MeshImplUtils.hpp"
#include "stk_mesh/base/FEMHelpers.hpp"
#include "stk_mesh/base/Field.hpp"
#include <stk_unit_test_utils/MeshFixture.hpp>
#include <Ioss_EntityType.h>            // for EntityType
#include <Ioss_ElementBlock.h>
#include <Ioss_SideSet.h>
#include <Ioss_SideBlock.h>
#include <Ioss_ConcreteVariableType.h>
#include <Ionit_Initializer.h>
#include <Ioss_IOFactory.h>
#include <Ioss_DBUsage.h>
#include <stk_unit_test_utils/TextMesh.hpp>
#include <stk_unit_test_utils/BuildMesh.hpp>

using stk::unit_test_util::build_mesh;
using stk::unit_test_util::build_mesh_no_simple_fields;

TEST(MeshGroupingEntity, universalSideset_get_entity_rank)
{
  Ioss::Init::Initializer init;
  Ioss::PropertyManager properties;
  Ioss::DatabaseIO *db_io = Ioss::IOFactory::create("generated", "1x1x32", Ioss::READ_MODEL, MPI_COMM_WORLD, properties);
  Ioss::Region* ioRegion = new Ioss::Region(db_io, "IO_Region");

  Ioss::SideSet* usideset = new Ioss::SideSet(db_io, "universal_sideset");
  Ioss::SideBlock* usideblk = new Ioss::SideBlock(db_io, "universal_sideset", "unknown", "unknown", 1);
  usideset->add(usideblk);
  stk::mesh::MetaData meta(3);
  EXPECT_EQ(stk::topology::FACE_RANK, stk::io::get_entity_rank(usideset, meta));
  EXPECT_EQ(stk::topology::FACE_RANK, stk::io::get_entity_rank(usideblk, meta));
  delete usideset;
  delete ioRegion;
}

TEST(MeshGroupingEntity, getElementBlockFromPart)
{
  if(stk::parallel_machine_size(MPI_COMM_WORLD) > 1) { return; }
  std::shared_ptr<stk::mesh::BulkData> bulk = build_mesh(MPI_COMM_WORLD);
  stk::mesh::MetaData& meta = bulk->mesh_meta_data();

  stk::io::StkMeshIoBroker ioBroker;
  stk::io::fill_mesh_preexisting(ioBroker, "generated:1x1x1", *bulk);
  std::shared_ptr<Ioss::Region> ioRegion = ioBroker.get_input_ioss_region();

  stk::mesh::Part* part = meta.get_part("block_1");
  EXPECT_TRUE(part != nullptr);
  EXPECT_EQ(stk::topology::ELEMENT_RANK, part->primary_entity_rank());

  Ioss::GroupingEntity* entity = stk::io::get_grouping_entity(*ioRegion, *part);
  EXPECT_TRUE(entity != nullptr);
  EXPECT_EQ(Ioss::ELEMENTBLOCK, entity->type());
}

TEST(MeshGroupingEntity, getElementBlockFromNonIoPart)
{
  if(stk::parallel_machine_size(MPI_COMM_WORLD) > 1) { return; }
  const unsigned spatialDim = 3;
  std::shared_ptr<stk::mesh::BulkData> bulk = build_mesh(spatialDim, MPI_COMM_WORLD);
  stk::mesh::MetaData& meta = bulk->mesh_meta_data();

  stk::io::StkMeshIoBroker ioBroker;
  stk::mesh::Part& dummyPart = meta.declare_part("dummy", stk::topology::ELEMENT_RANK);
  EXPECT_FALSE(stk::io::is_part_io_part(dummyPart));
  stk::io::fill_mesh_preexisting(ioBroker, "generated:1x1x1", *bulk);
  std::shared_ptr<Ioss::Region> ioRegion = ioBroker.get_input_ioss_region();

  Ioss::GroupingEntity* entity = stk::io::get_grouping_entity(*ioRegion, dummyPart);
  EXPECT_TRUE(entity == nullptr);
}

TEST(MeshGroupingEntity, getSideSetFromPart)
{
  if(stk::parallel_machine_size(MPI_COMM_WORLD) > 1) { return; }
  std::shared_ptr<stk::mesh::BulkData> bulk = build_mesh(MPI_COMM_WORLD);
  stk::mesh::MetaData& meta = bulk->mesh_meta_data();

  stk::io::StkMeshIoBroker ioBroker;
  stk::io::fill_mesh_preexisting(ioBroker, "generated:1x1x1|sideset:x", *bulk);
  std::shared_ptr<Ioss::Region> ioRegion = ioBroker.get_input_ioss_region();

  stk::mesh::Part* part = meta.get_part("surface_1");
  EXPECT_TRUE(part != nullptr);
  EXPECT_EQ(meta.side_rank(), part->primary_entity_rank());

  Ioss::GroupingEntity* entity = stk::io::get_grouping_entity(*ioRegion, *part);
  EXPECT_TRUE(entity != nullptr);
  EXPECT_EQ(Ioss::SIDESET, entity->type());
}

TEST(MeshGroupingEntity, getNodeSetFromPart)
{
  if(stk::parallel_machine_size(MPI_COMM_WORLD) > 1) { return; }
  std::shared_ptr<stk::mesh::BulkData> bulk = build_mesh(MPI_COMM_WORLD);
  stk::mesh::MetaData& meta = bulk->mesh_meta_data();

  stk::io::StkMeshIoBroker ioBroker;
  stk::io::fill_mesh_preexisting(ioBroker, "generated:1x1x1|nodeset:x", *bulk);
  std::shared_ptr<Ioss::Region> ioRegion = ioBroker.get_input_ioss_region();

  stk::mesh::Part* part = meta.get_part("nodelist_1");
  EXPECT_TRUE(part != nullptr);
  EXPECT_EQ(stk::topology::NODE_RANK, part->primary_entity_rank());

  Ioss::GroupingEntity* entity = stk::io::get_grouping_entity(*ioRegion, *part);
  EXPECT_TRUE(entity != nullptr);
  EXPECT_EQ(Ioss::NODESET, entity->type());
}

TEST(MeshGroupingEntity, nonInputIoPart)
{
  if(stk::parallel_machine_size(MPI_COMM_WORLD) > 1) { return; }
  const unsigned spatialDim = 3;
  std::shared_ptr<stk::mesh::BulkData> bulk = build_mesh(spatialDim, MPI_COMM_WORLD);
  stk::mesh::MetaData& meta = bulk->mesh_meta_data();

  stk::io::StkMeshIoBroker ioBroker;
  stk::mesh::Part& dummyPart = meta.declare_part("dummy", stk::topology::ELEMENT_RANK);
  EXPECT_FALSE(stk::io::is_part_io_part(dummyPart));
  stk::io::put_io_part_attribute(dummyPart);
  EXPECT_TRUE(stk::io::is_part_io_part(dummyPart));

  stk::io::fill_mesh_preexisting(ioBroker, "generated:1x1x1", *bulk);
  std::shared_ptr<Ioss::Region> ioRegion = ioBroker.get_input_ioss_region();

  Ioss::GroupingEntity* entity = stk::io::get_grouping_entity(*ioRegion, dummyPart);
  EXPECT_TRUE(entity == nullptr);
}

void create_2block_quad_mesh(const std::string& fileName)
{
  const unsigned spatialDim = 2;
  std::shared_ptr<stk::mesh::BulkData> bulk = build_mesh(spatialDim, MPI_COMM_WORLD);
  std::string meshDesc = "0,1,QUAD_4_2D,1,2,5,4,block_1\n"
                         "0,2,QUAD_4_2D,2,3,6,5,block_2";

  stk::unit_test_util::setup_text_mesh(*bulk, meshDesc);
  stk::io::write_mesh(fileName, *bulk);
}

TEST(MeshGroupingEntity, matchhingNameAndType)
{
  if(stk::parallel_machine_size(MPI_COMM_WORLD) > 1) { return; }
  const unsigned spatialDim = 2;
  std::shared_ptr<stk::mesh::BulkData> bulk = build_mesh(spatialDim, MPI_COMM_WORLD);
  stk::mesh::MetaData& meta = bulk->mesh_meta_data();

  stk::io::StkMeshIoBroker ioBroker;
  std::string fileName = "output.g";
  std::string skippedPartName = "block_2";

  create_2block_quad_mesh(fileName);

  ioBroker.set_bulk_data(bulk);
  size_t inputIndex = ioBroker.add_mesh_database(fileName, stk::io::READ_MESH);
  Ioss::DatabaseIO* dbIo = ioBroker.get_input_database(inputIndex);
  dbIo->set_block_omissions({skippedPartName});
  ioBroker.create_input_mesh();

  stk::mesh::Part* skippedPart = meta.get_part(skippedPartName);
  EXPECT_EQ(nullptr, skippedPart);

  skippedPart = &meta.declare_part(skippedPartName, stk::topology::NODE_RANK);
  EXPECT_FALSE(stk::io::is_part_io_part(*skippedPart));
  stk::io::put_io_part_attribute(*skippedPart);
  EXPECT_TRUE(stk::io::is_part_io_part(*skippedPart));

  ioBroker.add_all_mesh_fields_as_input_fields();
  ioBroker.populate_bulk_data();

  std::shared_ptr<Ioss::Region> ioRegion = ioBroker.get_input_ioss_region();

  Ioss::GroupingEntity* entity = stk::io::get_grouping_entity(*ioRegion, *skippedPart);
  EXPECT_TRUE(entity == nullptr);
}

TEST(MeshGroupingEntity, iossEntityTypesForElementRank)
{
  stk::mesh::MetaData meta(2);

  std::vector<Ioss::EntityType> entityTypeVec = stk::io::get_ioss_entity_types(meta, stk::topology::ELEMENT_RANK);

  EXPECT_EQ(2u, entityTypeVec.size());
  EXPECT_TRUE(entityTypeVec[0] == Ioss::ELEMENTBLOCK || entityTypeVec[1] == Ioss::ELEMENTBLOCK);
  EXPECT_TRUE(entityTypeVec[0] == Ioss::SUPERELEMENT || entityTypeVec[1] == Ioss::SUPERELEMENT);
}

TEST(MeshGroupingEntity, iossEntityTypesForFaceRank)
{
  stk::mesh::MetaData meta2D(2);
  stk::mesh::MetaData meta3D(3);

  std::vector<Ioss::EntityType> entityTypeVec = stk::io::get_ioss_entity_types(meta2D, stk::topology::FACE_RANK);

  EXPECT_EQ(0u, entityTypeVec.size());

  entityTypeVec = stk::io::get_ioss_entity_types(meta2D, stk::topology::EDGE_RANK);

  EXPECT_EQ(2u, entityTypeVec.size());
  EXPECT_TRUE(entityTypeVec[0] == Ioss::SIDESET || entityTypeVec[1] == Ioss::SIDESET);
  EXPECT_TRUE(entityTypeVec[0] == Ioss::SIDEBLOCK || entityTypeVec[1] == Ioss::SIDEBLOCK);

  entityTypeVec = stk::io::get_ioss_entity_types(meta3D, stk::topology::FACE_RANK);

  EXPECT_EQ(2u, entityTypeVec.size());
  EXPECT_TRUE(entityTypeVec[0] == Ioss::SIDESET || entityTypeVec[1] == Ioss::SIDESET);
  EXPECT_TRUE(entityTypeVec[0] == Ioss::SIDEBLOCK || entityTypeVec[1] == Ioss::SIDEBLOCK);

  entityTypeVec = stk::io::get_ioss_entity_types(meta3D, stk::topology::EDGE_RANK);
  EXPECT_EQ(0u, entityTypeVec.size());
}

TEST(MeshGroupingEntity, iossEntityTypesForConstraintRank)
{
  stk::mesh::MetaData meta(2);

  std::vector<Ioss::EntityType> entityTypeVec = stk::io::get_ioss_entity_types(meta, stk::topology::CONSTRAINT_RANK);

  EXPECT_EQ(0u, entityTypeVec.size());
}

TEST(MeshGroupingEntity, iossEntityTypesForNodeRank)
{
  stk::mesh::MetaData meta;

  std::vector<Ioss::EntityType> entityTypeVec = stk::io::get_ioss_entity_types(meta, stk::topology::NODE_RANK);

  EXPECT_EQ(2u, entityTypeVec.size());
  EXPECT_TRUE(entityTypeVec[0] == Ioss::NODESET || entityTypeVec[1] == Ioss::NODESET);
  EXPECT_TRUE(entityTypeVec[0] == Ioss::NODEBLOCK || entityTypeVec[1] == Ioss::NODEBLOCK);
}
