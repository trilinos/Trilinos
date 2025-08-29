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
#include <stk_mesh/base/GetEntities.hpp>  // for get_selected_entities
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
#include <stk_unit_test_utils/BuildMesh.hpp>

namespace {

using stk::unit_test_util::build_mesh;
using stk::unit_test_util::build_mesh_no_simple_fields;

void create_one_element_mesh_with_nodeset(stk::mesh::BulkData& bulk, const std::string & /*filename*/)
{
  if(bulk.parallel_size() == 1) {
    stk::mesh::MetaData& meta = bulk.mesh_meta_data();
    unsigned spatialDim = meta.spatial_dimension();

    stk::mesh::Part * hexPart   = &meta.declare_part_with_topology("block_1",   stk::topology::HEX_8);
    stk::io::put_io_part_attribute(*hexPart);

    stk::mesh::Part * nodesetPart   = &meta.declare_part_with_topology("nodelist_1",   stk::topology::NODE);
    stk::io::put_io_part_attribute(*nodesetPart);

    stk::mesh::Field<double> & coordsField = meta.declare_field<double>(stk::topology::NODE_RANK, "coordinates", 1);
    stk::mesh::put_field_on_mesh(coordsField, meta.universal_part(), spatialDim, nullptr);

    stk::mesh::Field<double> & nodeSetField = meta.declare_field<double>(stk::topology::NODE_RANK, "nodeSetField", 1);
    stk::mesh::put_field_on_mesh(nodeSetField, *nodesetPart, 1, nullptr);

    meta.commit();

    std::vector<stk::mesh::EntityIdVector> hexNodeIDs {
      { 1, 2, 3, 4, 5, 6, 7, 8 }
    };
    stk::mesh::EntityId hexElemIDs[] = { 1 };

    // Build the base hex mesh
    bulk.modification_begin();
    for (size_t i = 0; i < hexNodeIDs.size(); ++i) {
      stk::mesh::declare_element(bulk, *hexPart, hexElemIDs[i], hexNodeIDs[i]);
    }
    bulk.modification_end();

    stk::mesh::Entity node_1 = bulk.get_entity(stk::topology::NODE_RANK, 1);
    EXPECT_TRUE(bulk.is_valid(node_1));

    bulk.modification_begin();
    bulk.change_entity_parts(node_1, stk::mesh::PartVector{nodesetPart}, stk::mesh::PartVector{});

    stk::mesh::EntityVector nodes;
    stk::mesh::get_entities(bulk, stk::topology::NODE_RANK, nodes);
    EXPECT_EQ(8u, nodes.size());

    std::vector<double> coordinates = {0,0,0,
                                       1,0,0,
                                       1,1,0,
                                       0,1,0,
                                       0,0,1,
                                       1,0,1,
                                       1,1,1,
                                       0,1,1};

    auto coordData = coordsField.data<stk::mesh::ReadWrite>();
    for(size_t nodeIndex=0; nodeIndex < nodes.size(); nodeIndex++) {
      auto nodalCoords = coordData.entity_values(nodes[nodeIndex]);
      for(stk::mesh::ComponentIdx coordIndex(0); coordIndex < static_cast<int>(spatialDim); ++coordIndex){
        nodalCoords(coordIndex) = coordinates[nodeIndex*spatialDim+coordIndex];
      }
    }

    stk::mesh::EntityVector elements;
    stk::mesh::get_selected_entities(*hexPart, bulk.buckets(stk::topology::ELEM_RANK), elements);
    EXPECT_EQ(1u, elements.size());

    stk::mesh::Part * part = meta.get_part("nodelist_1");
    EXPECT_TRUE(stk::io::is_part_io_part(*part));
  }
}

void create_and_write_one_hex_mesh_with_nodeset(const std::string& filename)
{
  const unsigned spatialDim = 3;
  std::shared_ptr<stk::mesh::BulkData> bulk = build_mesh(spatialDim, MPI_COMM_WORLD, stk::mesh::BulkData::NO_AUTO_AURA);

  create_one_element_mesh_with_nodeset(*bulk, filename);

  stk::io::write_mesh_with_fields(filename, *bulk, 1, 1.0);
}

TEST(StkMeshIoBroker, readNodesetWithDistributionFactor) {
  if(stk::parallel_machine_size(MPI_COMM_WORLD) != 1) { return; }

  std::string inputFile = "readNodesetWithDistributionFactor.g";
  create_and_write_one_hex_mesh_with_nodeset(inputFile);

  const unsigned spatialDim = 3;
  std::shared_ptr<stk::mesh::BulkData> bulk = build_mesh(spatialDim, MPI_COMM_WORLD, stk::mesh::BulkData::NO_AUTO_AURA);
  stk::mesh::MetaData& meta = bulk->mesh_meta_data();

  stk::io::StkMeshIoBroker ioBroker;
  stk::io::fill_mesh_preexisting(ioBroker, inputFile, *bulk);

  stk::mesh::Part* nodeSet = meta.get_part("nodelist_1");
  EXPECT_TRUE(nodeSet != nullptr);

  stk::mesh::FieldBase* nodeSetDf = meta.get_field(stk::topology::NODE_RANK, "distribution_factors_nodelist_1");
  EXPECT_TRUE(nodeSetDf != nullptr);
  unlink(inputFile.c_str());
}

TEST(StkMeshIoBroker, readNodesetWithoutDistributionFactor) {
  if(stk::parallel_machine_size(MPI_COMM_WORLD) != 1) { return; }

  std::string inputFile = "readNodesetWithoutDistributionFactor.g";
  create_and_write_one_hex_mesh_with_nodeset(inputFile);

  const unsigned spatialDim = 3;
  std::shared_ptr<stk::mesh::BulkData> bulk = build_mesh(spatialDim, MPI_COMM_WORLD, stk::mesh::BulkData::NO_AUTO_AURA);
  stk::mesh::MetaData& meta = bulk->mesh_meta_data();

  stk::io::StkMeshIoBroker ioBroker;
  ioBroker.set_auto_load_distribution_factor_per_nodeset(false);
  stk::io::fill_mesh_preexisting(ioBroker, inputFile, *bulk);

  stk::mesh::Part* nodeSet = meta.get_part("nodelist_1");
  EXPECT_TRUE(nodeSet != nullptr);

  stk::mesh::FieldBase* nodeSetDf = meta.get_field(stk::topology::NODE_RANK, "distribution_factors_nodelist_1");
  EXPECT_TRUE(nodeSetDf == nullptr);
  unlink(inputFile.c_str());
}

void setup_field_data(const std::string& inputFile, const std::string& nodeSetFieldName, stk::mesh::BulkData& bulk)
{
  stk::mesh::MetaData& meta = bulk.mesh_meta_data();
  create_and_write_one_hex_mesh_with_nodeset(inputFile);

  meta.declare_part_with_topology("nodelist_1", stk::topology::NODE);
  meta.declare_field<double>(stk::topology::NODE_RANK, nodeSetFieldName, 1);
}

void load_field_data(const std::string& inputFile, stk::mesh::BulkData& bulk, stk::io::StkMeshIoBroker& ioBroker, stk::io::MeshField& meshField)
{
  ioBroker.set_bulk_data(bulk);
  ioBroker.add_mesh_database(inputFile, stk::io::READ_MESH);
  ioBroker.create_input_mesh();
  ioBroker.populate_bulk_data();

  ioBroker.add_input_field(meshField);
}

void test_field_data_no_throw(const std::string& inputFile, stk::mesh::BulkData& bulk, stk::io::MeshField& meshField)
{
  stk::io::StkMeshIoBroker ioBroker;
  load_field_data(inputFile, bulk, ioBroker, meshField);

  EXPECT_NO_THROW(ioBroker.read_defined_input_fields(0.0));
  unlink(inputFile.c_str());
}

void test_field_data_throw(const std::string& inputFile, stk::mesh::BulkData& bulk, stk::io::MeshField& meshField)
{
  stk::io::StkMeshIoBroker ioBroker;
  load_field_data(inputFile, bulk, ioBroker, meshField);

  EXPECT_THROW(ioBroker.read_defined_input_fields(0.0), std::runtime_error);
  unlink(inputFile.c_str());
}

TEST(StkMeshIoBroker, readSubsetFieldData) {
  if(stk::parallel_machine_size(MPI_COMM_WORLD) != 1) { return; }

  std::string inputFile = "readSubsetFieldData.g";
  std::string nodeSetString = "nodeSetField";

  const unsigned spatialDim = 3;
  std::shared_ptr<stk::mesh::BulkData> bulk = build_mesh(spatialDim, MPI_COMM_WORLD, stk::mesh::BulkData::NO_AUTO_AURA);
  stk::mesh::MetaData& meta = bulk->mesh_meta_data();

  setup_field_data(inputFile, nodeSetString, *bulk);

  stk::mesh::FieldBase* nodeSetField = meta.get_field(stk::topology::NODE_RANK, nodeSetString);
  stk::mesh::Part* nodesetPart = meta.get_part("nodelist_1");

  stk::mesh::put_field_on_mesh(*nodeSetField, *nodesetPart, 1, nullptr);

  stk::io::MeshField meshField(nodeSetField, nodeSetString);
  meshField.add_subset(*nodesetPart);

  test_field_data_no_throw(inputFile, *bulk, meshField);
}

TEST(StkMeshIoBroker, readFieldDataOnUniversalSetButNotDefinedOnSubset) {
  if(stk::parallel_machine_size(MPI_COMM_WORLD) != 1) { return; }

  std::string inputFile = "readFieldDataOnUniversalSetButNotDefinedOnSubset.g";
  std::string nodeSetString = "nodeSetField";

  const unsigned spatialDim = 3;
  std::shared_ptr<stk::mesh::BulkData> bulk = build_mesh(spatialDim, MPI_COMM_WORLD, stk::mesh::BulkData::NO_AUTO_AURA);
  stk::mesh::MetaData& meta = bulk->mesh_meta_data();

  setup_field_data(inputFile, nodeSetString, *bulk);

  stk::mesh::FieldBase* nodeSetField = meta.get_field(stk::topology::NODE_RANK, nodeSetString);
  stk::mesh::put_field_on_mesh(*nodeSetField, meta.universal_part(), 1, nullptr);

  stk::io::MeshField meshField(nodeSetField, nodeSetString);

  test_field_data_throw(inputFile, *bulk, meshField);
}

TEST(StkMeshIoBroker, readFieldDataOnUniversalSetAndDefinedOnSubset) {
  if(stk::parallel_machine_size(MPI_COMM_WORLD) != 1) { return; }

  std::string inputFile = "readFieldDataOnUniversalSetAndDefinedOnSubset.g";
  std::string nodeSetString = "nodeSetField";

  const unsigned spatialDim = 3;
  std::shared_ptr<stk::mesh::BulkData> bulk = build_mesh(spatialDim, MPI_COMM_WORLD, stk::mesh::BulkData::NO_AUTO_AURA);
  stk::mesh::MetaData& meta = bulk->mesh_meta_data();

  setup_field_data(inputFile, nodeSetString, *bulk);

  stk::mesh::FieldBase* nodeSetField = meta.get_field(stk::topology::NODE_RANK, nodeSetString);
  stk::mesh::Part* nodesetPart = meta.get_part("nodelist_1");
  stk::mesh::put_field_on_mesh(*nodeSetField, meta.universal_part(), 1, nullptr);

  stk::io::MeshField meshField(nodeSetField, nodeSetString);
  meshField.add_subset(*nodesetPart);

  test_field_data_no_throw(inputFile, *bulk, meshField);
}

}
