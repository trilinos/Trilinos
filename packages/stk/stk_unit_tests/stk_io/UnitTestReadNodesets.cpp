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

namespace {

void create_one_element_mesh_with_nodeset(stk::mesh::BulkData& bulk, const std::string & filename)
{
  if(bulk.parallel_size() == 1) {
    stk::mesh::MetaData& meta = bulk.mesh_meta_data();
    unsigned spatialDim = meta.spatial_dimension();

    stk::mesh::Part * hexPart   = &meta.declare_part_with_topology("block_1",   stk::topology::HEX_8);
    stk::io::put_io_part_attribute(*hexPart);

    stk::mesh::Part * nodesetPart   = &meta.declare_part_with_topology("nodelist_1",   stk::topology::NODE);
    stk::io::put_io_part_attribute(*nodesetPart);

    stk::mesh::Field<double, stk::mesh::Cartesian> & coordsField = meta.declare_field<stk::mesh::Field<double, stk::mesh::Cartesian>>(stk::topology::NODE_RANK, "coordinates", 1);
    stk::mesh::put_field_on_mesh(coordsField, meta.universal_part(), spatialDim, static_cast<double*>(nullptr));

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

    for(size_t nodeIndex=0; nodeIndex < nodes.size(); nodeIndex++)
    {
      double * nodalCoords = stk::mesh::field_data(coordsField, nodes[nodeIndex]);
      for(unsigned coordIndex=0; coordIndex < spatialDim; coordIndex++)
         nodalCoords[coordIndex] = coordinates[nodeIndex*spatialDim+coordIndex];
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
  stk::mesh::MetaData meta(spatialDim);
  stk::mesh::BulkData bulk(meta, MPI_COMM_WORLD, stk::mesh::BulkData::NO_AUTO_AURA);

  create_one_element_mesh_with_nodeset(bulk, filename);

  stk::io::write_mesh(filename, bulk);
}

TEST(StkMeshIoBroker, readNodesetWithDistributionFactor) {
  if(stk::parallel_machine_size(MPI_COMM_WORLD) != 1) { return; }

  std::string inputFile = "readNodesetWithDistributionFactor.g";
  create_and_write_one_hex_mesh_with_nodeset(inputFile);

  const unsigned spatialDim = 3;
  stk::mesh::MetaData meta(spatialDim);
  stk::mesh::BulkData bulk(meta, MPI_COMM_WORLD, stk::mesh::BulkData::NO_AUTO_AURA);

  stk::io::StkMeshIoBroker ioBroker;
  stk::io::fill_mesh_preexisting(ioBroker, inputFile, bulk);

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
  stk::mesh::MetaData meta(spatialDim);
  stk::mesh::BulkData bulk(meta, MPI_COMM_WORLD, stk::mesh::BulkData::NO_AUTO_AURA);

  stk::io::StkMeshIoBroker ioBroker;
  ioBroker.set_auto_load_distribution_factor_per_nodeset(false);
  stk::io::fill_mesh_preexisting(ioBroker, inputFile, bulk);

  stk::mesh::Part* nodeSet = meta.get_part("nodelist_1");
  EXPECT_TRUE(nodeSet != nullptr);

  stk::mesh::FieldBase* nodeSetDf = meta.get_field(stk::topology::NODE_RANK, "distribution_factors_nodelist_1");
  EXPECT_TRUE(nodeSetDf == nullptr);
  unlink(inputFile.c_str());
}
}
