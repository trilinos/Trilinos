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
#include <stk_io/FillMesh.hpp>
#include <stk_io/WriteMesh.hpp>
#include <stk_mesh/base/BulkData.hpp>   // for BulkData
#include <stk_mesh/base/MetaData.hpp>   // for MetaData
#include <stk_unit_test_utils/BuildMesh.hpp>
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
#include "stk_mesh/base/Comm.hpp"
#include "stk_unit_test_utils/ioUtils.hpp"
#include "stk_unit_test_utils/getOption.h"
#include "Ioss_Region.h"
#include "Ioss_ElementBlock.h"
#include "Ioss_NodeSet.h"

namespace {

using stk::unit_test_util::build_mesh;
using stk::unit_test_util::build_mesh_no_simple_fields;

//BeginDocTest1
TEST(StkMeshIoBroker, outputEqualsInput)
{
    // A simple test for reading and writing an exodus file, to make sure
    // that the output file is the same as the input file.

    std::string input_filename = stk::unit_test_util::get_option("--input-mesh", "no-mesh-specified");

    if (input_filename != "no-mesh-specified")
    {
        stk::ParallelMachine comm = MPI_COMM_WORLD;
        std::shared_ptr<stk::mesh::BulkData> bulk = build_mesh(comm);

        stk::io::StkMeshIoBroker ioBroker(comm);
        ioBroker.set_bulk_data(bulk);
        ioBroker.add_mesh_database(input_filename, stk::io::READ_MESH);
        ioBroker.create_input_mesh();
        ioBroker.populate_bulk_data();

        std::vector<size_t> counts;
        stk::mesh::comm_mesh_counts(*bulk, counts);

        if (bulk->parallel_rank() == 0)
        {
            std::cerr << "input-mesh = " << input_filename << std::endl;
            std::cerr << "      num-nodes = " << counts[stk::topology::NODE_RANK]
                      << ", num-elems = " << counts[stk::topology::ELEM_RANK] << std::endl;
        }

        std::cerr << "writing mesh with same io-broker as input" << std::endl;
        std::string same_iobroker_output_filename = std::string("output-same-iobroker-")+input_filename;
        size_t output_file_index = ioBroker.create_output_mesh(same_iobroker_output_filename, stk::io::WRITE_RESULTS);
        ioBroker.write_output_mesh(output_file_index);

        std::string new_iobroker_output_filename = std::string("output-")+input_filename;
        std::cerr << "writing mesh with new io-broker" << std::endl;

        stk::io::StkMeshIoBroker new_iobroker(comm);

        new_iobroker.set_bulk_data(bulk);
        size_t new_iobroker_output_file_index = new_iobroker.create_output_mesh(new_iobroker_output_filename, stk::io::WRITE_RESULTS);
        new_iobroker.write_output_mesh(new_iobroker_output_file_index);

        //Ideally we could use exodiff to compare the two output meshes to
        //verify they are the same.
        //But exodiff won't say two mesh files are different if the difference is only
        //a block-ordering difference. So let's instead verify that the two ioss-regions
        //have the same element-block grouping-entities.
        const Ioss::ElementBlockContainer& input_broker_elem_blocks = ioBroker.get_output_ioss_region(output_file_index)->get_element_blocks();
        const Ioss::ElementBlockContainer& new_iobroker_elem_blocks = new_iobroker.get_output_ioss_region(new_iobroker_output_file_index)->get_element_blocks();

        EXPECT_EQ(input_broker_elem_blocks.size(), new_iobroker_elem_blocks.size());
        for(size_t i=0; i<input_broker_elem_blocks.size(); ++i)
        {
            std::string input_name = input_broker_elem_blocks[i]->name();
            std::string new_name = new_iobroker_elem_blocks[i]->name();
            EXPECT_EQ(input_name, new_name);
        }

        //now also compare names/ids of nodesets
        const Ioss::NodeSetContainer input_broker_node_sets = ioBroker.get_output_ioss_region(output_file_index)->get_nodesets();
        const Ioss::NodeSetContainer new_iobroker_node_sets = new_iobroker.get_output_ioss_region(output_file_index)->get_nodesets();
        EXPECT_EQ(input_broker_node_sets.size(), new_iobroker_node_sets.size());
        for(size_t i=0; i<input_broker_node_sets.size(); ++i)
        {
            std::string input_name = input_broker_node_sets[i]->name();
            std::string new_name = new_iobroker_node_sets[i]->name();
            EXPECT_EQ(input_name, new_name);
            int input_id = input_broker_node_sets[i]->get_property("id").get_int();
            int new_id = new_iobroker_node_sets[i]->get_property("id").get_int();
            EXPECT_EQ(input_id, new_id);
        }
    }
    else
    {
        std::cout << "test StkMeshIoBroker.outputEqualsInput, no mesh specified" << std::endl;
    }
}
//BeginDocTest1

void assert_read_mesh_has_original_topology(const std::string& fileName,
                                            const std::string& expectedTopologyType)
{
  const unsigned spatialDim = 3;
  std::shared_ptr<stk::mesh::BulkData> bulk = build_mesh(spatialDim, MPI_COMM_WORLD);
  stk::io::fill_mesh(fileName, *bulk);

  stk::mesh::Part* block1 = bulk->mesh_meta_data().get_part("block_1");
  ASSERT_TRUE(block1 != nullptr);
  ASSERT_TRUE(stk::io::is_part_io_part(*block1));
  ASSERT_TRUE(stk::io::has_original_topology_type(*block1));
  ASSERT_EQ(expectedTopologyType, stk::io::get_original_topology_type(*block1));
}

void read_then_write_mesh(const std::string& inputFileName,
                          const std::string& outputFileName)
{
  const unsigned spatialDim = 3;
  std::shared_ptr<stk::mesh::BulkData> bulk = build_mesh(spatialDim, MPI_COMM_WORLD);
  stk::io::fill_mesh(inputFileName, *bulk);
  stk::io::write_mesh(outputFileName, *bulk);
}

TEST(StkMeshIoBroker, originalTopologyType)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) { return; }

  std::string fileName = stk::unit_test_util::get_option("--mesh", "no-mesh-specified");
  if (fileName == "no-mesh-specified") { return; }

  const bool isBar = (fileName == "bar.exo");
  const bool isBeam = (fileName == "beam.exo");
  ASSERT_TRUE(isBar || isBeam);
  std::string expectedTopologyType("beam2");
  if (isBar) {
    expectedTopologyType = "bar2";
  }

  assert_read_mesh_has_original_topology(fileName, expectedTopologyType);

  std::string tmpFileName("tmpMesh.exo");
  read_then_write_mesh(fileName, tmpFileName);

  assert_read_mesh_has_original_topology(tmpFileName, expectedTopologyType);

  unlink(tmpFileName.c_str());
}

void read_mesh_with_omitted_block(std::shared_ptr<stk::mesh::BulkData> bulk,
                                   const std::string& fileName,
                                   bool createEmptyOmittedBlock)
{
  stk::io::StkMeshIoBroker broker;
  const std::vector<std::string> omittedBlocks{"block_1"};

  broker.set_bulk_data(bulk);
  broker.set_create_empty_block_for_omitted_block(createEmptyOmittedBlock);
  size_t inputIndex = broker.add_mesh_database(fileName, stk::io::READ_MESH);
  Ioss::DatabaseIO* dbIo = broker.get_input_database(inputIndex);
  dbIo->set_block_omissions(omittedBlocks);
  broker.create_input_mesh();
  broker.add_all_mesh_fields_as_input_fields();
  broker.populate_bulk_data();
  int numSteps = broker.get_num_time_steps();
  if(numSteps>0) {
    broker.read_defined_input_fields(numSteps);
  }
}

void test_read_with_omitted_block(std::shared_ptr<stk::mesh::BulkData> bulk,
                                   const std::string& fileName,
                                   const std::string& elemFieldName,
                                   bool createEmptyOmittedBlock)
{
  read_mesh_with_omitted_block(bulk, fileName, createEmptyOmittedBlock);

  const stk::mesh::MetaData& meta = bulk->mesh_meta_data();

  const stk::mesh::FieldBase* coordinates = meta.coordinate_field();
  EXPECT_TRUE(nullptr != coordinates);

  const stk::mesh::FieldBase* elemField = meta.get_field(stk::topology::ELEM_RANK, elemFieldName);
  EXPECT_TRUE(nullptr != elemField);

  stk::mesh::Part* block1 = meta.get_part("block_1");
  if(createEmptyOmittedBlock) {
    EXPECT_TRUE(block1 != nullptr);
    EXPECT_FALSE(stk::io::is_part_io_part(*block1));
    EXPECT_TRUE(coordinates->defined_on(*block1));
    EXPECT_FALSE(elemField->defined_on(*block1));
    EXPECT_EQ(stk::topology::HEX_8, block1->topology());
    const stk::mesh::BucketVector& block1Buckets = bulk->get_buckets(stk::topology::ELEM_RANK, *block1);
    EXPECT_EQ(0u, block1Buckets.size());
  } else {
    EXPECT_TRUE(block1 == nullptr);
  }

  stk::mesh::Part* block2 = meta.get_part("block_2");
  EXPECT_TRUE(block2 != nullptr);
  EXPECT_TRUE(stk::io::is_part_io_part(*block2));
  EXPECT_TRUE(coordinates->defined_on(*block2));
  EXPECT_TRUE(elemField->defined_on(*block2));
  EXPECT_EQ(stk::topology::HEX_8, block2->topology());
  const stk::mesh::BucketVector& block2Buckets = bulk->get_buckets(stk::topology::ELEM_RANK, *block2);
  EXPECT_EQ(1u, block2Buckets.size());

  stk::mesh::Part* surface1 = meta.get_part("surface_1");
  EXPECT_TRUE(surface1 != nullptr);
  EXPECT_TRUE(bulk->does_sideset_exist(*surface1));
  const stk::mesh::SideSet& sideset = bulk->get_sideset(*surface1);
  EXPECT_EQ(1u, sideset.size());
  stk::mesh::Selector selector(meta.universal_part());
  stk::mesh::EntityVector faces;
  stk::mesh::get_entities(*bulk, meta.side_rank(), selector, faces);
  EXPECT_EQ(1u, faces.size());
}

void create_file_with_element_field(MPI_Comm communicator,
                                    unsigned spatialDim,
                                    const std::string& meshDesc,
                                    const std::string& fileName,
                                    const std::string& elemFieldName)
{
  stk::mesh::MeshBuilder builder(communicator);
  builder.set_spatial_dimension(spatialDim);
  std::shared_ptr<stk::mesh::BulkData> outBulk = builder.create();

  stk::mesh::MetaData& meta = outBulk->mesh_meta_data();
  stk::mesh::Field<double> & field = meta.declare_field<double>(stk::topology::ELEM_RANK, elemFieldName);
  stk::mesh::put_field_on_mesh(field, meta.universal_part(), nullptr);

  stk::io::fill_mesh(meshDesc, *outBulk);
  int step = 1;
  double time = 1.0;
  stk::io::write_mesh_with_specified_fields(fileName, *outBulk, std::vector<std::string>{elemFieldName}, step, time);
}

void run_mesh_with_omitted_block(MPI_Comm communicator, const std::string& fileName, bool createEmptyOmittedBlock)
{
  const unsigned spatialDim = 3;

  stk::mesh::MeshBuilder builder(communicator);
  builder.set_spatial_dimension(spatialDim);
  std::shared_ptr<stk::mesh::BulkData> bulk = builder.create();

  // 2 hexes in separate blocks with spanning sideset across both elements
  std::string meshDesc =
      "textmesh:\
       0,1,HEX_8, 1, 2, 3, 4, 5, 6, 7, 8,block_1\n\
       0,2,HEX_8, 5, 6, 7, 8, 9,10,11,12,block_2|sideset:data=1,1, 2,1";

  const std::string elemFieldName{"elemField"};

  create_file_with_element_field(communicator, spatialDim, meshDesc, fileName, elemFieldName);

  test_read_with_omitted_block(bulk, fileName, elemFieldName, createEmptyOmittedBlock);

  unlink(fileName.c_str());
}

TEST(StkMeshIoBroker, createMeshWithOmittedBlock_withEmptyBlock)
{
  MPI_Comm communicator = MPI_COMM_WORLD;
  size_t nProc = stk::parallel_machine_size(communicator);

  if(nProc > 1) {
    GTEST_SKIP();
  }

  bool createEmptyOmittedBlock = true;
  const std::string fileName{"emptyOmittedBlock.g"};
  run_mesh_with_omitted_block(communicator, fileName, createEmptyOmittedBlock);
}

TEST(StkMeshIoBroker, createMeshWithOmittedBlock_withNoEmptyBlock)
{
  MPI_Comm communicator = MPI_COMM_WORLD;
  size_t nProc = stk::parallel_machine_size(communicator);

  if(nProc > 1) {
    GTEST_SKIP();
  }

  bool createEmptyOmittedBlock = false;
  const std::string fileName{"noEmptyOmittedBlock.g"};
  run_mesh_with_omitted_block(communicator, fileName, createEmptyOmittedBlock);
}

}
