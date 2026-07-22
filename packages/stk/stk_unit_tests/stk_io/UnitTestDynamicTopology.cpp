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
#include "gtest/gtest.h"                // for AssertHelper, ASSERT_TRUE

#include <stddef.h>                     // for size_t
#include <string>                       // for string
#include <vector>                       // for vector, etc
#include <numeric>
#include <algorithm>
#include <fstream>
#include <functional>
#include <iostream>
#include <unistd.h> // for unlink

#include "stk_io/DatabasePurpose.hpp"   // for DatabasePurpose::READ_MESH, etc
#include <stk_io/IossBridge.hpp>        // for is_part_io_part
#include <stk_io/StkMeshIoBroker.hpp>   // for StkMeshIoBroker
#include <stk_io/FillMesh.hpp>
#include <stk_io/WriteMesh.hpp>

#include <stk_mesh/base/BulkData.hpp>   // for BulkData
#include <stk_mesh/base/GetEntities.hpp>  // for get_selected_entities
#include <stk_mesh/base/MetaData.hpp>   // for MetaData
#include "stk_mesh/base/Entity.hpp"     // for Entity
#include "stk_mesh/base/Part.hpp"       // for Part
#include "stk_mesh/base/Selector.hpp"   // for Selector, operator&
#include "stk_mesh/base/Types.hpp"      // for PartVector, BucketVector, etc
#include "stk_mesh/base/FEMHelpers.hpp"
#include "stk_mesh/base/Field.hpp"

#include "stk_topology/topology.hpp"    // for topology, etc

#include "stk_util/parallel/Parallel.hpp"  // for ParallelMachine
#include <stk_unit_test_utils/GetMeshSpec.hpp>
#include <stk_unit_test_utils/BuildMesh.hpp>
#include <stk_unit_test_utils/TextMesh.hpp>

#include <Ionit_Initializer.h>

namespace {

auto setup_mesh(const std::string & meshSpec,
                stk::ParallelMachine comm,
                stk::mesh::BulkData::AutomaticAuraOption auraOption = stk::mesh::BulkData::NO_AUTO_AURA,
                unsigned spatialDim = 3,
                unsigned initialBucketCapacity = stk::mesh::get_default_initial_bucket_capacity(),
                unsigned maximumBucketCapacity = stk::mesh::get_default_maximum_bucket_capacity())
{
  stk::mesh::MeshBuilder builder(comm);
  builder.set_spatial_dimension(spatialDim);
  builder.set_aura_option(auraOption);
  builder.set_initial_bucket_capacity(initialBucketCapacity);
  builder.set_maximum_bucket_capacity(maximumBucketCapacity);

  std::shared_ptr<stk::mesh::BulkData> bulk = builder.create();

  stk::mesh::Field<int> & field = bulk->mesh_meta_data().declare_field<int>(stk::topology::NODE_RANK, "nodal_field");
  stk::io::set_field_role(field, Ioss::Field::TRANSIENT);
  stk::mesh::put_field_on_mesh(field, bulk->mesh_meta_data().universal_part(), nullptr);

  stk::io::fill_mesh("textmesh:"+meshSpec, *bulk);

  return bulk;
}

std::string get_output_filename(const std::string &fileBase, const int modCount = 1)
{
  Ioss::ParallelUtils util(Ioss::ParallelUtils::comm_world());
  std::string file = fileBase;

  if(modCount > 1) {
    std::ostringstream oss;
    oss << fileBase;
    oss << "-s" << std::setw(4) << std::setfill('0') << modCount;
    file = oss.str();
  }

  return file;
}

std::string get_decoded_output_filename(const std::string &fileBase, const int modCount = 1)
{
  Ioss::ParallelUtils util(Ioss::ParallelUtils::comm_world());
  std::string file = get_output_filename(fileBase, modCount);
  return Ioss::Utils::decode_filename(file, util.parallel_rank(), util.parallel_size());
}

void cleanup_files(const std::string &outFile, const int numOutputs = 1)
{
  for(int i=1; i<=numOutputs; i++) {
    std::string file = get_decoded_output_filename(outFile, i);
    unlink(file.c_str());
  }
}

void setup_dynamic_topology(stk::ParallelMachine comm,
                            const std::string& outFile,
                            const int numSteps,
                            const stk::io::FileOption fileOption)
{
  int numBlocks = stk::parallel_machine_size(comm);

  cleanup_files(outFile, numSteps);

  std::vector<double> coords = stk::unit_test_util::get_many_block_coordinates(numBlocks);
  std::string meshDesc = stk::unit_test_util::get_many_block_mesh_desc(numBlocks, stk::parallel_machine_size(comm));
  std::string fullMeshDesc = stk::unit_test_util::get_full_text_mesh_desc(meshDesc, coords);
  std::shared_ptr<stk::mesh::BulkData> outBulk = setup_mesh(fullMeshDesc, comm);

  stk::io::StkMeshIoBroker outBroker;
  outBroker.set_bulk_data(outBulk);
  outBroker.enable_dynamic_topology(fileOption);

  size_t outputFileIndex = outBroker.create_output_mesh(outFile, stk::io::WRITE_RESULTS);
  EXPECT_EQ(Ioss::TOPOLOGY_SAME, outBroker.get_topology_modification(outputFileIndex));

  const stk::mesh::FieldVector fields = outBulk->mesh_meta_data().get_fields();
  for(stk::mesh::FieldBase* field : fields)
  {
    const Ioss::Field::RoleType* fieldRole = stk::io::get_field_role(*field);
    if(fieldRole == nullptr || *fieldRole == Ioss::Field::TRANSIENT)
      outBroker.add_field(outputFileIndex, *field);
  }

  outBroker.write_output_mesh(outputFileIndex);

  for(int i=0; i<numSteps; i++) {
    double time = i;

    outBroker.begin_output_step(outputFileIndex, time);
    outBroker.write_defined_output_fields(outputFileIndex);
    outBroker.end_output_step(outputFileIndex);
    EXPECT_EQ(Ioss::TOPOLOGY_SAME, outBroker.get_topology_modification(outputFileIndex));
    outBroker.set_topology_modification(outputFileIndex, Ioss::TOPOLOGY_UNKNOWN);
  }

  outBroker.close_output_mesh(outputFileIndex);
}

void test_legacy_file(stk::ParallelMachine comm,
                      const std::string& outFile,
                      const int numSteps)
{
  stk::mesh::MeshBuilder builder(comm);
  std::shared_ptr<stk::mesh::BulkData> inBulk = builder.create();

  stk::io::StkMeshIoBroker inBroker;
  stk::io::fill_mesh_with_fields(outFile, inBroker, *inBulk, stk::io::READ_MESH);

  EXPECT_EQ(0, inBroker.num_mesh_groups());

  for(int i=1; i<=numSteps; i++) {
    double goldTime = i-1;
    double time = inBroker.read_defined_input_fields(i);
    EXPECT_NEAR(goldTime, time, 1.0e-6);
  }

  cleanup_files(outFile);
}

TEST(StkIoDynamicTopology, legacyOutputFileBehavior)
{
  stk::ParallelMachine comm = MPI_COMM_WORLD;
  std::string outFile("legacyOutputFile.g");

  int numSteps = 5;
  setup_dynamic_topology(comm, outFile, numSteps, stk::io::FileOption::NO_DYNAMIC_TOPOLOGY_FILE_CONTROL);
  test_legacy_file(comm, outFile, numSteps);
}

void test_multi_file(stk::ParallelMachine comm,
                     const std::string& outFile,
                     const int numSteps)
{
  for(int i=1; i<=numSteps; i++) {
    std::string inFile = get_output_filename(outFile, i);

    double goldTime = i-1;

    stk::mesh::MeshBuilder builder(comm);
    std::shared_ptr<stk::mesh::BulkData> inBulk = builder.create();

    stk::io::StkMeshIoBroker inBroker;
    stk::io::fill_mesh_with_fields(inFile, inBroker, *inBulk, stk::io::READ_MESH);

    EXPECT_EQ(0, inBroker.num_mesh_groups());

    int stepsPerFile = 1;
    double time = inBroker.read_defined_input_fields(stepsPerFile);
    EXPECT_NEAR(goldTime, time, 1.0e-6);
  }

  cleanup_files(outFile, numSteps);
}

TEST(StkIoDynamicTopology, multiFileOutputFileBehavior)
{
  stk::ParallelMachine comm = MPI_COMM_WORLD;
  std::string outFile("multiFileOutputFile.g");

  int numSteps = 5;
  setup_dynamic_topology(comm, outFile, numSteps, stk::io::FileOption::USE_DYNAMIC_TOPOLOGY_MULTI_FILE);
  test_multi_file(comm, outFile, numSteps);
}

void test_dynamic_file(stk::ParallelMachine comm,
                       const std::string& outFile,
                       const int numSteps)
{
  stk::mesh::MeshBuilder builder(comm);
  std::shared_ptr<stk::mesh::BulkData> inBulk = builder.create();

  stk::io::StkMeshIoBroker inBroker;
  stk::io::fill_mesh_with_fields(outFile, inBroker, *inBulk, stk::io::READ_MESH);

  EXPECT_EQ(numSteps, inBroker.num_mesh_groups());
  std::vector<std::string> names = inBroker.mesh_group_names();

  for(int i=0; i<numSteps; i++) {
    EXPECT_TRUE(inBroker.load_mesh_group(names[i]));

    int stepsPerFile = 1;
    double goldTime = i;
    double time = inBroker.read_defined_input_fields(stepsPerFile);
    EXPECT_NEAR(goldTime, time, 1.0e-6);
  }

  cleanup_files(outFile);
}

TEST(StkIoDynamicTopology, dynamicOutputFileBehavior)
{
  stk::ParallelMachine comm = MPI_COMM_WORLD;
  std::string outFile("dynamicOutputFile.g");

  int numSteps = 5;
  setup_dynamic_topology(comm, outFile, numSteps, stk::io::FileOption::USE_DYNAMIC_TOPOLOGY_GROUP_FILE);
  test_dynamic_file(comm, outFile, numSteps);
}

}
