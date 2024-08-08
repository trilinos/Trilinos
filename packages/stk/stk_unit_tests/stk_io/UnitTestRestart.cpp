

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
#include <stk_util/diag/StringUtil.hpp> // for make_lower
#include <stk_io/StkMeshIoBroker.hpp>   // for StkMeshIoBroker
#include <stk_mesh/base/Field.hpp>      // for Field
#include <stk_mesh/base/MetaData.hpp>   // for MetaData, put_field
#include "Ioss_DBUsage.h"               // for DatabaseUsage::READ_MODEL, etc
#include "Ioss_ElementTopology.h"       // for NameList
#include "Ioss_Field.h"                 // for Field, etc
#include "Ioss_IOFactory.h"             // for IOFactory
#include "Ioss_NodeBlock.h"             // for NodeBlock
#include "Ioss_Region.h"                // for Region, NodeBlockContainer
#include "Ioss_Utils.h"                 // for Utils
#include "stk_io/DatabasePurpose.hpp"   // for DatabasePurpose::READ_MESH, etc
#include "stk_io/WriteMesh.hpp"
#include "stk_topology/topology.hpp"    // for topology, etc
#include "stk_unit_test_utils/BuildMesh.hpp"
#include <stk_unit_test_utils/getOption.h>
#include <stk_unit_test_utils/TextMesh.hpp>
#include <stddef.h>
#include <unistd.h>
#include <string>
#include <algorithm>
#include <cctype>
#include <cstdio>

namespace {

void test_read_corrupt_restart(MPI_Comm communicator, const std::string& fileName,
                               const std::string& fieldName, double expectedMaxTime, double restartTime)
{
  std::shared_ptr<stk::mesh::BulkData> bulk = stk::unit_test_util::build_mesh(communicator);
  stk::mesh::MetaData& meta = bulk->mesh_meta_data();

  stk::io::StkMeshIoBroker stkIo;
  stkIo.set_bulk_data(*bulk);
  stkIo.add_mesh_database(fileName, stk::io::READ_RESTART);
  stkIo.create_input_mesh();
  stkIo.add_all_mesh_fields_as_input_fields();

  stk::mesh::FieldBase* field0 = meta.get_field(stk::topology::NODE_RANK, fieldName);
  ASSERT_TRUE(field0 != nullptr);
  stkIo.populate_bulk_data();

  stk::io::MeshField meshField(*field0, field0->name());
  stk::io::FieldReadStatus readStatus;

  meshField.set_read_time(restartTime);
  meshField.set_read_once(false);
  stkIo.read_input_field(meshField, readStatus);

  EXPECT_NE(expectedMaxTime, stkIo.get_max_time());
  EXPECT_TRUE(readStatus.possiblyCorrupt);

  bool saveFile = stk::unit_test_util::has_option("-s");
  if(!saveFile) {
    std::string parallelFilename = stkIo.get_input_ioss_region()->get_database()->decoded_filename();
    unlink(parallelFilename.c_str());
  } else {
    if(stk::parallel_machine_rank(communicator) == 0) {
      std::cout << "Saving corrupt restart files" << std::endl;
    }
  }
}

std::string get_corrupt_restart_mesh_filename(MPI_Comm communicator)
{
  std::ostringstream oss;
  oss << "generated:1x1x";
  oss << stk::parallel_machine_size(communicator);
  return oss.str();
}

void create_corrupt_restart(MPI_Comm communicator,
                            const std::string& restartFilename,
                            const std::string& internalClientFieldName,
                            const int nSteps, const int skipStep)
{
  std::string parallelFilename;

  stk::io::StkMeshIoBroker stkIo(communicator);
  const std::string exodusFileName = get_corrupt_restart_mesh_filename(communicator);
  size_t index = stkIo.add_mesh_database(exodusFileName, stk::io::READ_MESH);
  stkIo.set_active_mesh(index);
  stkIo.create_input_mesh();

  stk::mesh::MetaData &stkMeshMetaData = stkIo.meta_data();
  const int numberOfStates = 2;
  stk::mesh::Field<double> &field0 = stkMeshMetaData.declare_field<double>(stk::topology::NODE_RANK,
                                                                           internalClientFieldName,
                                                                           numberOfStates);

  stk::mesh::put_field_on_mesh(field0, stkMeshMetaData.universal_part(), nullptr);

  stkIo.populate_bulk_data();

  stkIo.property_add(Ioss::Property("FLUSH_INTERVAL", 1));
  size_t fileIndex = stkIo.create_output_mesh(restartFilename, stk::io::WRITE_RESTART);
  stkIo.add_field(fileIndex, field0);

  std::shared_ptr<Ioss::Region> region = stkIo.get_output_ioss_region(fileIndex);

  for(int i=0; i<=nSteps; i++) {
    if(stk::parallel_machine_rank(communicator) == 1 && (i == skipStep)) {
      continue;
    }

    double time = i;
    stkIo.begin_output_step(fileIndex, time);
    stkIo.write_defined_output_fields(fileIndex);
    stkIo.end_output_step(fileIndex);
  }
}

TEST(StkIO, CorruptRestart)
{
  MPI_Comm communicator = MPI_COMM_WORLD;
  if (stk::parallel_machine_size(communicator) == 1) { GTEST_SKIP(); }

#ifndef NDEBUG
  // Skip debug due to parallel collective calls we want to avoid
  GTEST_SKIP();
#endif

  std::string restartFilename = "output.rst";
  const std::string internalClientFieldName = "Field0";

  int nSteps = 5;
  int skipStep = 5;

  nSteps = stk::unit_test_util::get_command_line_option("-nSteps", nSteps);
  skipStep = stk::unit_test_util::get_command_line_option("-skipStep", skipStep);

  if(stk::parallel_machine_rank(communicator) == 0) {
    std::cout << "Creating corrupt restart files" << std::endl;
  }
  create_corrupt_restart(communicator, restartFilename, internalClientFieldName, nSteps, skipStep);

  if(stk::parallel_machine_rank(communicator) == 0) {
    std::cout << "Testing corrupt restart" << std::endl;
  }

  double expectedMaxTime = std::min(nSteps, skipStep);
  double restartTime = double(skipStep);
  test_read_corrupt_restart(communicator, restartFilename, internalClientFieldName, expectedMaxTime, restartTime);
}

TEST(StkIo, EmptyLocalBlock_beam2)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 2) { GTEST_SKIP(); }

  std::unique_ptr<stk::mesh::BulkData> bulk = stk::mesh::MeshBuilder(MPI_COMM_WORLD).set_spatial_dimension(3).create();

  const std::string meshDesc =
       "0,1,SHELL_QUAD_4, 1,4,5,2, block_1\n\
        1,2,SHELL_QUAD_4, 2,5,6,3, block_1\n\
        1,3,BEAM_2, 6,7, block_nodes";
  std::vector<double> coords = {0,0,0,  0,1,0,  0,2,0,
                                1,0,0,  1,1,0, 1,2,0,
                                0,3,0};

  stk::unit_test_util::setup_text_mesh(*bulk, stk::unit_test_util::get_full_text_mesh_desc(meshDesc, coords));

  stk::io::write_mesh("shellq4_beam.g", *bulk, stk::io::WRITE_RESTART);
  stk::parallel_machine_barrier(MPI_COMM_WORLD);

  stk::io::StkMeshIoBroker ioBroker(MPI_COMM_SELF);
  ioBroker.set_mesh_builder(std::make_shared<stk::mesh::MeshBuilder>());
  std::string pllFileName = "shellq4_beam.g.2." + std::to_string(stk::parallel_machine_rank(MPI_COMM_WORLD));

  ioBroker.add_mesh_database(pllFileName, stk::io::READ_MESH);
  EXPECT_NO_THROW(ioBroker.create_input_mesh());
  ioBroker.populate_bulk_data();

  std::remove(pllFileName.c_str());
}

}
