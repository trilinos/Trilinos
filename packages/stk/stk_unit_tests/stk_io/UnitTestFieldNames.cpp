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
#include <unistd.h>                     // for unlink
#include <stk_util/diag/StringUtil.hpp> // for make_lower
#include <stk_io/StkMeshIoBroker.hpp>   // for StkMeshIoBroker
#include <stk_mesh/base/Field.hpp>      // for Field
#include <stk_mesh/base/MetaData.hpp>   // for MetaData, put_field
#include <gtest/gtest.h>
#include <string>                       // for string, basic_string
#include <algorithm>
#include <cctype>
#include "Ioss_DBUsage.h"               // for DatabaseUsage::READ_MODEL, etc
#include "Ioss_ElementTopology.h"       // for NameList
#include "Ioss_Field.h"                 // for Field, etc
#include "Ioss_IOFactory.h"             // for IOFactory
#include "Ioss_NodeBlock.h"             // for NodeBlock
#include "Ioss_Region.h"                // for Region, NodeBlockContainer
#include "Ioss_Utils.h"                 // for Utils
#include "stk_io/DatabasePurpose.hpp"   // for DatabasePurpose::READ_MESH, etc
#include "stk_topology/topology.hpp"    // for topology, etc


namespace Ioss { class DatabaseIO; }
namespace stk { namespace mesh { class FieldBase; } }
namespace {

void createNamedFieldOnMesh(stk::mesh::MetaData &stkMeshMetaData, const std::string &internalClientFieldName)
{
  const int numberOfStates = 1;
  stk::mesh::Field<double> &field0 = stkMeshMetaData.declare_field<double>(stk::topology::NODE_RANK, internalClientFieldName, numberOfStates);
  stk::mesh::put_field_on_mesh(field0, stkMeshMetaData.universal_part(), nullptr);
}

void testFieldNamedCorrectly(Ioss::Region &ioRegion, MPI_Comm communicator, std::vector<std::string> goldFieldNames)
{
  Ioss::NodeBlock *nodeBlockAssociatedWithField0 = ioRegion.get_node_blocks()[0];
  Ioss::NameList fieldNames;
  nodeBlockAssociatedWithField0->field_describe(Ioss::Field::TRANSIENT, &fieldNames);
  for (std::string & name : goldFieldNames) {
    for (char & c : name) {
      c = tolower(c);
    }
  }

  ASSERT_EQ(goldFieldNames.size(), fieldNames.size());
  for (auto goldFieldName : goldFieldNames) {
    auto entry = std::find(fieldNames.begin(), fieldNames.end(), goldFieldName.c_str());
    if (entry == fieldNames.end()) {
      EXPECT_TRUE(false) << "Field " << goldFieldName << " not found in file" << std::endl;
    }
  }
}

TEST(FieldNamesTest, FieldNameRenameTwice)
{
  const std::string outputFilename = "ourSillyOwlput.exo";
  MPI_Comm communicator = MPI_COMM_WORLD;
  const std::string internalClientFieldName = "Field0";
  std::string requestedFieldNameForResultsOutput("NotjeSSe");
  std::vector<std::string> outputFieldNames;
  {
    stk::io::StkMeshIoBroker stkIo(communicator);
    const std::string exodusFileName = "generated:1x1x8";
    size_t index = stkIo.add_mesh_database(exodusFileName, stk::io::READ_MESH);
    stkIo.set_active_mesh(index);
    stkIo.create_input_mesh();

    stk::mesh::MetaData &stkMeshMetaData = stkIo.meta_data();
    createNamedFieldOnMesh(stkMeshMetaData, internalClientFieldName);
    stkIo.populate_bulk_data();

    size_t results_output_index = stkIo.create_output_mesh(outputFilename, stk::io::WRITE_RESULTS);

    stk::mesh::FieldBase *field0 = stkMeshMetaData.get_field(stk::topology::NODE_RANK, internalClientFieldName);
    stkIo.add_field(results_output_index, *field0, requestedFieldNameForResultsOutput);
    outputFieldNames.push_back(requestedFieldNameForResultsOutput);


    requestedFieldNameForResultsOutput = "jeSSe";
    stkIo.add_field(results_output_index, *field0, requestedFieldNameForResultsOutput);
    outputFieldNames.push_back(requestedFieldNameForResultsOutput);

    stkIo.process_output_request(results_output_index, 0.0);
  }

  Ioss::DatabaseIO *iossDb = Ioss::IOFactory::create("exodus", outputFilename, Ioss::READ_MODEL, communicator);
  Ioss::Region ioRegion(iossDb);
  testFieldNamedCorrectly(ioRegion, communicator, outputFieldNames);

  unlink(outputFilename.c_str());
}

TEST(FieldNamesTest, FieldNameWithRestart)
{
  std::string restartFilename = "output.restart";
  MPI_Comm communicator = MPI_COMM_WORLD;
  const std::string internalClientFieldName = "Field0";
  {
    stk::io::StkMeshIoBroker stkIo(communicator);
    const std::string exodusFileName = "generated:1x1x8";
    size_t index = stkIo.add_mesh_database(exodusFileName, stk::io::READ_MESH);
    stkIo.set_active_mesh(index);
    stkIo.create_input_mesh();

    stk::mesh::MetaData &stkMeshMetaData = stkIo.meta_data();
    createNamedFieldOnMesh(stkMeshMetaData, internalClientFieldName);
    stkIo.populate_bulk_data();

    stk::mesh::FieldBase *field0 = stkMeshMetaData.get_field(stk::topology::NODE_RANK, internalClientFieldName);

    size_t fileIndex = stkIo.create_output_mesh(restartFilename, stk::io::WRITE_RESTART);
    stkIo.add_field(fileIndex, *field0);

    double time = 0.0;
    stkIo.begin_output_step(fileIndex, time);
    stkIo.write_defined_output_fields(fileIndex);
    stkIo.end_output_step(fileIndex);
  }

  Ioss::DatabaseIO *iossDb = Ioss::IOFactory::create("exodus", restartFilename, Ioss::READ_RESTART, communicator);
  Ioss::Region iossRegion(iossDb);
  std::string goldFieldName = sierra::make_lower(internalClientFieldName);
  testFieldNamedCorrectly(iossRegion, communicator, {goldFieldName});

  unlink(restartFilename.c_str());
}

TEST(FieldNamesTest, FieldNameWithResultsAndRestart)
{
  const std::string restartFilename = "FieldNameWithResultsAndRestart.restart";
  const std::string outputFileName = "FieldNameWithResultsAndRestart.exo";
  MPI_Comm communicator = MPI_COMM_WORLD;
  const std::string internalClientFieldName = "Field0";
  {
    stk::io::StkMeshIoBroker stkIo(communicator);
    const std::string exodusFileName = "generated:1x1x8";
    size_t index = stkIo.add_mesh_database(exodusFileName, stk::io::READ_MESH);
    stkIo.set_active_mesh(index);
    stkIo.create_input_mesh();

    stk::mesh::MetaData &stkMeshMetaData = stkIo.meta_data();
    createNamedFieldOnMesh(stkMeshMetaData, internalClientFieldName);
    stkIo.populate_bulk_data();

    size_t results_output_index = stkIo.create_output_mesh(outputFileName, stk::io::WRITE_RESULTS);
    stk::mesh::FieldBase *field0 = stkMeshMetaData.get_field(stk::topology::NODE_RANK, internalClientFieldName);
    std::string requestedFieldNameForResultsOutput("jeSSe");
    stkIo.add_field(results_output_index, *field0, requestedFieldNameForResultsOutput);

    size_t restartFileIndex = stkIo.create_output_mesh(restartFilename, stk::io::WRITE_RESTART);
    stkIo.add_field(restartFileIndex, *field0);

    double time = 0.0;
    stkIo.begin_output_step(results_output_index, time);
    stkIo.write_defined_output_fields(results_output_index);
    stkIo.end_output_step(results_output_index);

    stkIo.begin_output_step(restartFileIndex, time);
    stkIo.write_defined_output_fields(restartFileIndex);
    stkIo.end_output_step(restartFileIndex);
  }
  Ioss::DatabaseIO *iossResultDb = Ioss::IOFactory::create("exodus", restartFilename, Ioss::READ_MODEL, communicator);
  Ioss::Region resultRegion(iossResultDb);
  std::string goldFieldName = sierra::make_lower(internalClientFieldName);
  testFieldNamedCorrectly(resultRegion, communicator, {goldFieldName});

  Ioss::DatabaseIO *iossRestartDb = Ioss::IOFactory::create("exodus", restartFilename, Ioss::READ_RESTART, communicator);
  Ioss::Region restartRegion(iossRestartDb);
  testFieldNamedCorrectly(restartRegion, communicator, {goldFieldName});

  unlink(outputFileName.c_str());
  unlink(restartFilename.c_str());
}

}
