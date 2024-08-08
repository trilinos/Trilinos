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
#include <exception>                    // for exception
#include <stk_io/StkMeshIoBroker.hpp>   // for StkMeshIoBroker
#include <stk_mesh/base/BulkData.hpp>   // for BulkData
#include <stk_mesh/base/Field.hpp>      // for Field
#include <stk_mesh/base/GetEntities.hpp>  // for get_entities
#include <stk_mesh/base/MetaData.hpp>   // for MetaData, put_field
#include <gtest/gtest.h>
#include <string>                       // for string, basic_string
#include <vector>                       // for vector
#include "Ioss_DBUsage.h"               // for DatabaseUsage::READ_RESTART
#include "Ioss_Field.h"                 // for Field, etc
#include "Ioss_IOFactory.h"             // for IOFactory
#include "Ioss_NodeBlock.h"             // for NodeBlock
#include "Ioss_Region.h"                // for Region
#include "gtest/gtest.h"                // for AssertHelper, ASSERT_TRUE, etc
#include "stk_io/DatabasePurpose.hpp"   // for DatabasePurpose::READ_MESH, etc
#include "stk_mesh/base/Entity.hpp"     // for Entity
#include "stk_mesh/base/FieldBase.hpp"  // for field_data
#include "stk_topology/topology.hpp"    // for topology, etc
#include "stk_util/parallel/Parallel.hpp"
namespace Ioss { class DatabaseIO; }

namespace {

Ioss::Field::BasicType iossBasicType(double)
{
  return Ioss::Field::REAL;
}

Ioss::Field::BasicType iossBasicType(int)
{
  return Ioss::Field::INTEGER;
}

template <typename DataType>
void testGlobalVarOnFile(const std::string &outputFileName, const int stepNumber, const std::vector<std::string> &goldGlobalVarName,
                         const std::vector<DataType> goldGlobalVarValue, DataType goldGlobalScale, MPI_Comm comm)
{
  stk::io::StkMeshIoBroker stkIo(comm);
  stkIo.add_mesh_database(outputFileName, stk::io::READ_MESH);
  stkIo.create_input_mesh();
  stkIo.populate_bulk_data();
  stkIo.read_defined_input_fields(stepNumber);
  std::vector<std::string> globalVarNames;
  stkIo.get_global_variable_names(globalVarNames);
  ASSERT_EQ(goldGlobalVarName.size(), globalVarNames.size());
  for(size_t i=0; i<goldGlobalVarName.size(); i++)
  {
    EXPECT_STRCASEEQ(goldGlobalVarName[i].c_str(), globalVarNames[i].c_str());
    std::vector<DataType> globalVar;
    ASSERT_TRUE(stkIo.get_global(globalVarNames[i], globalVar));
    for(size_t j=0; j<globalVar.size(); j++)
    {
      EXPECT_NEAR(goldGlobalVarValue[i]+goldGlobalScale*j, globalVar[j], std::numeric_limits<double>::epsilon());
    }
  }
}

void testNodalFieldOnFile(const std::string &outputFileName, const int stepNumber, const std::string &goldNodalFieldName,
                          const std::vector<double> goldNodalFieldValues, MPI_Comm comm)
{
  Ioss::DatabaseIO *iossDb = Ioss::IOFactory::create("exodus", outputFileName, Ioss::READ_RESTART, comm);
  Ioss::Region inputRegion(iossDb);
  Ioss::NodeBlock *nodeBlock = inputRegion.get_node_block("nodeblock_1");
  ASSERT_TRUE(nodeBlock->field_exists(goldNodalFieldName));

  inputRegion.begin_state(stepNumber);

  std::vector<double> fieldValues;
  nodeBlock->get_field_data(goldNodalFieldName, fieldValues);
  ASSERT_EQ(goldNodalFieldValues.size(), fieldValues.size());
  for(size_t i=0; i<goldNodalFieldValues.size(); i++)
  {
    EXPECT_NEAR(goldNodalFieldValues[i], fieldValues[i], std::numeric_limits<double>::epsilon());
  }
}

TEST(GlobalVariablesTest, OneGlobalDouble)
{
  const std::string outputFileName = "OneGlobalDouble.exo";
  //                                 0        1         2         3         4         5
  //                                 12345678901234567890123456789012345678901234567890
  const std::string globalVarName = "a_global_variable_with_a_longer_name_than_the_default";
  const double globalVarValue = 13.0;
  MPI_Comm communicator = MPI_COMM_WORLD;
  {
    stk::io::StkMeshIoBroker stkIo(communicator);
    stkIo.property_add(Ioss::Property("MAXIMUM_NAME_LENGTH", 64));
    const std::string exodusFileName = "generated:1x1x8";
    stkIo.add_mesh_database(exodusFileName, stk::io::READ_MESH);
    stkIo.create_input_mesh();
    stkIo.populate_bulk_data();

    size_t result_file_index = stkIo.create_output_mesh(outputFileName, stk::io::WRITE_RESULTS);

    stkIo.add_global(result_file_index, globalVarName, Ioss::Field::REAL);

    const double time = 1.0;
    stkIo.begin_output_step(result_file_index, time);

    stkIo.write_global(result_file_index, globalVarName, globalVarValue);

    stkIo.end_output_step(result_file_index);
  }

  const int stepNumber = 1;
  std::vector<std::string> globalVarNames(1, globalVarName);
  std::vector<double> globalVarValues(1,globalVarValue);
  double goldGlobalScaleFactor = 0.0;
  testGlobalVarOnFile(outputFileName, stepNumber, globalVarNames, globalVarValues, goldGlobalScaleFactor, communicator);
  unlink(outputFileName.c_str());
}

TEST(GlobalVariablesTest, InvalidGlobalRequest)
{
  const std::string outputFileName = "InvalidGlobalRequest.exo";
  const std::string globalVarName = "testGlobal";
  const double globalVarValue = 13.0;
  MPI_Comm communicator = MPI_COMM_WORLD;
  {
    stk::io::StkMeshIoBroker stkIo(communicator);
    const std::string exodusFileName = "generated:1x1x8";
    stkIo.add_mesh_database(exodusFileName, stk::io::READ_MESH);
    stkIo.create_input_mesh();
    stkIo.populate_bulk_data();

    size_t result_file_index = stkIo.create_output_mesh(outputFileName, stk::io::WRITE_RESULTS);

    stkIo.add_global(result_file_index, globalVarName, Ioss::Field::REAL);

    const double time = 1.0;
    stkIo.begin_output_step(result_file_index, time);

    stkIo.write_global(result_file_index, globalVarName, globalVarValue);

    stkIo.end_output_step(result_file_index);
  }

  {
    double global_value = 0.0;
    stk::io::StkMeshIoBroker stkIo(communicator);
    stkIo.add_mesh_database(outputFileName, stk::io::READ_MESH);
    stkIo.create_input_mesh();
    stkIo.populate_bulk_data();

    bool abort_if_not_exist = true;
    EXPECT_ANY_THROW(stkIo.get_global("does_not_exist", global_value));
    EXPECT_ANY_THROW(stkIo.get_global("does_not_exist", global_value, abort_if_not_exist));
    ASSERT_TRUE(stkIo.get_global(globalVarName, global_value, abort_if_not_exist));

    abort_if_not_exist = false;
    ASSERT_FALSE(stkIo.get_global("does_not_exist", global_value, abort_if_not_exist));
    ASSERT_TRUE(stkIo.get_global(globalVarName, global_value, abort_if_not_exist));
  }
  unlink(outputFileName.c_str());
}

TEST(GlobalVariablesTest, OneGlobalDoubleVector3)
{
  const std::string outputFileName = "OneGlobalDoubleVector3.exo";
  const std::string globalVarName = "testGlobal";
  double goldGlobalScaleFactor = 0.1;
  std::vector<double> globalVarValue;
  for (int i=0; i < 3; i++) {
    globalVarValue.push_back(13.0 + i*goldGlobalScaleFactor);
  }

  MPI_Comm communicator = MPI_COMM_WORLD;
  {
    stk::io::StkMeshIoBroker stkIo(communicator);
    const std::string exodusFileName = "generated:1x1x8";
    stkIo.add_mesh_database(exodusFileName, stk::io::READ_MESH);
    stkIo.create_input_mesh();
    stkIo.populate_bulk_data();

    size_t result_file_index = stkIo.create_output_mesh(outputFileName, stk::io::WRITE_RESULTS);

    stkIo.add_global(result_file_index, globalVarName, "vector_3d", Ioss::Field::REAL);

    const double time = 1.0;
    stkIo.begin_output_step(result_file_index, time);

    stkIo.write_global(result_file_index, globalVarName, globalVarValue);

    stkIo.end_output_step(result_file_index);
  }

  const int stepNumber = 1;
  std::vector<std::string> globalVarNames(1, globalVarName);
  testGlobalVarOnFile(outputFileName, stepNumber, globalVarNames, globalVarValue, goldGlobalScaleFactor, communicator);
  unlink(outputFileName.c_str());
}

TEST(GlobalVariablesTest, OneGlobalIntegerVector3)
{
  const std::string outputFileName = "OneGlobalIntegerVector3.exo";
  const std::string globalVarName = "testGlobal";
  std::vector<int> globalVarValue;
  int goldGlobalScaleFactor = 10;
  for (int i=0; i < 3; i++) {
    globalVarValue.push_back(13 + i*goldGlobalScaleFactor);
  }

  MPI_Comm communicator = MPI_COMM_WORLD;
  {
    stk::io::StkMeshIoBroker stkIo(communicator);
    const std::string exodusFileName = "generated:1x1x8";
    stkIo.add_mesh_database(exodusFileName, stk::io::READ_MESH);
    stkIo.create_input_mesh();
    stkIo.populate_bulk_data();

    size_t result_file_index = stkIo.create_output_mesh(outputFileName, stk::io::WRITE_RESULTS);

    stkIo.add_global(result_file_index, globalVarName, "vector_3d", Ioss::Field::INTEGER);

    const double time = 1.0;
    stkIo.begin_output_step(result_file_index, time);

    stkIo.write_global(result_file_index, globalVarName, globalVarValue);

    stkIo.end_output_step(result_file_index);
  }

  const int stepNumber = 1;
  std::vector<std::string> globalVarNames(1, globalVarName);
  testGlobalVarOnFile(outputFileName, stepNumber, globalVarNames, globalVarValue, goldGlobalScaleFactor, communicator);
  unlink(outputFileName.c_str());
}

TEST(GlobalVariablesTest, OneGlobalDouble10)
{
  const std::string outputFileName = "OneGlobalDouble10.exo";
  const std::string globalVarName = "testGlobal";
  std::vector<double> globalVarValue;
  double goldGlobalScaleFactor = 0.1;
  for (int i=0; i < 10; i++) {
    globalVarValue.push_back(3.14159 + i*goldGlobalScaleFactor);
  }

  MPI_Comm communicator = MPI_COMM_WORLD;
  {
    stk::io::StkMeshIoBroker stkIo(communicator);
    const std::string exodusFileName = "generated:1x1x8";
    stkIo.add_mesh_database(exodusFileName, stk::io::READ_MESH);
    stkIo.create_input_mesh();
    stkIo.populate_bulk_data();

    size_t result_file_index = stkIo.create_output_mesh(outputFileName, stk::io::WRITE_RESULTS);

    stkIo.add_global(result_file_index, globalVarName, globalVarValue.size(), Ioss::Field::REAL);

    const double time = 1.0;
    stkIo.begin_output_step(result_file_index, time);

    stkIo.write_global(result_file_index, globalVarName, globalVarValue);

    stkIo.end_output_step(result_file_index);
  }

  const int stepNumber = 1;
  std::vector<std::string> globalVarNames(1, globalVarName);
  testGlobalVarOnFile(outputFileName, stepNumber, globalVarNames, globalVarValue, goldGlobalScaleFactor, communicator);
  unlink(outputFileName.c_str());
}

template <typename DataType>
void testTwoGlobals(const std::string &outputFileName, const std::vector<std::string> &globalVarNames)
{
  MPI_Comm communicator = MPI_COMM_WORLD;
  std::vector<DataType> globalVarValues = {13, 14};
  {
    stk::io::StkMeshIoBroker stkIo(communicator);
    const std::string exodusFileName = "generated:1x1x8";
    stkIo.add_mesh_database(exodusFileName, stk::io::READ_MESH);
    stkIo.create_input_mesh();
    stkIo.populate_bulk_data();

    size_t resultOuputIndex = stkIo.create_output_mesh(outputFileName, stk::io::WRITE_RESULTS);

    Ioss::Field::BasicType iossDataType = iossBasicType(DataType());
    stkIo.add_global(resultOuputIndex, globalVarNames[0], iossDataType);
    stkIo.add_global(resultOuputIndex, globalVarNames[1], iossDataType);

    const double time = 1.0;
    stkIo.begin_output_step(resultOuputIndex, time);

    stkIo.write_global(resultOuputIndex, globalVarNames[0], globalVarValues[0]);
    stkIo.write_global(resultOuputIndex, globalVarNames[1], globalVarValues[1]);

    stkIo.end_output_step(resultOuputIndex);
  }

  const int stepNumber = 1;
  DataType goldGlobalScaleFactor = 0;
  testGlobalVarOnFile(outputFileName, stepNumber, globalVarNames, globalVarValues, goldGlobalScaleFactor, communicator);
  unlink(outputFileName.c_str());
}

TEST(GlobalVariablesTest, TwoGlobalIntegers)
{
  std::vector<std::string> globalVarNames = {"testGlobal", "testGlobal2"};
  testTwoGlobals<int>("TwoGlobalIntegers.exo", globalVarNames);
}

TEST(GlobalVariablesTest, TwoGlobalDoubles)
{
  std::vector<std::string> globalVarNames = {"testGlobal", "testGlobal2"};
  testTwoGlobals<double>("TwoGlobalDoubles.exo", globalVarNames);
}

TEST(GlobalVariablesTest, TwoGlobalDoublesSameName)
{
  std::vector<std::string> globalVarNames = {"testGlobal", "testGlobal"};
  EXPECT_ANY_THROW(testTwoGlobals<double>("TwoGlobalDoublesSameName.exo", globalVarNames));
  unlink("TwoGlobalDoublesSameName.exo");
}

stk::mesh::Field<double> &createNodalTestField(stk::mesh::MetaData &stkMeshMetaData, const std::string &fieldName)
{
  const int numberOfStates = 1;
  stk::mesh::Field<double> &field0 = stkMeshMetaData.declare_field<double>(stk::topology::NODE_RANK, fieldName, numberOfStates);
  stk::mesh::put_field_on_mesh(field0, stkMeshMetaData.universal_part(), nullptr);
  return field0;
}
void putDataOnTestField(stk::mesh::BulkData &stkMeshBulkData, stk::mesh::Field<double> &field0, std::vector<double> &nodalFieldValues)
{
  std::vector<stk::mesh::Entity> nodes;
  stk::mesh::get_entities(stkMeshBulkData, stk::topology::NODE_RANK, nodes);
  for(size_t i=0; i<nodes.size(); i++)
  {
    double *fieldDataForNode = stk::mesh::field_data(field0, nodes[i]);
    *fieldDataForNode = static_cast<double>(stkMeshBulkData.identifier(nodes[i]));
    nodalFieldValues.push_back(*fieldDataForNode);
  }
}

TEST(GlobalVariablesTest, GlobalDoubleWithFieldMultipleTimeSteps)
{
  const std::string outputFileName = "GlobalDoubleWithFieldMultipleTimeSteps.exo";
  const std::string fieldName = "field0";
  std::vector<double> nodalFieldValues;
  const std::string globalVarName = "testGlobal";
  std::vector<double> globalVarValuesOverTime;
  const int numTimeSteps = 5;
  MPI_Comm communicator = MPI_COMM_WORLD;
  int numProcs = stk::parallel_machine_size(communicator);
  if (numProcs != 1) {
    return;
  }
  {
    stk::io::StkMeshIoBroker stkIo(communicator);
    const std::string exodusFileName = "generated:1x1x8";
    stkIo.add_mesh_database(exodusFileName, stk::io::READ_MESH);
    stkIo.create_input_mesh();

    stk::mesh::Field<double> &field0 = createNodalTestField(stkIo.meta_data(), fieldName);

    stkIo.populate_bulk_data();

    putDataOnTestField(stkIo.bulk_data(), field0, nodalFieldValues);

    size_t result_file_index = stkIo.create_output_mesh(outputFileName, stk::io::WRITE_RESULTS);
    stkIo.add_field(result_file_index, field0);
    stkIo.add_global(result_file_index, globalVarName, Ioss::Field::REAL);

    double time = 1.0;
    const double stepSize = 1.0;
    for(int i=0; i<numTimeSteps; i++)
    {
      stkIo.begin_output_step(result_file_index, time);

      const double globalVarValue = time;
      stkIo.write_global(result_file_index, globalVarName, globalVarValue);
      globalVarValuesOverTime.push_back(globalVarValue);

      stkIo.write_defined_output_fields(result_file_index);

      stkIo.end_output_step(result_file_index);
      time += stepSize;
    }
  }

  for(int i=0; i<numTimeSteps; i++)
  {
    std::vector<std::string> globalVarNames(1, globalVarName);
    std::vector<double> globalVarValues(1,globalVarValuesOverTime[i]);
    double goldGlobalScaleFactor = 0.0;
    testGlobalVarOnFile(outputFileName, i+1, globalVarNames, globalVarValues, goldGlobalScaleFactor, communicator);
    testNodalFieldOnFile(outputFileName, i+1, fieldName, nodalFieldValues, communicator);
  }
  unlink(outputFileName.c_str());
}

TEST(GlobalVariablesTest, OneGlobalDoubleRestart)
{
  const std::string restartFileName = "OneGlobalDouble.restart";
  const std::string globalVarName = "testGlobal";
  const double globalVarValue = 13.0;
  const double time = 1.0;
  MPI_Comm communicator = MPI_COMM_WORLD;
  {
    stk::io::StkMeshIoBroker stkIo(communicator);
    const std::string exodusFileName = "generated:1x1x8";
    stkIo.add_mesh_database(exodusFileName, stk::io::READ_MESH);
    stkIo.create_input_mesh();
    stkIo.populate_bulk_data();

    size_t fileIndex = stkIo.create_output_mesh(restartFileName, stk::io::WRITE_RESTART);

    stkIo.add_global(fileIndex, globalVarName, Ioss::Field::REAL);

    stkIo.begin_output_step(fileIndex, time);

    stkIo.write_global(fileIndex, globalVarName, globalVarValue);

    stkIo.end_output_step(fileIndex);
  }

  {
    stk::io::StkMeshIoBroker stkIo(communicator);
    stkIo.add_mesh_database(restartFileName, stk::io::READ_RESTART);
    stkIo.create_input_mesh();
    stkIo.populate_bulk_data();
    stkIo.read_defined_input_fields(time);
    std::vector<std::string> globalVarNames;
    stkIo.get_global_variable_names(globalVarNames);
    ASSERT_EQ(1u, globalVarNames.size());
    EXPECT_STRCASEEQ(globalVarName.c_str(), globalVarNames[0].c_str());
    double globalVar = 0.0;
    ASSERT_TRUE(stkIo.get_global(globalVarNames[0], globalVar));
    EXPECT_NEAR(globalVarValue, globalVar, std::numeric_limits<double>::epsilon());
  }
  const int stepNumber = 1;
  std::vector<std::string> globalVarNames(1, globalVarName);
  std::vector<double> globalVarValues(1,globalVarValue);
  double goldGlobalScaleFactor = 0.0;
  testGlobalVarOnFile(restartFileName, stepNumber, globalVarNames, globalVarValues, goldGlobalScaleFactor,communicator);
  unlink(restartFileName.c_str());
}

TEST(GlobalVariablesTest, OneGlobalDoubleWithFieldRestart)
{
  const std::string outputFileName = "GlobalDoubleWithFieldMultipleTimeSteps.restart";
  const std::string fieldName = "field0";
  std::vector<double> nodalFieldValues;
  const std::string globalVarName = "testGlobal";
  std::vector<double> globalVarValuesOverTime;
  const int numTimeSteps = 5;
  MPI_Comm communicator = MPI_COMM_WORLD;
  int numProcs = stk::parallel_machine_size(communicator);
  if (numProcs != 1) {
    return;
  }

  {
    stk::io::StkMeshIoBroker stkIo(communicator);
    const std::string exodusFileName = "generated:1x1x8";
    stkIo.add_mesh_database(exodusFileName, stk::io::READ_MESH);
    stkIo.create_input_mesh();

    stk::mesh::Field<double> &field0 = createNodalTestField(stkIo.meta_data(), fieldName);

    stkIo.populate_bulk_data();

    putDataOnTestField(stkIo.bulk_data(), field0, nodalFieldValues);

    size_t fileIndex = stkIo.create_output_mesh(outputFileName, stk::io::WRITE_RESTART);
    stkIo.add_field(fileIndex, field0);

    stkIo.add_global(fileIndex, globalVarName, Ioss::Field::REAL);

    double time = 1.0;
    const double stepSize = 1.0;
    for(int i=0; i<numTimeSteps; i++)
    {
      stkIo.begin_output_step(fileIndex, time);

      const double globalVarValue = time;
      stkIo.write_global(fileIndex, globalVarName, globalVarValue);
      globalVarValuesOverTime.push_back(globalVarValue);

      stkIo.write_defined_output_fields(fileIndex);

      stkIo.end_output_step(fileIndex);
      time += stepSize;
    }
  }

  for(int i=0; i<numTimeSteps; i++)
  {
    std::vector<std::string> globalVarNames(1, globalVarName);
    std::vector<double> globalVarValues(1,globalVarValuesOverTime[i]);
    double goldGlobalScaleFactor = 0.0;
    testGlobalVarOnFile(outputFileName, i+1, globalVarNames, globalVarValues, goldGlobalScaleFactor, communicator);
    testNodalFieldOnFile(outputFileName, i+1, fieldName, nodalFieldValues, communicator);
  }
  unlink(outputFileName.c_str());
}

}
