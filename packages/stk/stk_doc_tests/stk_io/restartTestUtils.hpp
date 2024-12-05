#ifndef STK_RESTART_TEST_UTILS_H
#define STK_RESTART_TEST_UTILS_H

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
#include <vector>
#include <string>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/Field.hpp>
#include <stk_mesh/base/Types.hpp>
#include <stk_mesh/base/GetEntities.hpp>

namespace stk { namespace mesh { class BulkData; } }

inline void checkFileForNodalVarNames(const std::string &exodusFilename, const std::vector<std::string>& nodalVarNames)
{
  Ioss::DatabaseIO *iossDb = Ioss::IOFactory::create("exodus", exodusFilename, Ioss::READ_MODEL, MPI_COMM_WORLD);
  Ioss::Region ioRegion(iossDb);
  Ioss::NodeBlock *nodeBlockAssociatedWithField0 = ioRegion.get_node_blocks()[0];
  Ioss::NameList fieldNames;
  nodeBlockAssociatedWithField0->field_describe(Ioss::Field::TRANSIENT, &fieldNames);

  ASSERT_EQ(nodalVarNames.size(), fieldNames.size());
  for (size_t i=0;i<nodalVarNames.size();i++)
  {
    EXPECT_TRUE(nodeBlockAssociatedWithField0->field_exists(nodalVarNames[i])) << nodalVarNames[i];
  }
}

inline void checkFileForGlobal(const std::string &exodusFilename, const std::string &globalVarName, const double value)
{
  Ioss::DatabaseIO *iossDb = Ioss::IOFactory::create("exodus", exodusFilename, Ioss::READ_MODEL, MPI_COMM_WORLD);
  Ioss::Region ioRegion(iossDb);
  ioRegion.begin_state(1);
  ASSERT_TRUE(ioRegion.field_exists(globalVarName));
  double valueOnFile = 0.0;
  ioRegion.get_field_data(globalVarName, &valueOnFile, sizeof(double));
  EXPECT_EQ(value, valueOnFile);
}

inline stk::mesh::Field<double> & declareNodalField(stk::mesh::MetaData &stkMeshMetaData, const std::string &fieldName,
                                                    int numberOfStates)
{
  stk::mesh::Field<double> &multiStateField = stkMeshMetaData.declare_field<double>(stk::topology::NODE_RANK, fieldName, numberOfStates);
  stk::mesh::put_field_on_mesh(multiStateField, stkMeshMetaData.universal_part(), nullptr);
  return multiStateField;
}

inline stk::mesh::Field<double> & declareTriStateNodalField(stk::mesh::MetaData &stkMeshMetaData, const std::string &fieldName)
{
  const int numberOfStates = 3;
  return declareNodalField(stkMeshMetaData, fieldName, numberOfStates);
}

inline void putDataOnTestField(stk::mesh::BulkData &stkMeshBulkData, const double value, stk::mesh::FieldBase &field)
{
  std::vector<stk::mesh::Entity> nodes;
  stk::mesh::get_entities(stkMeshBulkData, stk::topology::NODE_RANK, nodes);
  for(size_t i=0; i<nodes.size(); i++)
  {
    double *fieldDataForNode = reinterpret_cast<double*>(stk::mesh::field_data(field, nodes[i]));
    *fieldDataForNode = value;
  }
}

inline void testDataOnField(stk::mesh::BulkData &stkMeshBulkData, const double goldValue, stk::mesh::FieldBase &field)
{
  std::vector<stk::mesh::Entity> nodes;
  stk::mesh::get_entities(stkMeshBulkData, stk::topology::NODE_RANK, nodes);
  for(size_t i=0; i<nodes.size(); i++)
  {
    double *fieldDataForNode = reinterpret_cast<double*>(stk::mesh::field_data(field, nodes[i]));
    EXPECT_DOUBLE_EQ(goldValue, *fieldDataForNode);
  }
}

inline void setupMeshAndFieldsForTest(stk::io::StkMeshIoBroker &stkMeshIoBroker,
                                      const std::string &displacementFieldName,
                                      const std::string &velocityFieldName)
{
  const double displacementValue = 1.0;
  const double velocityValue = 2.0;

  const std::string exodusFileName = "generated:1x1x8";
  size_t index = stkMeshIoBroker.add_mesh_database(exodusFileName, stk::io::READ_MESH);
  stkMeshIoBroker.set_active_mesh(index);
  stkMeshIoBroker.create_input_mesh();

  const int numberOfStates = 1;
  stk::mesh::Field<double> &displacementField = declareNodalField(stkMeshIoBroker.meta_data(), displacementFieldName, numberOfStates);
  stk::mesh::Field<double> &velocityField = declareNodalField(stkMeshIoBroker.meta_data(), velocityFieldName, numberOfStates);

  stkMeshIoBroker.populate_bulk_data();

  putDataOnTestField(stkMeshIoBroker.bulk_data(), displacementValue, displacementField);
  putDataOnTestField(stkMeshIoBroker.bulk_data(), velocityValue, velocityField);
}

//inline void testMultistateFieldWroteCorrectlyToRestart(const std::string &restartFilename,
//        const double time,
//        const std::string &fieldName,
//        const double stateNp1Value,
//        const double stateNValue)
//{
//    MPI_Comm communicator = MPI_COMM_WORLD;
//    stk::io::StkMeshIoBroker stkIo(communicator);
//    size_t index = stkIo.add_mesh_database(restartFilename, stk::io::READ_RESTART);
//    stkIo.set_active_mesh(index);
//    stkIo.create_input_mesh();
//
//    stk::mesh::MetaData &restartedMetaData = stkIo.meta_data();
//    stk::mesh::Field<double> &triStateField =
//            declareTriStateNodalField(restartedMetaData, fieldName);
//
//    stkIo.add_input_field(stk::io::MeshField(triStateField));
//    stkIo.populate_bulk_data();
//    stkIo.read_defined_input_fields(time);
//
//    stk::mesh::FieldBase *statedFieldNp1 =
//            triStateField.field_state(stk::mesh::StateNP1);
//    testDataOnField(stkIo.bulk_data(), stateNp1Value, *statedFieldNp1);
//    stk::mesh::FieldBase *statedFieldN =
//            triStateField.field_state(stk::mesh::StateN);
//    testDataOnField(stkIo.bulk_data(), stateNValue, *statedFieldN);
//}

inline void testMultistateFieldWroteCorrectly(const std::string &resultsFilename,
                                              const double time,
                                              const std::string &np1Name,
                                              const std::string &nName,
                                              const std::string &nm1Name,
                                              const double stateNp1Value,
                                              const double stateNValue,
                                              const double stateNm1Value)
{
  MPI_Comm communicator = MPI_COMM_WORLD;
  stk::io::StkMeshIoBroker stkIo(communicator);
  size_t index = stkIo.add_mesh_database(resultsFilename, stk::io::READ_RESTART);
  stkIo.set_active_mesh(index);
  stkIo.create_input_mesh();

  stk::mesh::MetaData &resultsedMetaData = stkIo.meta_data();
  stk::mesh::Field<double> &FieldNp1 = declareNodalField(resultsedMetaData, np1Name, 1);
  stk::mesh::Field<double> &FieldN   = declareNodalField(resultsedMetaData, nName, 1);
  stk::mesh::Field<double> &FieldNm1 = declareNodalField(resultsedMetaData, nm1Name, 1);

  stkIo.add_input_field(stk::io::MeshField(FieldNp1, np1Name));
  stkIo.add_input_field(stk::io::MeshField(FieldN,   nName));
  stkIo.add_input_field(stk::io::MeshField(FieldNm1, nm1Name));

  stkIo.populate_bulk_data();
  stkIo.read_defined_input_fields(time);

  testDataOnField(stkIo.bulk_data(), stateNp1Value, FieldNp1);
  testDataOnField(stkIo.bulk_data(), stateNValue,   FieldN);
  testDataOnField(stkIo.bulk_data(), stateNm1Value, FieldNm1);
}

#endif
