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
#include <stk_io/IossBridge.hpp>        // for is_part_io_part
#include <stk_io/StkMeshIoBroker.hpp>   // for StkMeshIoBroker
#include <stk_mesh/base/BulkData.hpp>   // for BulkData
#include "stk_unit_test_utils/GeneratedMeshToFile.hpp"
#include "stk_unit_test_utils/ioUtils.hpp"
#include <stddef.h>                     // for size_t
#include <stdlib.h>                     // for rand, srand, RAND_MAX

namespace {

void write_mesh_with_transient_field_data(const std::string & fileName,
                                          const std::vector<double> & transientTimeSteps,
                                          const std::string & transientFieldName)
{
  std::string globalVariableName = "global_variable";
  stk::unit_test_util::GeneratedMeshToFileWithTransientFields gMesh(MPI_COMM_WORLD,
                                                                                   stk::mesh::BulkData::AUTO_AURA,
                                                                                   transientFieldName,
                                                                                   stk::topology::NODE_RANK);
  gMesh.setup_mesh("1x1x8", fileName);
  gMesh.write_mesh_with_field(transientTimeSteps,
                              stk::unit_test_util::IdAndTimeFieldValueSetter(),
                              globalVariableName);
}

void verify_transient_field_data(stk::unit_test_util::MeshFromFile & mesh,
                                 const std::vector<double> & transientTimeSteps,
                                 const std::string & transientFieldName)
{
  stk::unit_test_util::TransientVerifier verifier(MPI_COMM_WORLD);
  verifier.verify_time_steps(mesh, transientTimeSteps);
  verifier.verify_num_transient_fields(mesh, 2);
  verifier.verify_transient_field_names(mesh, transientFieldName);
  verifier.verify_transient_fields(mesh);
}


TEST(StkMeshIoBroker, readTransientFieldData) {
  const std::string fieldDataFile = "meshWithFieldData.e";
  std::vector<double> transientTimeSteps = {0.0, 1.0, 2.0};
  std::string transientFieldName = "transient_field";

  write_mesh_with_transient_field_data(fieldDataFile, transientTimeSteps, transientFieldName);

  stk::unit_test_util::MeshFromFile meshWithFieldData(MPI_COMM_WORLD);
  meshWithFieldData.fill_from_parallel(fieldDataFile);

  verify_transient_field_data(meshWithFieldData, transientTimeSteps, transientFieldName);
}

TEST(StkMeshIoBroker, readTransientFieldData_withCache) {
  const std::string fieldDataFile = "meshWithFieldData.e";
  std::vector<double> transientTimeSteps = {0.0, 1.0, 2.0};
  std::string transientFieldName = "transient_field";

  write_mesh_with_transient_field_data(fieldDataFile, transientTimeSteps, transientFieldName);

  stk::unit_test_util::MeshFromFile meshWithFieldData(MPI_COMM_WORLD);
  meshWithFieldData.broker.cache_entity_list_for_transient_steps(true);
  meshWithFieldData.fill_from_parallel(fieldDataFile);

  verify_transient_field_data(meshWithFieldData, transientTimeSteps, transientFieldName);
}

TEST(StkMeshIoBroker, missingInputField) {
  const std::string fieldDataFile = "meshWithMissingFieldData.e";
  std::vector<double> transientTimeSteps = {0.0, 1.0, 2.0};
  std::string transientFieldName = "transient_field";

  write_mesh_with_transient_field_data(fieldDataFile, transientTimeSteps, transientFieldName);

  std::unique_ptr<stk::mesh::BulkData> bulk = stk::mesh::MeshBuilder(MPI_COMM_WORLD).set_spatial_dimension(3).create();
  stk::mesh::MetaData& meta = bulk->mesh_meta_data();

  const stk::mesh::EntityRank rank = stk::topology::NODE_RANK;
  const std::string   fieldName =  transientFieldName+"_scalar";
  const std::string dbFieldName =  fieldName + "_missingField";

  stk::mesh::Field<double> &scalarField = meta.declare_field<double>(rank, fieldName, 1);
  stk::mesh::put_field_on_mesh(scalarField, meta.universal_part(), nullptr);

  stk::io::MeshField meshField(&scalarField, dbFieldName);
  stk::io::StkMeshIoBroker broker(MPI_COMM_WORLD);

  broker.set_throw_on_missing_input_fields(true);
  broker.set_bulk_data(*bulk);
  broker.add_mesh_database(fieldDataFile, stk::io::READ_MESH);
  broker.create_input_mesh();
  broker.add_input_field(meshField);
  broker.populate_bulk_data();

  std::vector<stk::io::MeshField> missingFields;
  EXPECT_THROW(broker.read_defined_input_fields(0.0, &missingFields), std::logic_error);

  unlink(fieldDataFile.c_str());
}

}
