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

#include "gtest/gtest.h"
#include "stk_transfer_util/TransferMainIO.hpp"
#include "stk_unit_test_utils/ioUtils.hpp"
#include "stk_unit_test_utils/GeneratedMeshToFile.hpp"
#include "stk_unit_test_utils/CommandLineArgs.hpp"
#include <stk_util/parallel/OutputStreams.hpp>
#include "stk_io/WriteMesh.hpp"

namespace {

void build_serial_mesh_with_transient_data(const std::string& fileName, const std::string& fieldName,
                                           const std::vector<double>& timeSteps, const double fieldValueScaleFactor)
{
  stk::unit_test_util::generated_mesh_with_transient_data_to_file_in_serial("1x1x1", fileName,
                                                          fieldName, stk::topology::NODE_RANK,
                                                          "global_var", timeSteps,
                                                          stk::unit_test_util::IdAndTimeFieldValueSetter(),
                                                          fieldValueScaleFactor);
}

void perform_mock_transfer(std::shared_ptr<stk::mesh::BulkData> recvBulk, const double& lastTimeStep, double fieldScaleValue)
{
  stk::mesh::EntityRank fieldRank = stk::topology::NODE_RANK;

  std::vector<stk::mesh::Entity> entities;
  stk::mesh::get_entities(*recvBulk, fieldRank, entities);

  stk::mesh::FieldVector allTransientFields = stk::io::get_transient_fields(recvBulk->mesh_meta_data());

  for (stk::mesh::FieldBase* transientField : allTransientFields) {
    stk::mesh::field_data_execute<double, stk::mesh::ReadWrite>(*transientField,
      [&](auto& fieldData) {
        for (size_t i = 0; i < entities.size(); ++i) {
          double value = 100.0 * fieldScaleValue * lastTimeStep + static_cast<double>(recvBulk->identifier(entities[i]));
          auto data = fieldData.entity_values(entities[i]);
          for (stk::mesh::ComponentIdx j : data.components()) {
            data(j) = value + static_cast<int>(j);
          }
        }
      }
    );
  }
}

void verify_transferred_transient_field_data(stk::unit_test_util::MeshFromFile & mesh,
                                  stk::ParallelMachine& comm,
                                 const std::vector<double> & transientTimeSteps,
                                 const std::string & transientFieldName,
                                 const double fieldValueScaleFactor)
{
  stk::unit_test_util::TransientVerifier verifier(comm);
  verifier.verify_time_steps(mesh, transientTimeSteps);
  verifier.verify_num_transient_fields(mesh, 2);
  verifier.verify_transient_field_names(mesh, transientFieldName);
  verifier.verify_transient_fields(mesh, fieldValueScaleFactor);
}

TEST(TransferMainIO, mockTransfer)
{
  stk::ParallelMachine comm = stk::parallel_machine_world();
  if(stk::parallel_machine_size(comm) != 1) { GTEST_SKIP(); }

  std::string sendMesh = "sendInput.exo";
  std::string recvMesh = "receiveInput.exo";

  std::string fieldName = "field1";
  std::vector<double> timeSteps = {0.0, 1.0};
  double lastTimeStep = timeSteps.back();

  double sendFieldScaleValue = 1.0;
  double recvFieldScaleValue = 10.0;

  build_serial_mesh_with_transient_data(sendMesh, fieldName, timeSteps, sendFieldScaleValue);
  build_serial_mesh_with_transient_data(recvMesh, fieldName, timeSteps, recvFieldScaleValue);

  stk::transfer_util::TransferMainIO io(comm, sendMesh, recvMesh);
  io.load_meshes();

  std::shared_ptr<stk::mesh::BulkData> recvBulk = io.get_recvBulkData();
  perform_mock_transfer(recvBulk, lastTimeStep, sendFieldScaleValue);

  std::string transferOutput = "transferredReceive.exo";
  io.write_transfer_output({"field1_scalar", "field1_vector"}, transferOutput);

  stk::unit_test_util::MeshFromFile transferredFieldData(comm);
  transferredFieldData.fill_from_parallel(transferOutput);

  verify_transferred_transient_field_data(transferredFieldData, comm, {lastTimeStep}, fieldName, sendFieldScaleValue);

  unlink(sendMesh.c_str());
  unlink(recvMesh.c_str());
  unlink(transferOutput.c_str());
}

}
