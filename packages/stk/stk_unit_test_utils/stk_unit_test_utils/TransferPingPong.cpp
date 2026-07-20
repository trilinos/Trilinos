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

#include <stk_unit_test_utils/TransferPingPong.hpp>
#include <stk_unit_test_utils/meshCreationHelpers.hpp>
#include <stk_io/WriteMesh.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/FieldBLAS.hpp>

#include <utility>

namespace stk {
namespace unit_test_util {

TransferPingPong::TransferPingPong(MPI_Comm comm,
                   const std::string& sendMeshName,
                   const std::string& sendFieldName,
                   const FieldEvaluator& sendFieldEval,
                   const std::string& recvMeshName,
                   const std::string& recvFieldName,
                   const FieldEvaluator& recvFieldEval,
                   std::shared_ptr<stk::search::MasterElementProviderInterface> masterElemProvider)
 : m_pingSettings(),
   m_pongSettings(),
   m_pingBroker(),
   m_pongBroker(),
   m_sendFieldEval(sendFieldEval),
   m_errFieldName("error")
{
  std::shared_ptr<stk::mesh::BulkData> sendMesh =
      create_mesh_with_field(comm, sendFieldName, sendFieldEval, stk::topology::NODE_RANK,
                             Ioss::Field::TRANSIENT, sendMeshName);
  std::shared_ptr<stk::mesh::BulkData> recvMesh =
      create_mesh_with_field(comm, recvFieldName, recvFieldEval, stk::topology::NODE_RANK,
                             Ioss::Field::TRANSIENT, recvMeshName);

  ConstantFieldEvaluator fieldEvalConst(999.99);
  add_field_to_mesh(*sendMesh, m_errFieldName, fieldEvalConst);

  m_pingSettings.set_transfer_field(std::make_pair(sendFieldName,recvFieldName));
  m_pongSettings.set_transfer_field(std::make_pair(recvFieldName,sendFieldName));

  m_pingBroker = std::make_shared<transfer_util::TransferMainBroker>(comm, sendMesh, recvMesh, m_pingSettings, masterElemProvider);
  m_pongBroker = std::make_shared<transfer_util::TransferMainBroker>(comm, recvMesh, sendMesh, m_pongSettings, masterElemProvider);

  m_pingBroker->check_and_create_fields();
  m_pongBroker->check_and_create_fields();
}
                   
void TransferPingPong::run_steps(int numSteps)
{
  const stk::mesh::BulkData& sendBulk = *m_pingBroker->get_send_bulk();
  const stk::mesh::MetaData& sendMeta = sendBulk.mesh_meta_data();
  auto sendFieldNames = m_pingBroker->get_send_fields_for_transfer();
  STK_ThrowRequireMsg(sendFieldNames.size() == 1, "Ping-pong-demo only supports 1 send-field.");
  const stk::mesh::FieldBase* sendField = sendMeta.get_field(stk::topology::NODE_RANK, sendFieldNames[0]);
  STK_ThrowRequireMsg(sendField != nullptr,"failed to retrieve sendField.");
  const stk::mesh::FieldBase* errField = sendMeta.get_field(stk::topology::NODE_RANK, m_errFieldName);
  STK_ThrowRequireMsg(errField != nullptr,"failed to retrieve errField.");

  for(int step=0; step<numSteps; ++step) {
    m_pingBroker->transfer();
    m_pongBroker->transfer();
    set_error_field(sendBulk, *sendField, *errField, m_sendFieldEval);
  }
}

double TransferPingPong::get_max_err() const
{
  double maxErr = 0.0;
  const stk::mesh::BulkData& sendBulk = *m_pingBroker->get_send_bulk();
  const stk::mesh::MetaData& sendMeta = sendBulk.mesh_meta_data();
  const stk::mesh::FieldBase* errField = sendMeta.get_field(stk::topology::NODE_RANK, m_errFieldName);
  STK_ThrowRequireMsg(errField != nullptr,"failed to retrieve errField.");

  stk::mesh::field_amax<stk::ngp::HostSpace>(maxErr, *errField);

  return maxErr;
}

double run_ping_pong_transfer(MPI_Comm comm,
                              int numSteps,
                              const FieldEvaluator& sendFieldEval,
                              std::shared_ptr<stk::search::MasterElementProviderInterface> masterElemProvider,
                              const std::string& sendMeshName,
                              const std::string& recvMeshName)
{
  ConstantFieldEvaluator recvFieldEval(999.99);
  TransferPingPong xfer(comm, sendMeshName, "field_1", sendFieldEval,
                              recvMeshName, "field_1", recvFieldEval,
                              masterElemProvider);
  xfer.run_steps(numSteps);
  return xfer.get_max_err();
}

}}

