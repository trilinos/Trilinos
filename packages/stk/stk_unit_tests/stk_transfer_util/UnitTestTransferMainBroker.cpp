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
#include "stk_transfer_util/TransferMainBroker.hpp"
#include "stk_unit_test_utils/ioUtils.hpp"
#include "stk_unit_test_utils/GeneratedMeshToFile.hpp"
#include "stk_unit_test_utils/CommandLineArgs.hpp"
#include "stk_unit_test_utils/BuildMesh.hpp"
#include "stk_mesh/base/ForEachEntity.hpp"

namespace {

std::shared_ptr<stk::mesh::BulkData>
create_mesh_with_field(const stk::ParallelMachine comm,
                       const std::string& fieldName,
                       const double fieldValue,
                       const stk::mesh::EntityRank fieldRank = stk::topology::NODE_RANK,
                       const Ioss::Field::RoleType fieldRole = Ioss::Field::TRANSIENT,
                       const std::string generatedMeshString = "generated:4x1x1")
{
  std::shared_ptr<stk::mesh::BulkData> bulk = stk::unit_test_util::build_mesh(3u, comm);
  stk::mesh::Field<double> &field=
          bulk->mesh_meta_data().declare_field<double>(fieldRank, fieldName, 1);
  stk::io::set_field_role(field, fieldRole);
  stk::mesh::put_field_on_mesh(field, bulk->mesh_meta_data().universal_part(), &fieldValue);
  stk::io::fill_mesh(generatedMeshString, *bulk);

  return bulk;
}

void check_field_values(const std::shared_ptr<stk::mesh::BulkData> bulk, 
                        const std::string& fieldName,
                        const double expectedValue,
                        const stk::mesh::EntityRank fieldRank = stk::topology::NODE_RANK)
{
  auto field = bulk->mesh_meta_data().get_field(fieldRank, fieldName);
  auto fieldData = field->data<double, stk::mesh::ReadOnly>();
  
  auto checkFieldValues = [&](const stk::mesh::BulkData& , stk::mesh::Entity node) {
    auto dataForNode = fieldData.entity_values(node);
    for(stk::mesh::ComponentIdx i=0_comp; i<dataForNode.num_components(); ++i) {
      EXPECT_EQ(expectedValue, dataForNode(i));
    }
  };

  stk::mesh::for_each_entity_run(*bulk, fieldRank, checkFieldValues);
  
}

void check_recv_tranfer_fields(const std::vector<std::string>& recvFieldNames,
                               const std::vector<std::string>& expectedRecvFieldNames)
{
  EXPECT_EQ(expectedRecvFieldNames.size(), recvFieldNames.size());
  EXPECT_EQ(expectedRecvFieldNames, recvFieldNames);
}

TEST(TransferMainBroker, basicTransferSameSendRecvFieldName)
{
  stk::ParallelMachine comm = stk::parallel_machine_world();
  if(stk::parallel_machine_size(comm) != 1) { GTEST_SKIP(); }

  std::string sendFieldName = "field_1";
  double sendFieldValue = 1.234;
  auto sendBulk = create_mesh_with_field(comm, sendFieldName, sendFieldValue);

  std::string recvFieldName = "field_1";
  double recvFieldValue = 4.321;
  auto recvBulk = create_mesh_with_field(comm, recvFieldName, recvFieldValue);

  stk::transfer_util::TransferMainSettings settings;
  std::vector<std::pair<std::string, std::string>> fieldList = {{sendFieldName, recvFieldName}};
  settings.set_transfer_field(fieldList[0]);

  stk::transfer_util::TransferMainBroker broker(comm, sendBulk, recvBulk, settings);

  broker.check_fields();
  broker.check_parts();

  check_field_values(sendBulk, sendFieldName, sendFieldValue);
  check_field_values(recvBulk, recvFieldName, recvFieldValue);

  std::vector<std::string> recvFieldsForTransfer = broker.get_recv_fields_for_transfer();
  std::vector<std::string> expectedFields = {recvFieldName};
  check_recv_tranfer_fields(recvFieldsForTransfer, expectedFields);

  broker.transfer();

  check_field_values(sendBulk, sendFieldName, sendFieldValue);
  check_field_values(recvBulk, recvFieldName, sendFieldValue);

}

TEST(TransferMainBroker, basicTransferDifferentSendRecvFieldName)
{
  stk::ParallelMachine comm = stk::parallel_machine_world();
  if(stk::parallel_machine_size(comm) != 1) { GTEST_SKIP(); }

  std::string sendFieldName = "send_field_1";
  double sendFieldValue = 1.234;
  auto sendBulk = create_mesh_with_field(comm, sendFieldName, sendFieldValue);

  std::string recvFieldName = "recv_field_1";
  double recvFieldValue = 4.321;
  auto recvBulk = create_mesh_with_field(comm, recvFieldName, recvFieldValue);

  stk::transfer_util::TransferMainSettings settings;
  std::vector<std::pair<std::string, std::string>> fieldList = {{sendFieldName, recvFieldName}};
  settings.set_transfer_field(fieldList[0]);

  stk::transfer_util::TransferMainBroker broker(comm, sendBulk, recvBulk, settings);

  broker.check_fields();
  broker.check_parts();

  check_field_values(sendBulk, sendFieldName, sendFieldValue);
  check_field_values(recvBulk, recvFieldName, recvFieldValue);

  std::vector<std::string> recvFieldsForTransfer = broker.get_recv_fields_for_transfer();
  std::vector<std::string> expectedFields = {recvFieldName};
  check_recv_tranfer_fields(recvFieldsForTransfer, expectedFields);

  broker.transfer();

  check_field_values(sendBulk, sendFieldName, sendFieldValue);
  check_field_values(recvBulk, recvFieldName, sendFieldValue);

}

TEST(TransferMainBroker, basicTransferSendFieldOnly)
{
  stk::ParallelMachine comm = stk::parallel_machine_world();
  if(stk::parallel_machine_size(comm) != 1) { GTEST_SKIP(); }

  std::string sendFieldName = "send_field_1";
  double sendFieldValue = 1.234;
  auto sendBulk = create_mesh_with_field(comm, sendFieldName, sendFieldValue);

  std::shared_ptr<stk::mesh::BulkData> recvBulk = stk::unit_test_util::build_mesh(3u, comm);
  stk::io::fill_mesh("generated:4x1x1", *recvBulk);

  stk::transfer_util::TransferMainSettings settings;
  std::vector<std::pair<std::string, std::string>> fieldList = {{sendFieldName, sendFieldName}};
  settings.set_transfer_field(fieldList[0]);

  stk::transfer_util::TransferMainBroker broker(comm, sendBulk, recvBulk, settings);

  broker.check_fields();
  broker.check_parts();

  check_field_values(sendBulk, sendFieldName, sendFieldValue);

  std::vector<std::string> recvFieldsForTransfer = broker.get_recv_fields_for_transfer();
  std::vector<std::string> expectedFields = {sendFieldName};
  check_recv_tranfer_fields(recvFieldsForTransfer, expectedFields);

  broker.transfer();

  check_field_values(sendBulk, sendFieldName, sendFieldValue);
  check_field_values(recvBulk, sendFieldName, sendFieldValue);

}

TEST(TransferMainBroker, basicTransferNoFieldListSendRecvFieldName)
{
  stk::ParallelMachine comm = stk::parallel_machine_world();
  if(stk::parallel_machine_size(comm) != 1) { GTEST_SKIP(); }

  std::string sendFieldName = "field_1";
  double sendFieldValue = 1.234;
  auto sendBulk = create_mesh_with_field(comm, sendFieldName, sendFieldValue);

  std::string recvFieldName = "field_1";
  double recvFieldValue = 4.321;
  auto recvBulk = create_mesh_with_field(comm, recvFieldName, recvFieldValue);

  stk::transfer_util::TransferMainSettings settings;

  stk::transfer_util::TransferMainBroker broker(comm, sendBulk, recvBulk, settings);

  broker.check_fields();
  broker.check_parts();

  check_field_values(sendBulk, sendFieldName, sendFieldValue);
  check_field_values(recvBulk, recvFieldName, recvFieldValue);

  std::vector<std::string> recvFieldsForTransfer = broker.get_recv_fields_for_transfer();
  std::vector<std::string> expectedFields = {recvFieldName};
  check_recv_tranfer_fields(recvFieldsForTransfer, expectedFields);

  broker.transfer();

  check_field_values(sendBulk, sendFieldName, sendFieldValue);
  check_field_values(recvBulk, recvFieldName, sendFieldValue);

}

TEST(TransferMainBroker, basicTransferNoFieldListSendFieldOnly)
{
  stk::ParallelMachine comm = stk::parallel_machine_world();
  if(stk::parallel_machine_size(comm) != 1) { GTEST_SKIP(); }

  std::string sendFieldName = "field_1";
  double sendFieldValue = 1.234;
  auto sendBulk = create_mesh_with_field(comm, sendFieldName, sendFieldValue);

  std::shared_ptr<stk::mesh::BulkData> recvBulk = stk::unit_test_util::build_mesh(3u, comm);
  stk::io::fill_mesh("generated:4x1x1", *recvBulk);

  stk::transfer_util::TransferMainSettings settings;

  stk::transfer_util::TransferMainBroker broker(comm, sendBulk, recvBulk, settings);

  broker.check_fields();
  broker.check_parts();

  check_field_values(sendBulk, sendFieldName, sendFieldValue);

  std::vector<std::string> recvFieldsForTransfer = broker.get_recv_fields_for_transfer();
  std::vector<std::string> expectedFields = {sendFieldName};
  check_recv_tranfer_fields(recvFieldsForTransfer, expectedFields);

  broker.transfer();

  check_field_values(sendBulk, sendFieldName, sendFieldValue);
  check_field_values(recvBulk, sendFieldName, sendFieldValue);

}

TEST(TransferMainBroker, sendFieldWrongRank)
{
  stk::ParallelMachine comm = stk::parallel_machine_world();
  if(stk::parallel_machine_size(comm) != 1) { GTEST_SKIP(); }

  std::string sendFieldName = "field_1";
  double sendFieldValue = 1.234;
  auto sendBulk = create_mesh_with_field(comm, sendFieldName, sendFieldValue, stk::topology::ELEM_RANK);

  std::string recvFieldName = "field_1";
  double recvFieldValue = 4.321;
  auto recvBulk = create_mesh_with_field(comm, recvFieldName, recvFieldValue);

  stk::transfer_util::TransferMainSettings settings;
  std::vector<std::pair<std::string, std::string>> fieldList = {{sendFieldName, recvFieldName}};
  settings.set_transfer_field(fieldList[0]);

  stk::transfer_util::TransferMainBroker broker(comm, sendBulk, recvBulk, settings);

  EXPECT_ANY_THROW(broker.check_fields());

}

TEST(TransferMainBroker, sendFieldWrongRankNoFieldList)
{
  stk::ParallelMachine comm = stk::parallel_machine_world();
  if(stk::parallel_machine_size(comm) != 1) { GTEST_SKIP(); }

  std::string sendFieldName = "field_1";
  double sendFieldValue = 1.234;
  auto sendBulk = create_mesh_with_field(comm, sendFieldName, sendFieldValue, stk::topology::ELEM_RANK);

  std::string recvFieldName = "field_1";
  double recvFieldValue = 4.321;
  auto recvBulk = create_mesh_with_field(comm, recvFieldName, recvFieldValue);

  stk::transfer_util::TransferMainSettings settings;

  stk::transfer_util::TransferMainBroker broker(comm, sendBulk, recvBulk, settings);

  EXPECT_ANY_THROW(broker.check_fields());

}

TEST(TransferMainBroker, recvFieldWrongRank)
{
  stk::ParallelMachine comm = stk::parallel_machine_world();
  if(stk::parallel_machine_size(comm) != 1) { GTEST_SKIP(); }

  std::string sendFieldName = "field_1";
  double sendFieldValue = 1.234;
  auto sendBulk = create_mesh_with_field(comm, sendFieldName, sendFieldValue);

  std::string recvFieldName = "field_1";
  double recvFieldValue = 4.321;
  auto recvBulk = create_mesh_with_field(comm, recvFieldName, recvFieldValue, stk::topology::ELEM_RANK);

  stk::transfer_util::TransferMainSettings settings;
  std::vector<std::pair<std::string, std::string>> fieldList = {{sendFieldName, recvFieldName}};
  settings.set_transfer_field(fieldList[0]);

  stk::transfer_util::TransferMainBroker broker(comm, sendBulk, recvBulk, settings);

  broker.check_fields();
  broker.check_parts();

  check_field_values(sendBulk, sendFieldName, sendFieldValue);

  std::vector<std::string> recvFieldsForTransfer = broker.get_recv_fields_for_transfer();
  std::vector<std::string> expectedFields = {recvFieldName};
  check_recv_tranfer_fields(recvFieldsForTransfer, expectedFields);

  broker.transfer();

  check_field_values(sendBulk, sendFieldName, sendFieldValue);
  check_field_values(recvBulk, recvFieldName, sendFieldValue);
  check_field_values(recvBulk, recvFieldName, recvFieldValue, stk::topology::ELEM_RANK);

}

TEST(TransferMainBroker, recvFieldWrongRankNoFieldList)
{
  stk::ParallelMachine comm = stk::parallel_machine_world();
  if(stk::parallel_machine_size(comm) != 1) { GTEST_SKIP(); }

  std::string sendFieldName = "field_1";
  double sendFieldValue = 1.234;
  auto sendBulk = create_mesh_with_field(comm, sendFieldName, sendFieldValue);

  std::string recvFieldName = "field_1";
  double recvFieldValue = 4.321;
  auto recvBulk = create_mesh_with_field(comm, recvFieldName, recvFieldValue, stk::topology::ELEM_RANK);

  stk::transfer_util::TransferMainSettings settings;

  stk::transfer_util::TransferMainBroker broker(comm, sendBulk, recvBulk, settings);

  broker.check_fields();
  broker.check_parts();

  check_field_values(sendBulk, sendFieldName, sendFieldValue);

  std::vector<std::string> recvFieldsForTransfer = broker.get_recv_fields_for_transfer();
  std::vector<std::string> expectedFields = {recvFieldName};
  check_recv_tranfer_fields(recvFieldsForTransfer, expectedFields);

  broker.transfer();

  check_field_values(sendBulk, sendFieldName, sendFieldValue);
  check_field_values(recvBulk, recvFieldName, sendFieldValue);
  check_field_values(recvBulk, recvFieldName, recvFieldValue, stk::topology::ELEM_RANK);

}

TEST(TransferMainBroker, recvFieldWrongRankDifferentSendRecv)
{
  stk::ParallelMachine comm = stk::parallel_machine_world();
  if(stk::parallel_machine_size(comm) != 1) { GTEST_SKIP(); }

  std::string sendFieldName = "send_field_1";
  double sendFieldValue = 1.234;
  auto sendBulk = create_mesh_with_field(comm, sendFieldName, sendFieldValue);

  std::string recvFieldName = "recv_field_1";
  double recvFieldValue = 4.321;
  auto recvBulk = create_mesh_with_field(comm, recvFieldName, recvFieldValue, stk::topology::ELEM_RANK);

  stk::transfer_util::TransferMainSettings settings;
  std::vector<std::pair<std::string, std::string>> fieldList = {{sendFieldName, recvFieldName}};
  settings.set_transfer_field(fieldList[0]);

  stk::transfer_util::TransferMainBroker broker(comm, sendBulk, recvBulk, settings);

  broker.check_fields();
  broker.check_parts();

  check_field_values(sendBulk, sendFieldName, sendFieldValue);

  std::vector<std::string> recvFieldsForTransfer = broker.get_recv_fields_for_transfer();
  std::vector<std::string> expectedFields = {recvFieldName};
  check_recv_tranfer_fields(recvFieldsForTransfer, expectedFields);

  broker.transfer();

  check_field_values(sendBulk, sendFieldName, sendFieldValue);
  check_field_values(recvBulk, recvFieldName, sendFieldValue);
  check_field_values(recvBulk, recvFieldName, recvFieldValue, stk::topology::ELEM_RANK);

}

TEST(TransferMainBroker, basicTransferDifferentSendRecvFieldNameOnlySendFieldInList)
{
  stk::ParallelMachine comm = stk::parallel_machine_world();
  if(stk::parallel_machine_size(comm) != 1) { GTEST_SKIP(); }

  std::string sendFieldName = "send_field_1";
  double sendFieldValue = 1.234;
  auto sendBulk = create_mesh_with_field(comm, sendFieldName, sendFieldValue);

  std::string recvFieldName = "recv_field_1";
  double recvFieldValue = 4.321;
  auto recvBulk = create_mesh_with_field(comm, recvFieldName, recvFieldValue);

  stk::transfer_util::TransferMainSettings settings;
  std::vector<std::pair<std::string, std::string>> fieldList = {{sendFieldName, sendFieldName}};
  settings.set_transfer_field(fieldList[0]);

  stk::transfer_util::TransferMainBroker broker(comm, sendBulk, recvBulk, settings);

  broker.check_fields();
  broker.check_parts();

  check_field_values(sendBulk, sendFieldName, sendFieldValue);
  check_field_values(recvBulk, recvFieldName, recvFieldValue);

  std::vector<std::string> recvFieldsForTransfer = broker.get_recv_fields_for_transfer();
  std::vector<std::string> expectedFields = {sendFieldName};
  check_recv_tranfer_fields(recvFieldsForTransfer, expectedFields);

  broker.transfer();

  check_field_values(sendBulk, sendFieldName, sendFieldValue);
  check_field_values(recvBulk, sendFieldName, sendFieldValue);
  check_field_values(recvBulk, recvFieldName, recvFieldValue);

}

TEST(TransferMainBroker, sendFieldNotInFieldList)
{
  stk::ParallelMachine comm = stk::parallel_machine_world();
  if(stk::parallel_machine_size(comm) != 1) { GTEST_SKIP(); }

  std::string sendFieldName = "send_field_1";
  double sendFieldValue = 1.234;
  auto sendBulk = create_mesh_with_field(comm, sendFieldName, sendFieldValue);

  std::string recvFieldName = "recv_field_1";
  double recvFieldValue = 4.321;
  auto recvBulk = create_mesh_with_field(comm, recvFieldName, recvFieldValue);

  stk::transfer_util::TransferMainSettings settings;
  std::vector<std::pair<std::string, std::string>> fieldList = {{recvFieldName, recvFieldName}};
  settings.set_transfer_field(fieldList[0]);

  stk::transfer_util::TransferMainBroker broker(comm, sendBulk, recvBulk, settings);

  EXPECT_ANY_THROW(broker.check_fields());

}

TEST(TransferMainBroker, sendFieldNotTransient)
{
  stk::ParallelMachine comm = stk::parallel_machine_world();
  if(stk::parallel_machine_size(comm) != 1) { GTEST_SKIP(); }

  std::string sendFieldName = "field_1";
  double sendFieldValue = 1.234;
  auto sendBulk = create_mesh_with_field(comm, sendFieldName, sendFieldValue,
                                         stk::topology::NODE_RANK, Ioss::Field::MESH);

  std::string recvFieldName = "field_1";
  double recvFieldValue = 4.321;
  auto recvBulk = create_mesh_with_field(comm, recvFieldName, recvFieldValue);

  stk::transfer_util::TransferMainSettings settings;
  std::vector<std::pair<std::string, std::string>> fieldList = {{sendFieldName, recvFieldName}};
  settings.set_transfer_field(fieldList[0]);

  stk::transfer_util::TransferMainBroker broker(comm, sendBulk, recvBulk, settings);

  EXPECT_ANY_THROW(broker.check_fields());

}

TEST(TransferMainBroker, sendFieldNotTransientNoFieldList)
{
  stk::ParallelMachine comm = stk::parallel_machine_world();
  if(stk::parallel_machine_size(comm) != 1) { GTEST_SKIP(); }

  std::string sendFieldName = "field_1";
  double sendFieldValue = 1.234;
  auto sendBulk = create_mesh_with_field(comm, sendFieldName, sendFieldValue,
                                         stk::topology::NODE_RANK, Ioss::Field::MESH);

  std::string recvFieldName = "field_1";
  double recvFieldValue = 4.321;
  auto recvBulk = create_mesh_with_field(comm, recvFieldName, recvFieldValue);

  stk::transfer_util::TransferMainSettings settings;

  stk::transfer_util::TransferMainBroker broker(comm, sendBulk, recvBulk, settings);

  EXPECT_ANY_THROW(broker.check_fields());

}

TEST(TransferMainBroker, recvFieldNotTransientSameSendandReceive)
{
  stk::ParallelMachine comm = stk::parallel_machine_world();
  if(stk::parallel_machine_size(comm) != 1) { GTEST_SKIP(); }

  std::string sendFieldName = "field_1";
  double sendFieldValue = 1.234;
  auto sendBulk = create_mesh_with_field(comm, sendFieldName, sendFieldValue);

  std::string recvFieldName = "field_1";
  double recvFieldValue = 4.321;
  auto recvBulk = create_mesh_with_field(comm, recvFieldName, recvFieldValue,
                                         stk::topology::NODE_RANK, Ioss::Field::MESH);

  stk::transfer_util::TransferMainSettings settings;
  std::vector<std::pair<std::string, std::string>> fieldList = {{sendFieldName, recvFieldName}};
  settings.set_transfer_field(fieldList[0]);

  stk::transfer_util::TransferMainBroker broker(comm, sendBulk, recvBulk, settings);

  EXPECT_ANY_THROW(broker.check_fields());

}

TEST(TransferMainBroker, recvFieldNotTransientSameSendandReceiveNoFieldList)
{
  stk::ParallelMachine comm = stk::parallel_machine_world();
  if(stk::parallel_machine_size(comm) != 1) { GTEST_SKIP(); }

  std::string sendFieldName = "field_1";
  double sendFieldValue = 1.234;
  auto sendBulk = create_mesh_with_field(comm, sendFieldName, sendFieldValue);

  std::string recvFieldName = "field_1";
  double recvFieldValue = 4.321;
  auto recvBulk = create_mesh_with_field(comm, recvFieldName, recvFieldValue,
                                         stk::topology::NODE_RANK, Ioss::Field::MESH);

  stk::transfer_util::TransferMainSettings settings;

  stk::transfer_util::TransferMainBroker broker(comm, sendBulk, recvBulk, settings);

  EXPECT_ANY_THROW(broker.check_fields());

}

TEST(TransferMainBroker, recvFieldNotTransientDifferentSendandReceive)
{
  stk::ParallelMachine comm = stk::parallel_machine_world();
  if(stk::parallel_machine_size(comm) != 1) { GTEST_SKIP(); }

  std::string sendFieldName = "send_field_1";
  double sendFieldValue = 1.234;
  auto sendBulk = create_mesh_with_field(comm, sendFieldName, sendFieldValue);

  std::string recvFieldName = "recv_field_1";
  double recvFieldValue = 4.321;
  auto recvBulk = create_mesh_with_field(comm, recvFieldName, recvFieldValue,
                                         stk::topology::NODE_RANK, Ioss::Field::MESH);

  stk::transfer_util::TransferMainSettings settings;
  std::vector<std::pair<std::string, std::string>> fieldList = {{sendFieldName, recvFieldName}};
  settings.set_transfer_field(fieldList[0]);

  stk::transfer_util::TransferMainBroker broker(comm, sendBulk, recvBulk, settings);

  EXPECT_ANY_THROW(broker.check_fields());

}

TEST(TransferMainBroker, recvFieldNotTransientDifferentSendandReceiveNoFieldList)
{
  stk::ParallelMachine comm = stk::parallel_machine_world();
  if(stk::parallel_machine_size(comm) != 1) { GTEST_SKIP(); }

  std::string sendFieldName = "send_field_1";
  double sendFieldValue = 1.234;
  auto sendBulk = create_mesh_with_field(comm, sendFieldName, sendFieldValue);

  std::string recvFieldName = "recv_field_1";
  double recvFieldValue = 4.321;
  auto recvBulk = create_mesh_with_field(comm, recvFieldName, recvFieldValue,
                                         stk::topology::NODE_RANK, Ioss::Field::MESH);

  stk::transfer_util::TransferMainSettings settings;

  stk::transfer_util::TransferMainBroker broker(comm, sendBulk, recvBulk, settings);

  broker.check_fields();
  broker.check_parts();

  check_field_values(sendBulk, sendFieldName, sendFieldValue);
  check_field_values(recvBulk, recvFieldName, recvFieldValue);

  std::vector<std::string> recvFieldsForTransfer = broker.get_recv_fields_for_transfer();
  std::vector<std::string> expectedFields = {sendFieldName};
  check_recv_tranfer_fields(recvFieldsForTransfer, expectedFields);

  broker.transfer();

  check_field_values(sendBulk, sendFieldName, sendFieldValue);
  check_field_values(recvBulk, sendFieldName, sendFieldValue);
  check_field_values(recvBulk, recvFieldName, recvFieldValue);

}

TEST(TransferMainBroker, sendMeshTets)
{
  stk::ParallelMachine comm = stk::parallel_machine_world();
  if(stk::parallel_machine_size(comm) != 1) { GTEST_SKIP(); }

  std::string sendFieldName = "field_1";
  double sendFieldValue = 1.234;
  auto sendBulk = create_mesh_with_field(comm, sendFieldName, sendFieldValue,
                                         stk::topology::NODE_RANK, Ioss::Field::TRANSIENT,
                                         "generated:4x1x1|tets");

  std::string recvFieldName = "field_1";
  double recvFieldValue = 4.321;
  auto recvBulk = create_mesh_with_field(comm, recvFieldName, recvFieldValue);

  stk::transfer_util::TransferMainSettings settings;

  stk::transfer_util::TransferMainBroker broker(comm, sendBulk, recvBulk, settings);

  broker.check_fields();
  broker.check_parts();

  check_field_values(sendBulk, sendFieldName, sendFieldValue);
  check_field_values(recvBulk, recvFieldName, recvFieldValue);

  EXPECT_ANY_THROW(broker.transfer());

}

TEST(TransferMainBroker, recvMeshTets)
{
  stk::ParallelMachine comm = stk::parallel_machine_world();
  if(stk::parallel_machine_size(comm) != 1) { GTEST_SKIP(); }

  std::string sendFieldName = "send_field_1";
  double sendFieldValue = 1.234;
  auto sendBulk = create_mesh_with_field(comm, sendFieldName, sendFieldValue);

  std::string recvFieldName = "recv_field_1";
  double recvFieldValue = 4.321;
  auto recvBulk = create_mesh_with_field(comm, recvFieldName, recvFieldValue,
                                         stk::topology::NODE_RANK, Ioss::Field::TRANSIENT,
                                         "generated:4x1x1|tets");

  stk::transfer_util::TransferMainSettings settings;

  stk::transfer_util::TransferMainBroker broker(comm, sendBulk, recvBulk, settings);

  broker.check_fields();
  broker.check_parts();

  check_field_values(sendBulk, sendFieldName, sendFieldValue);
  check_field_values(recvBulk, recvFieldName, recvFieldValue);

  std::vector<std::string> recvFieldsForTransfer = broker.get_recv_fields_for_transfer();
  std::vector<std::string> expectedFields = {sendFieldName};
  check_recv_tranfer_fields(recvFieldsForTransfer, expectedFields);

  broker.transfer();
  
  check_field_values(sendBulk, sendFieldName, sendFieldValue);
  check_field_values(recvBulk, recvFieldName, recvFieldValue);
  check_field_values(recvBulk, sendFieldName, sendFieldValue);

}

std::shared_ptr<stk::mesh::BulkData>
create_mesh_with_fields(const stk::ParallelMachine comm,
                       const std::vector<std::string>& fieldNames,
                       const std::vector<double>& fieldValues,
                       const stk::mesh::EntityRank fieldRank = stk::topology::NODE_RANK,
                       const Ioss::Field::RoleType fieldRole = Ioss::Field::TRANSIENT)
{
  std::shared_ptr<stk::mesh::BulkData> bulk = stk::unit_test_util::build_mesh(3u, comm);
  for (size_t i = 0; i < fieldNames.size(); i++) {
    stk::mesh::Field<double> &field=
          bulk->mesh_meta_data().declare_field<double>(fieldRank, fieldNames[i], 1);
    stk::io::set_field_role(field, fieldRole);
    stk::mesh::put_field_on_mesh(field, bulk->mesh_meta_data().universal_part(), &fieldValues[i]);
  }
  stk::io::fill_mesh("generated:4x1x1", *bulk);

  return bulk;
}

TEST(TransferMainBroker, basicTransferTwoFieldsBothTransferred)
{
  stk::ParallelMachine comm = stk::parallel_machine_world();
  if(stk::parallel_machine_size(comm) != 1) { GTEST_SKIP(); }

  std::vector<std::string> sendFieldNames = {"field_1", "field_2"};
  std::vector<double> sendFieldValues = {1.234, 5.678};
  auto sendBulk = create_mesh_with_fields(comm, sendFieldNames, sendFieldValues);

  std::vector<std::string> recvFieldNames = {"field_1", "field_2"};
  std::vector<double> recvFieldValues = {4.321, 8.765};
  auto recvBulk = create_mesh_with_fields(comm, recvFieldNames, recvFieldValues);

  stk::transfer_util::TransferMainSettings settings;
  std::vector<std::pair<std::string, std::string>> fieldList = {{sendFieldNames[0], recvFieldNames[0]},
                                                                {sendFieldNames[1], recvFieldNames[1]}};
  settings.set_transfer_field(fieldList[0]);
  settings.set_transfer_field(fieldList[1]);

  stk::transfer_util::TransferMainBroker broker(comm, sendBulk, recvBulk, settings);

  broker.check_fields();
  broker.check_parts();

  check_field_values(sendBulk, sendFieldNames[0], sendFieldValues[0]);
  check_field_values(sendBulk, sendFieldNames[1], sendFieldValues[1]);
  check_field_values(recvBulk, recvFieldNames[0], recvFieldValues[0]);
  check_field_values(recvBulk, recvFieldNames[1], recvFieldValues[1]);

  std::vector<std::string> recvFieldsForTransfer = broker.get_recv_fields_for_transfer();
  check_recv_tranfer_fields(recvFieldsForTransfer, recvFieldNames);

  broker.transfer();

  check_field_values(sendBulk, sendFieldNames[0], sendFieldValues[0]);
  check_field_values(sendBulk, sendFieldNames[1], sendFieldValues[1]);
  check_field_values(recvBulk, recvFieldNames[0], sendFieldValues[0]);
  check_field_values(recvBulk, recvFieldNames[1], sendFieldValues[1]);

}

TEST(TransferMainBroker, basicTransferTwoFieldsBothTransferredNoFieldList)
{
  stk::ParallelMachine comm = stk::parallel_machine_world();
  if(stk::parallel_machine_size(comm) != 1) { GTEST_SKIP(); }

  std::vector<std::string> sendFieldNames = {"field_1", "field_2"};
  std::vector<double> sendFieldValues = {1.234, 5.678};
  auto sendBulk = create_mesh_with_fields(comm, sendFieldNames, sendFieldValues);

  std::vector<std::string> recvFieldNames = {"field_1", "field_2"};
  std::vector<double> recvFieldValues = {4.321, 8.765};
  auto recvBulk = create_mesh_with_fields(comm, recvFieldNames, recvFieldValues);

  stk::transfer_util::TransferMainSettings settings;

  stk::transfer_util::TransferMainBroker broker(comm, sendBulk, recvBulk, settings);

  broker.check_fields();
  broker.check_parts();

  check_field_values(sendBulk, sendFieldNames[0], sendFieldValues[0]);
  check_field_values(sendBulk, sendFieldNames[1], sendFieldValues[1]);
  check_field_values(recvBulk, recvFieldNames[0], recvFieldValues[0]);
  check_field_values(recvBulk, recvFieldNames[1], recvFieldValues[1]);

  std::vector<std::string> recvFieldsForTransfer = broker.get_recv_fields_for_transfer();
  check_recv_tranfer_fields(recvFieldsForTransfer, recvFieldNames);

  broker.transfer();

  check_field_values(sendBulk, sendFieldNames[0], sendFieldValues[0]);
  check_field_values(sendBulk, sendFieldNames[1], sendFieldValues[1]);
  check_field_values(recvBulk, recvFieldNames[0], sendFieldValues[0]);
  check_field_values(recvBulk, recvFieldNames[1], sendFieldValues[1]);

}

TEST(TransferMainBroker, basicTransferTwoFieldsOneTransferred)
{
  stk::ParallelMachine comm = stk::parallel_machine_world();
  if(stk::parallel_machine_size(comm) != 1) { GTEST_SKIP(); }

  std::vector<std::string> sendFieldNames = {"send_field_1", "send_field_2"};
  std::vector<double> sendFieldValues = {1.234, 5.678};
  auto sendBulk = create_mesh_with_fields(comm, sendFieldNames, sendFieldValues);

  std::vector<std::string> recvFieldNames = {"recv_field_1", "recv_field_2"};
  std::vector<double> recvFieldValues = {4.321, 8.765};
  auto recvBulk = create_mesh_with_fields(comm, recvFieldNames, recvFieldValues);

  stk::transfer_util::TransferMainSettings settings;
  std::vector<std::pair<std::string, std::string>> fieldList = {{sendFieldNames[1], recvFieldNames[1]}};
  settings.set_transfer_field(fieldList[0]);

  stk::transfer_util::TransferMainBroker broker(comm, sendBulk, recvBulk, settings);

  broker.check_fields();
  broker.check_parts();

  check_field_values(sendBulk, sendFieldNames[0], sendFieldValues[0]);
  check_field_values(sendBulk, sendFieldNames[1], sendFieldValues[1]);
  check_field_values(recvBulk, recvFieldNames[0], recvFieldValues[0]);
  check_field_values(recvBulk, recvFieldNames[1], recvFieldValues[1]);

  std::vector<std::string> recvFieldsForTransfer = broker.get_recv_fields_for_transfer();
  std::vector<std::string> expectedFields = {recvFieldNames[1]};
  check_recv_tranfer_fields(recvFieldsForTransfer, expectedFields);

  broker.transfer();

  check_field_values(sendBulk, sendFieldNames[0], sendFieldValues[0]);
  check_field_values(sendBulk, sendFieldNames[1], sendFieldValues[1]);
  check_field_values(recvBulk, recvFieldNames[0], recvFieldValues[0]);
  check_field_values(recvBulk, recvFieldNames[1], sendFieldValues[1]);

}


}
