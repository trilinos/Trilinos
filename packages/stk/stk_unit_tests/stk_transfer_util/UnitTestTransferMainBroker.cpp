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
#include "stk_unit_test_utils/meshCreationHelpers.hpp"
#include "stk_unit_test_utils/CommandLineArgs.hpp"
#include "stk_unit_test_utils/BuildMesh.hpp"
#include "stk_mesh/base/ForEachEntity.hpp"

namespace {

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
      EXPECT_NEAR(expectedValue, dataForNode(i), 1.0e-6);
    }
  };

  stk::mesh::for_each_entity_run(*bulk, fieldRank, checkFieldValues);
}

void check_field_values(const std::shared_ptr<stk::mesh::BulkData> bulk, 
                        const std::string& fieldName,
                        const stk::unit_test_util::FieldEvaluator& expectedValue,
                        const stk::mesh::EntityRank fieldRank = stk::topology::NODE_RANK)
{
  auto field = bulk->mesh_meta_data().get_field(fieldRank, fieldName);
  const auto* coordField = bulk->mesh_meta_data().coordinate_field();
  auto fieldData = field->data<double, stk::mesh::ReadOnly>();
  auto coordFieldData = coordField->data<double, stk::mesh::ReadOnly>();
  
  auto checkFieldValues = [&](const stk::mesh::BulkData& , stk::mesh::Entity node) {
    auto dataForNode = fieldData.entity_values(node);
    auto nodeCoord = coordFieldData.entity_values(node);
    for(stk::mesh::ComponentIdx i=0_comp; i<dataForNode.num_components(); ++i) {
      double x = nodeCoord(0_comp);
      double y = nodeCoord(1_comp);
      double z = nodeCoord.num_components()==3 ? nodeCoord(2_comp) : 0.0;
      EXPECT_NEAR(expectedValue(node,x,y,z), dataForNode(i), 1.0e-6);
    }
  };

  stk::mesh::for_each_entity_run(*bulk, fieldRank, checkFieldValues);
  
}

template<class CoordCondition>
void check_field_values(const std::shared_ptr<stk::mesh::BulkData> bulk, 
                        const std::string& fieldName,
                        const stk::unit_test_util::FieldEvaluator& expectedValue,
                        const CoordCondition& condition,
                        const stk::mesh::EntityRank fieldRank = stk::topology::NODE_RANK)
{
  auto* field = bulk->mesh_meta_data().get_field(fieldRank, fieldName);
  const auto* coordField = bulk->mesh_meta_data().coordinate_field();
  auto fieldData = field->data<double, stk::mesh::ReadOnly>();
  auto coordFieldData = coordField->data<double, stk::mesh::ReadOnly>();
  
  auto checkFieldValues = [&](const stk::mesh::BulkData& mesh, stk::mesh::Entity node) {
    auto dataForNode = fieldData.entity_values(node);
    auto nodeCoord = coordFieldData.entity_values(node);
    for(stk::mesh::ComponentIdx i=0_comp; i<dataForNode.num_components(); ++i) {
      double x = nodeCoord(0_comp);
      double y = nodeCoord(1_comp);
      double z = nodeCoord.num_components()==3 ? nodeCoord(2_comp) : 0.0;
      if (condition(x, y, z)) {
        EXPECT_NEAR(expectedValue(node, x,y,z), dataForNode(i), 1.0e-6)<<"nodeId="<<mesh.identifier(node)<<" x="<<x<<" y="<<y<<" z="<<z;
      }
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
  auto sendBulk = stk::unit_test_util::create_mesh_with_field(comm, sendFieldName, sendFieldValue);

  std::string recvFieldName = "field_1";
  double recvFieldValue = 4.321;
  auto recvBulk = stk::unit_test_util::create_mesh_with_field(comm, recvFieldName, recvFieldValue);

  stk::transfer_util::TransferMainSettings settings;
  std::vector<std::pair<std::string, std::string>> fieldList = {{sendFieldName, recvFieldName}};
  settings.set_transfer_field(fieldList[0]);

  stk::transfer_util::TransferMainBroker broker(comm, sendBulk, recvBulk, settings);

  EXPECT_TRUE(broker.get_send_fields_for_transfer().empty());
  EXPECT_TRUE(broker.get_recv_fields_for_transfer().empty());

  broker.check_and_create_fields();

  auto sendFieldNames = broker.get_send_fields_for_transfer();
  ASSERT_EQ(1u, sendFieldNames.size());
  EXPECT_EQ(sendFieldName, sendFieldNames[0]);

  auto recvFieldNames = broker.get_recv_fields_for_transfer();
  ASSERT_EQ(1u, recvFieldNames.size());
  EXPECT_EQ(recvFieldName, recvFieldNames[0]);

  broker.check_parts();

  check_field_values(sendBulk, sendFieldName, sendFieldValue);
  check_field_values(recvBulk, recvFieldName, recvFieldValue);

  std::vector<std::string> recvFieldsForTransfer = broker.get_recv_fields_for_transfer();
  std::vector<std::string> expectedFields = {recvFieldName};
  check_recv_tranfer_fields(recvFieldsForTransfer, expectedFields);

  broker.transfer_initialize();
  broker.transfer_apply();

  check_field_values(sendBulk, sendFieldName, sendFieldValue);
  check_field_values(recvBulk, recvFieldName, sendFieldValue);

}

TEST(TransferMainBroker, basicTransfer_3D_to_2D)
{
  stk::ParallelMachine comm = stk::parallel_machine_world();
  if(stk::parallel_machine_size(comm) != 1) { GTEST_SKIP(); }

  std::string meshSpec3D = "generated:2x2x8|bbox:0,0,0,2,2,8";
  std::string sendFieldName = "field_1";
  constexpr int numCopies = 1;
  constexpr int spatialDim = 3;
  stk::unit_test_util::PlanarLinearFieldEvaluator fieldEval(spatialDim);
  auto sendBulk = stk::unit_test_util::create_mesh_with_field(comm,
                        sendFieldName, fieldEval, stk::topology::NODE_RANK,
                        Ioss::Field::TRANSIENT, meshSpec3D, numCopies, spatialDim);

  std::string meshSpec2D = "textmesh:0,1,QUAD_4_2D,1,2,5,4\n"
                                    "0,2,QUAD_4_2D,2,3,6,5\n"
                                    "0,3,QUAD_4_2D,4,5,8,7\n"
                                    "0,4,QUAD_4_2D,5,6,9,8"
          "|coordinates: 0,0, 1,0, 2,0,  0,1, 1,1, 2,1,  0,2, 1,2, 2,2 "
          "|dimension:2";
  std::string recvFieldName = "field_1";
  double recvFieldValue = 4.321;
  auto recvBulk = stk::unit_test_util::create_mesh_with_field_2d(comm, recvFieldName, recvFieldValue, stk::topology::NODE_RANK, Ioss::Field::TRANSIENT, meshSpec2D);

  stk::transfer_util::TransferMainSettings settings;
  std::vector<std::pair<std::string, std::string>> fieldList = {{sendFieldName, recvFieldName}};
  settings.set_transfer_field(fieldList[0]);

  stk::transfer_util::TransferMainBroker broker(comm, sendBulk, recvBulk, settings);

  broker.check_and_create_fields();
  broker.check_parts();

  check_field_values(sendBulk, sendFieldName, fieldEval);
  check_field_values(recvBulk, recvFieldName, recvFieldValue);

  broker.transfer_initialize();
  broker.transfer_apply();

  check_field_values(sendBulk, sendFieldName, fieldEval);
  check_field_values(recvBulk, recvFieldName, fieldEval);
}

TEST(TransferMainBroker, basicTransfer_2D_to_3D_extrude)
{
  stk::ParallelMachine comm = stk::parallel_machine_world();
  if(stk::parallel_machine_size(comm) != 1) { GTEST_SKIP(); }

  std::string meshSpec2D = "textmesh:0,1,QUAD_4_2D,1,2,5,4\n"
                                    "0,2,QUAD_4_2D,2,3,6,5\n"
                                    "0,3,QUAD_4_2D,4,5,8,7\n"
                                    "0,4,QUAD_4_2D,5,6,9,8"
          "|coordinates: 0,0, 1,0, 2,0,  0,1, 1,1, 2,1,  0,2, 1,2, 2,2 "
          "|dimension:2";
  std::string sendFieldName = "field_1";
  constexpr int numCopies = 1;
  constexpr int spatialDim = 2;
  stk::unit_test_util::LinearFieldEvaluator linearFieldEval(spatialDim);
  auto sendBulk = stk::unit_test_util::create_mesh_with_field(comm,
                        sendFieldName, linearFieldEval, stk::topology::NODE_RANK,
                        Ioss::Field::TRANSIENT, meshSpec2D, numCopies, spatialDim);

  std::string meshSpec3D = "generated:2x2x8|bbox:0,0,0,2,2,8";
  std::string recvFieldName = "field_1";
  double recvFieldValue = 4.321;
  auto recvBulk = stk::unit_test_util::create_mesh_with_field(comm, recvFieldName, recvFieldValue, stk::topology::NODE_RANK, Ioss::Field::TRANSIENT, meshSpec3D);

  stk::transfer_util::TransferMainSettings settings;
  std::vector<std::pair<std::string, std::string>> fieldList = {{sendFieldName, recvFieldName}};
  settings.set_transfer_field(fieldList[0]);

  settings.set_2d_to_3d_mapping_type("EXTRUDE");

  stk::transfer_util::TransferMainBroker broker(comm, sendBulk, recvBulk, settings);

  broker.check_and_create_fields();
  broker.check_parts();

  check_field_values(sendBulk, sendFieldName, linearFieldEval);
  check_field_values(recvBulk, recvFieldName, recvFieldValue);

  broker.transfer_initialize();
  broker.transfer_apply();

  check_field_values(sendBulk, sendFieldName, linearFieldEval);
  check_field_values(recvBulk, recvFieldName, linearFieldEval);
}

TEST(TransferMainBroker, basicTransfer_2D_to_3D_zequal0)
{
  stk::ParallelMachine comm = stk::parallel_machine_world();
  if(stk::parallel_machine_size(comm) != 1) { GTEST_SKIP(); }

  std::string meshSpec2D = "textmesh:0,1,QUAD_4_2D,1,2,5,4\n"
                                    "0,2,QUAD_4_2D,2,3,6,5\n"
                                    "0,3,QUAD_4_2D,4,5,8,7\n"
                                    "0,4,QUAD_4_2D,5,6,9,8"
          "|coordinates: 0,0, 1,0, 2,0,  0,1, 1,1, 2,1,  0,2, 1,2, 2,2 "
          "|dimension:2";
  std::string sendFieldName = "field_1";
  constexpr int numCopies = 1;
  constexpr int spatialDim = 2;
  stk::unit_test_util::LinearFieldEvaluator linearFieldEval(spatialDim);
  auto sendBulk = stk::unit_test_util::create_mesh_with_field(comm,
                        sendFieldName, linearFieldEval, stk::topology::NODE_RANK,
                        Ioss::Field::TRANSIENT, meshSpec2D, numCopies, spatialDim);

  std::string meshSpec3D = "generated:2x2x8|bbox:0,0,0,2,2,8";
  std::string recvFieldName = "field_1";
  double recvFieldValue = 0.0;
  auto recvBulk = stk::unit_test_util::create_mesh_with_field(comm, recvFieldName, recvFieldValue, stk::topology::NODE_RANK, Ioss::Field::TRANSIENT, meshSpec3D);

  stk::transfer_util::TransferMainSettings settings;
  std::vector<std::pair<std::string, std::string>> fieldList = {{sendFieldName, recvFieldName}};
  settings.set_transfer_field(fieldList[0]);

  settings.set_2d_to_3d_mapping_type("ZPLANE");

  stk::transfer_util::TransferMainBroker broker(comm, sendBulk, recvBulk, settings);

  broker.check_and_create_fields();
  broker.check_parts();

  check_field_values(sendBulk, sendFieldName, linearFieldEval);
  check_field_values(recvBulk, recvFieldName, recvFieldValue);

  broker.transfer_initialize();
  broker.transfer_apply();

  check_field_values(sendBulk, sendFieldName, linearFieldEval);

  auto zEquals0 = [](double /*x*/, double /*y*/, double z) -> bool
  { return std::abs(z) < std::numeric_limits<double>::epsilon(); };

  check_field_values(recvBulk, recvFieldName, linearFieldEval, zEquals0);
}

TEST(TransferMainBroker, basicTransfer_2D_to_3D_zequalConstant)
{
  stk::ParallelMachine comm = stk::parallel_machine_world();
  if(stk::parallel_machine_size(comm) != 1) { GTEST_SKIP(); }

  std::string meshSpec2D = "textmesh:0,1,QUAD_4_2D,1,2,5,4\n"
                                    "0,2,QUAD_4_2D,2,3,6,5\n"
                                    "0,3,QUAD_4_2D,4,5,8,7\n"
                                    "0,4,QUAD_4_2D,5,6,9,8"
          "|coordinates: 0,0, 1,0, 2,0,  0,1, 1,1, 2,1,  0,2, 1,2, 2,2 "
          "|dimension:2";
  std::string sendFieldName = "field_1";
  constexpr int numCopies = 1;
  constexpr int spatialDim = 2;
  stk::unit_test_util::LinearFieldEvaluator linearFieldEval(spatialDim);
  auto sendBulk = stk::unit_test_util::create_mesh_with_field(comm,
                        sendFieldName, linearFieldEval, stk::topology::NODE_RANK,
                        Ioss::Field::TRANSIENT, meshSpec2D, numCopies, spatialDim);

  std::string meshSpec3D = "generated:2x2x8|bbox:0,0,0,2,2,8";
  std::string recvFieldName = "field_1";
  double recvFieldValue = 0.0;
  auto recvBulk = stk::unit_test_util::create_mesh_with_field(comm, recvFieldName, recvFieldValue, stk::topology::NODE_RANK, Ioss::Field::TRANSIENT, meshSpec3D);

  stk::transfer_util::TransferMainSettings settings;
  std::vector<std::pair<std::string, std::string>> fieldList = {{sendFieldName, recvFieldName}};
  settings.set_transfer_field(fieldList[0]);

  settings.set_2d_to_3d_mapping_type("ZPLANE");
  settings.set_coord_transf_z_expr("z=1");

  stk::transfer_util::TransferMainBroker broker(comm, sendBulk, recvBulk, settings);

  broker.check_and_create_fields();
  broker.check_parts();

  check_field_values(sendBulk, sendFieldName, linearFieldEval);
  check_field_values(recvBulk, recvFieldName, recvFieldValue);

  broker.transfer_initialize();
  broker.transfer_apply();

  check_field_values(sendBulk, sendFieldName, linearFieldEval);

  auto zEquals1 = [](double /*x*/, double /*y*/, double z) -> bool
  { return std::abs(z) == 1.0; };

  check_field_values(recvBulk, recvFieldName, linearFieldEval, zEquals1);
}

TEST(TransferMainBroker, basicTransferDifferentSendRecvFieldName)
{
  stk::ParallelMachine comm = stk::parallel_machine_world();
  if(stk::parallel_machine_size(comm) != 1) { GTEST_SKIP(); }

  std::string sendFieldName = "send_field_1";
  double sendFieldValue = 1.234;
  auto sendBulk = stk::unit_test_util::create_mesh_with_field(comm, sendFieldName, sendFieldValue);

  std::string recvFieldName = "recv_field_1";
  double recvFieldValue = 4.321;
  auto recvBulk = stk::unit_test_util::create_mesh_with_field(comm, recvFieldName, recvFieldValue);

  stk::transfer_util::TransferMainSettings settings;
  std::vector<std::pair<std::string, std::string>> fieldList = {{sendFieldName, recvFieldName}};
  settings.set_transfer_field(fieldList[0]);

  stk::transfer_util::TransferMainBroker broker(comm, sendBulk, recvBulk, settings);

  broker.check_and_create_fields();
  broker.check_parts();

  check_field_values(sendBulk, sendFieldName, sendFieldValue);
  check_field_values(recvBulk, recvFieldName, recvFieldValue);

  std::vector<std::string> recvFieldsForTransfer = broker.get_recv_fields_for_transfer();
  std::vector<std::string> expectedFields = {recvFieldName};
  check_recv_tranfer_fields(recvFieldsForTransfer, expectedFields);

  broker.transfer_initialize();
  broker.transfer_apply();

  check_field_values(sendBulk, sendFieldName, sendFieldValue);
  check_field_values(recvBulk, recvFieldName, sendFieldValue);

}

TEST(TransferMainBroker, basicTransferSendFieldOnly)
{
  stk::ParallelMachine comm = stk::parallel_machine_world();
  if(stk::parallel_machine_size(comm) != 1) { GTEST_SKIP(); }

  std::string sendFieldName = "send_field_1";
  double sendFieldValue = 1.234;
  auto sendBulk = stk::unit_test_util::create_mesh_with_field(comm, sendFieldName, sendFieldValue);

  std::shared_ptr<stk::mesh::BulkData> recvBulk = stk::unit_test_util::build_mesh(3u, comm);
  stk::io::fill_mesh("generated:4x1x1", *recvBulk);

  stk::transfer_util::TransferMainSettings settings;
  std::vector<std::pair<std::string, std::string>> fieldList = {{sendFieldName, sendFieldName}};
  settings.set_transfer_field(fieldList[0]);

  stk::transfer_util::TransferMainBroker broker(comm, sendBulk, recvBulk, settings);

  broker.check_and_create_fields();
  broker.check_parts();

  check_field_values(sendBulk, sendFieldName, sendFieldValue);

  std::vector<std::string> recvFieldsForTransfer = broker.get_recv_fields_for_transfer();
  std::vector<std::string> expectedFields = {sendFieldName};
  check_recv_tranfer_fields(recvFieldsForTransfer, expectedFields);

  broker.transfer_initialize();
  broker.transfer_apply();

  check_field_values(sendBulk, sendFieldName, sendFieldValue);
  check_field_values(recvBulk, sendFieldName, sendFieldValue);

}

TEST(TransferMainBroker, basicTransferNoFieldListSendRecvFieldName)
{
  stk::ParallelMachine comm = stk::parallel_machine_world();
  if(stk::parallel_machine_size(comm) != 1) { GTEST_SKIP(); }

  std::string sendFieldName = "field_1";
  double sendFieldValue = 1.234;
  auto sendBulk = stk::unit_test_util::create_mesh_with_field(comm, sendFieldName, sendFieldValue);

  std::string recvFieldName = "field_1";
  double recvFieldValue = 4.321;
  auto recvBulk = stk::unit_test_util::create_mesh_with_field(comm, recvFieldName, recvFieldValue);

  stk::transfer_util::TransferMainSettings settings;

  stk::transfer_util::TransferMainBroker broker(comm, sendBulk, recvBulk, settings);

  broker.check_and_create_fields();
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
  auto sendBulk = stk::unit_test_util::create_mesh_with_field(comm, sendFieldName, sendFieldValue);

  std::shared_ptr<stk::mesh::BulkData> recvBulk = stk::unit_test_util::build_mesh(3u, comm);
  stk::io::fill_mesh("generated:4x1x1", *recvBulk);

  stk::transfer_util::TransferMainSettings settings;

  stk::transfer_util::TransferMainBroker broker(comm, sendBulk, recvBulk, settings);

  broker.check_and_create_fields();
  broker.check_parts();

  check_field_values(sendBulk, sendFieldName, sendFieldValue);

  std::vector<std::string> recvFieldsForTransfer = broker.get_recv_fields_for_transfer();
  std::vector<std::string> expectedFields = {sendFieldName};
  check_recv_tranfer_fields(recvFieldsForTransfer, expectedFields);

  broker.transfer();

  check_field_values(sendBulk, sendFieldName, sendFieldValue);
  check_field_values(recvBulk, sendFieldName, sendFieldValue);

}

TEST(TransferMainBroker, basicTransferNoFieldListSendFieldOnlyRecvFieldDiffRank)
{
  stk::ParallelMachine comm = stk::parallel_machine_world();
  if(stk::parallel_machine_size(comm) != 1) { GTEST_SKIP(); }

  std::string sendFieldName = "field_1";
  double sendFieldValue = 1.234;
  auto sendBulk = stk::unit_test_util::create_mesh_with_field(comm, sendFieldName, sendFieldValue);

  std::shared_ptr<stk::mesh::BulkData> recvBulk = stk::unit_test_util::build_mesh(3u, comm);
  stk::io::fill_mesh("generated:4x1x1", *recvBulk);

  stk::transfer_util::TransferMainSettings settings;
  settings.set_recv_type("ELEMENT_CENTROID");

  stk::transfer_util::TransferMainBroker broker(comm, sendBulk, recvBulk, settings);

  EXPECT_NO_THROW(broker.check_and_create_fields());

}

TEST(TransferMainBroker, basicTransferNoFieldListSendFieldOnlyRecvFieldDiffRankGaussPoint)
{
  stk::ParallelMachine comm = stk::parallel_machine_world();
  if(stk::parallel_machine_size(comm) != 1) { GTEST_SKIP(); }

  std::string sendFieldName = "field_1";
  double sendFieldValue = 1.234;
  auto sendBulk = stk::unit_test_util::create_mesh_with_field(comm, sendFieldName, sendFieldValue);

  std::shared_ptr<stk::mesh::BulkData> recvBulk = stk::unit_test_util::build_mesh(3u, comm);
  stk::io::fill_mesh("generated:4x1x1", *recvBulk);

  stk::transfer_util::TransferMainSettings settings;
  settings.set_recv_type("ELEMENT_GAUSS_POINT");

  stk::transfer_util::TransferMainBroker broker(comm, sendBulk, recvBulk, settings);

  EXPECT_NO_THROW(broker.check_and_create_fields());

}

TEST(TransferMainBroker, sendPartInvalid)
{
  stk::ParallelMachine comm = stk::parallel_machine_world();
  if(stk::parallel_machine_size(comm) != 1) { GTEST_SKIP(); }

  std::string sendFieldName = "field_1";
  std::vector<std::string> sendPartNames = {"part_1234"};
  double sendFieldValue = 1.234;
  auto sendBulk = stk::unit_test_util::create_mesh_with_field(comm, sendFieldName, sendFieldValue);

  std::string recvFieldName = "field_1";
  std::vector<std::string> defaultRecvPartNames = {"block_1"};
  double recvFieldValue = 4.321;
  auto recvBulk = stk::unit_test_util::create_mesh_with_field(comm, recvFieldName, recvFieldValue);

  stk::transfer_util::TransferMainSettings settings;
  std::vector<std::pair<std::string, std::string>> fieldList = {{sendFieldName, recvFieldName}};
  settings.set_transfer_field(fieldList[0]);
  settings.set_transfer_send_parts(sendPartNames);

  stk::transfer_util::TransferMainBroker broker(comm, sendBulk, recvBulk, settings);

  auto expectedSendPartNames = broker.get_send_parts_for_transfer();
  auto expectedRecvPartNames = broker.get_recv_parts_for_transfer();
  EXPECT_EQ(sendPartNames, expectedSendPartNames);
  EXPECT_EQ(defaultRecvPartNames, expectedRecvPartNames);
  EXPECT_ANY_THROW(broker.check_parts());

}

TEST(TransferMainBroker, recvPartInvalid)
{
  stk::ParallelMachine comm = stk::parallel_machine_world();
  if(stk::parallel_machine_size(comm) != 1) { GTEST_SKIP(); }

  std::string sendFieldName = "field_1";
  std::vector<std::string> defaultSendPartNames = {"block_1"};
  double sendFieldValue = 1.234;
  auto sendBulk = stk::unit_test_util::create_mesh_with_field(comm, sendFieldName, sendFieldValue);

  std::string recvFieldName = "field_1";
  std::vector<std::string> recvPartNames = {"part_1234"};
  double recvFieldValue = 4.321;
  auto recvBulk = stk::unit_test_util::create_mesh_with_field(comm, recvFieldName, recvFieldValue);

  stk::transfer_util::TransferMainSettings settings;
  std::vector<std::pair<std::string, std::string>> fieldList = {{sendFieldName, recvFieldName}};
  settings.set_transfer_field(fieldList[0]);
  settings.set_transfer_recv_parts(recvPartNames);

  stk::transfer_util::TransferMainBroker broker(comm, sendBulk, recvBulk, settings);

  auto expectedSendPartNames = broker.get_send_parts_for_transfer();
  auto expectedRecvPartNames = broker.get_recv_parts_for_transfer();
  EXPECT_EQ(defaultSendPartNames, expectedSendPartNames);
  EXPECT_EQ(recvPartNames, expectedRecvPartNames);
  EXPECT_ANY_THROW(broker.check_parts());

}

TEST(TransferMainBroker, noRecvPart)
{
  stk::ParallelMachine comm = stk::parallel_machine_world();
  if(stk::parallel_machine_size(comm) != 1) { GTEST_SKIP(); }

  std::string sendFieldName = "field_1";
  double sendFieldValue = 1.234;
  auto sendBulk = stk::unit_test_util::create_mesh_with_field(comm, sendFieldName, sendFieldValue);

  std::shared_ptr<stk::mesh::BulkData> recvBulk = stk::unit_test_util::build_mesh(3u, comm);
  stk::io::fill_mesh("generated:4x1x1", *recvBulk);

  stk::transfer_util::TransferMainSettings settings;
  settings.set_recv_type("EDGE_GAUSS_POINT");

  EXPECT_ANY_THROW(stk::transfer_util::TransferMainBroker broker(comm, sendBulk, recvBulk, settings));

}

TEST(TransferMainBroker, basicTransferDifferentSendRecvFieldNameOnlySendFieldInList)
{
  stk::ParallelMachine comm = stk::parallel_machine_world();
  if(stk::parallel_machine_size(comm) != 1) { GTEST_SKIP(); }

  std::string sendFieldName = "send_field_1";
  double sendFieldValue = 1.234;
  auto sendBulk = stk::unit_test_util::create_mesh_with_field(comm, sendFieldName, sendFieldValue);

  std::string recvFieldName = "recv_field_1";
  double recvFieldValue = 4.321;
  auto recvBulk = stk::unit_test_util::create_mesh_with_field(comm, recvFieldName, recvFieldValue);

  stk::transfer_util::TransferMainSettings settings;
  std::vector<std::pair<std::string, std::string>> fieldList = {{sendFieldName, sendFieldName}};
  settings.set_transfer_field(fieldList[0]);

  stk::transfer_util::TransferMainBroker broker(comm, sendBulk, recvBulk, settings);

  broker.check_parts();
  broker.check_and_create_fields();

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
  auto sendBulk = stk::unit_test_util::create_mesh_with_field(comm, sendFieldName, sendFieldValue);

  std::string recvFieldName = "recv_field_1";
  double recvFieldValue = 4.321;
  auto recvBulk = stk::unit_test_util::create_mesh_with_field(comm, recvFieldName, recvFieldValue);

  stk::transfer_util::TransferMainSettings settings;
  std::vector<std::pair<std::string, std::string>> fieldList = {{recvFieldName, recvFieldName}};
  settings.set_transfer_field(fieldList[0]);

  stk::transfer_util::TransferMainBroker broker(comm, sendBulk, recvBulk, settings);

  EXPECT_ANY_THROW(broker.check_and_create_fields());

}

TEST(TransferMainBroker, sendFieldNotTransient)
{
  stk::ParallelMachine comm = stk::parallel_machine_world();
  if(stk::parallel_machine_size(comm) != 1) { GTEST_SKIP(); }

  std::string sendFieldName = "field_1";
  double sendFieldValue = 1.234;
  auto sendBulk = stk::unit_test_util::create_mesh_with_field(comm, sendFieldName, sendFieldValue,
                                         stk::topology::NODE_RANK, Ioss::Field::MESH);

  std::string recvFieldName = "field_1";
  double recvFieldValue = 4.321;
  auto recvBulk = stk::unit_test_util::create_mesh_with_field(comm, recvFieldName, recvFieldValue);

  stk::transfer_util::TransferMainSettings settings;
  std::vector<std::pair<std::string, std::string>> fieldList = {{sendFieldName, recvFieldName}};
  settings.set_transfer_field(fieldList[0]);

  stk::transfer_util::TransferMainBroker broker(comm, sendBulk, recvBulk, settings);

  EXPECT_ANY_THROW(broker.check_and_create_fields());

}

TEST(TransferMainBroker, sendFieldNotTransientNoFieldList)
{
  stk::ParallelMachine comm = stk::parallel_machine_world();
  if(stk::parallel_machine_size(comm) != 1) { GTEST_SKIP(); }

  std::string sendFieldName = "field_1";
  double sendFieldValue = 1.234;
  auto sendBulk = stk::unit_test_util::create_mesh_with_field(comm, sendFieldName, sendFieldValue,
                                         stk::topology::NODE_RANK, Ioss::Field::MESH);

  std::string recvFieldName = "field_1";
  double recvFieldValue = 4.321;
  auto recvBulk = stk::unit_test_util::create_mesh_with_field(comm, recvFieldName, recvFieldValue);

  stk::transfer_util::TransferMainSettings settings;

  stk::transfer_util::TransferMainBroker broker(comm, sendBulk, recvBulk, settings);

  EXPECT_ANY_THROW(broker.check_and_create_fields());

}

TEST(TransferMainBroker, recvFieldNotTransientSameSendandReceive)
{
  stk::ParallelMachine comm = stk::parallel_machine_world();
  if(stk::parallel_machine_size(comm) != 1) { GTEST_SKIP(); }

  std::string sendFieldName = "field_1";
  double sendFieldValue = 1.234;
  auto sendBulk = stk::unit_test_util::create_mesh_with_field(comm, sendFieldName, sendFieldValue);

  std::string recvFieldName = "field_1";
  double recvFieldValue = 4.321;
  auto recvBulk = stk::unit_test_util::create_mesh_with_field(comm, recvFieldName, recvFieldValue,
                                         stk::topology::NODE_RANK, Ioss::Field::MESH);

  stk::transfer_util::TransferMainSettings settings;
  std::vector<std::pair<std::string, std::string>> fieldList = {{sendFieldName, recvFieldName}};
  settings.set_transfer_field(fieldList[0]);

  stk::transfer_util::TransferMainBroker broker(comm, sendBulk, recvBulk, settings);

  EXPECT_ANY_THROW(broker.check_and_create_fields());

}

TEST(TransferMainBroker, recvFieldNotTransientSameSendandReceiveNoFieldList)
{
  stk::ParallelMachine comm = stk::parallel_machine_world();
  if(stk::parallel_machine_size(comm) != 1) { GTEST_SKIP(); }

  std::string sendFieldName = "field_1";
  double sendFieldValue = 1.234;
  auto sendBulk = stk::unit_test_util::create_mesh_with_field(comm, sendFieldName, sendFieldValue);

  std::string recvFieldName = "field_1";
  double recvFieldValue = 4.321;
  auto recvBulk = stk::unit_test_util::create_mesh_with_field(comm, recvFieldName, recvFieldValue,
                                         stk::topology::NODE_RANK, Ioss::Field::MESH);

  stk::transfer_util::TransferMainSettings settings;

  stk::transfer_util::TransferMainBroker broker(comm, sendBulk, recvBulk, settings);

  EXPECT_ANY_THROW(broker.check_and_create_fields());

}

TEST(TransferMainBroker, recvFieldNotTransientDifferentSendandReceive)
{
  stk::ParallelMachine comm = stk::parallel_machine_world();
  if(stk::parallel_machine_size(comm) != 1) { GTEST_SKIP(); }

  std::string sendFieldName = "send_field_1";
  double sendFieldValue = 1.234;
  auto sendBulk = stk::unit_test_util::create_mesh_with_field(comm, sendFieldName, sendFieldValue);

  std::string recvFieldName = "recv_field_1";
  double recvFieldValue = 4.321;
  auto recvBulk = stk::unit_test_util::create_mesh_with_field(comm, recvFieldName, recvFieldValue,
                                         stk::topology::NODE_RANK, Ioss::Field::MESH);

  stk::transfer_util::TransferMainSettings settings;
  std::vector<std::pair<std::string, std::string>> fieldList = {{sendFieldName, recvFieldName}};
  settings.set_transfer_field(fieldList[0]);

  stk::transfer_util::TransferMainBroker broker(comm, sendBulk, recvBulk, settings);

  EXPECT_ANY_THROW(broker.check_and_create_fields());

}

TEST(TransferMainBroker, recvFieldNotTransientDifferentSendandReceiveNoFieldList)
{
  stk::ParallelMachine comm = stk::parallel_machine_world();
  if(stk::parallel_machine_size(comm) != 1) { GTEST_SKIP(); }

  std::string sendFieldName = "send_field_1";
  double sendFieldValue = 1.234;
  auto sendBulk = stk::unit_test_util::create_mesh_with_field(comm, sendFieldName, sendFieldValue);

  std::string recvFieldName = "recv_field_1";
  double recvFieldValue = 4.321;
  auto recvBulk = stk::unit_test_util::create_mesh_with_field(comm, recvFieldName, recvFieldValue,
                                         stk::topology::NODE_RANK, Ioss::Field::MESH);

  stk::transfer_util::TransferMainSettings settings;

  stk::transfer_util::TransferMainBroker broker(comm, sendBulk, recvBulk, settings);

  broker.check_parts();
  broker.check_and_create_fields();

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
  auto sendBulk = stk::unit_test_util::create_mesh_with_field(comm, sendFieldName, sendFieldValue,
                                         stk::topology::NODE_RANK, Ioss::Field::TRANSIENT,
                                         "generated:4x1x1|tets");

  std::string recvFieldName = "field_1";
  double recvFieldValue = 4.321;
  auto recvBulk = stk::unit_test_util::create_mesh_with_field(comm, recvFieldName, recvFieldValue);

  stk::transfer_util::TransferMainSettings settings;

  stk::transfer_util::TransferMainBroker broker(comm, sendBulk, recvBulk, settings);

  broker.check_parts();
  broker.check_and_create_fields();

  check_field_values(sendBulk, sendFieldName, sendFieldValue);
  check_field_values(recvBulk, recvFieldName, recvFieldValue);

  broker.transfer();
  
  check_field_values(sendBulk, sendFieldName, sendFieldValue);
  check_field_values(recvBulk, sendFieldName, sendFieldValue);

}

TEST(TransferMainBroker, recvMeshTets)
{
  stk::ParallelMachine comm = stk::parallel_machine_world();
  if(stk::parallel_machine_size(comm) != 1) { GTEST_SKIP(); }

  std::string sendFieldName = "send_field_1";
  double sendFieldValue = 1.234;
  auto sendBulk = stk::unit_test_util::create_mesh_with_field(comm, sendFieldName, sendFieldValue);

  std::string recvFieldName = "recv_field_1";
  double recvFieldValue = 4.321;
  auto recvBulk = stk::unit_test_util::create_mesh_with_field(comm, recvFieldName, recvFieldValue,
                                         stk::topology::NODE_RANK, Ioss::Field::TRANSIENT,
                                         "generated:4x1x1|tets");

  stk::transfer_util::TransferMainSettings settings;

  stk::transfer_util::TransferMainBroker broker(comm, sendBulk, recvBulk, settings);

  broker.check_parts();
  broker.check_and_create_fields();

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

  broker.check_parts();
  broker.check_and_create_fields();

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

  broker.check_parts();
  broker.check_and_create_fields();

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

  broker.check_parts();
  broker.check_and_create_fields();

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
