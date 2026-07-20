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
#include "stk_search_util/MasterElementProviderIntrepid2.hpp"
#include "stk_unit_test_utils/ioUtils.hpp"
#include "stk_unit_test_utils/meshCreationHelpers.hpp"
#include "stk_unit_test_utils/CommandLineArgs.hpp"
#include "stk_unit_test_utils/BuildMesh.hpp"
#include "stk_mesh/base/ForEachEntity.hpp"
#include "stk_unit_test_utils/FieldEvaluator.hpp"
#include "stk_unit_test_utils/TransferTesters.hpp"
#include <memory>

namespace {

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~Send field must ALWAYS be node rank~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

class TestInterpolationSendFieldRank : public ::testing::Test
{
public:
  TestInterpolationSendFieldRank()
  : comm(stk::parallel_machine_world()),
    settings(),
    broker()
  {}

  void setup_interp_transfer_with_send_field_rank(stk::mesh::EntityRank sendFieldRank)
  {
    auto sendBulk = stk::unit_test_util::create_mesh_with_field(comm, sendFieldRank);
    auto recvBulk = stk::unit_test_util::create_mesh_with_field(comm);
  
    broker = std::make_shared<stk::transfer_util::TransferMainBroker>(comm, sendBulk, recvBulk, settings);
  }

  stk::ParallelMachine comm;
  stk::transfer_util::TransferMainSettings settings;
  std::shared_ptr<stk::transfer_util::TransferMainBroker> broker;
};

TEST_F(TestInterpolationSendFieldRank, nodeRank)
{
  if(stk::parallel_machine_size(comm) != 1) { GTEST_SKIP(); }

  setup_interp_transfer_with_send_field_rank(stk::topology::NODE_RANK);

  EXPECT_NO_THROW(broker->check_and_create_fields_impl());
}

TEST_F(TestInterpolationSendFieldRank, elemRank)
{
  if(stk::parallel_machine_size(comm) != 1) { GTEST_SKIP(); }

  setup_interp_transfer_with_send_field_rank(stk::topology::ELEM_RANK);

  EXPECT_ANY_THROW(broker->check_and_create_fields_impl());
}

TEST_F(TestInterpolationSendFieldRank, faceRank)
{
  if(stk::parallel_machine_size(comm) != 1) { GTEST_SKIP(); }

  setup_interp_transfer_with_send_field_rank(stk::topology::FACE_RANK);

  EXPECT_ANY_THROW(broker->check_and_create_fields_impl());
}

TEST_F(TestInterpolationSendFieldRank, edgeRank)
{
  if(stk::parallel_machine_size(comm) != 1) { GTEST_SKIP(); }

  setup_interp_transfer_with_send_field_rank(stk::topology::EDGE_RANK);

  EXPECT_ANY_THROW(broker->check_and_create_fields_impl());
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//~~~~~~~~~~~~~~~~~~~~~~~Send parts must be elem rank due to Intrepid2 limitation~~~~~~~~~~~~~~~~~~~~~~
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

class TestInterpolationSendParts : public ::testing::Test
{
public:
  TestInterpolationSendParts()
  : comm(stk::parallel_machine_world()),
    settings(),
    broker()
  {}

  void setup_interp_transfer_with_send_recv_parts(const std::vector<std::string>& sendPartNames,
                                                  const std::vector<std::string>& recvPartNames)

  {
    settings.set_transfer_send_parts(sendPartNames);
    settings.set_transfer_recv_parts(recvPartNames);

    auto sendBulk = stk::unit_test_util::create_mesh_with_field(comm);
    auto recvBulk = stk::unit_test_util::create_mesh_with_field(comm);
    broker = std::make_shared<stk::transfer_util::TransferMainBroker>(comm, sendBulk, recvBulk, settings);
  }

  stk::ParallelMachine comm;
  stk::transfer_util::TransferMainSettings settings;
  std::shared_ptr<stk::transfer_util::TransferMainBroker> broker;
};

TEST_F(TestInterpolationSendParts, sendPartNodeRank)
{
  if(stk::parallel_machine_size(stk::parallel_machine_world()) != 1) { GTEST_SKIP(); }

  setup_interp_transfer_with_send_recv_parts({"nodeset_1"}, {"block_1"});

  EXPECT_NO_THROW(broker->check_part_names());
  EXPECT_NO_THROW(broker->check_part_consistency());
  EXPECT_ANY_THROW(broker->check_part_ranks_for_interpolation());
}

TEST_F(TestInterpolationSendParts, sendPartElemRank)
{
  if(stk::parallel_machine_size(stk::parallel_machine_world()) != 1) { GTEST_SKIP(); }

  setup_interp_transfer_with_send_recv_parts({"block_1"}, {"block_1"});

  EXPECT_NO_THROW(broker->check_part_names());
  EXPECT_NO_THROW(broker->check_part_consistency());
  EXPECT_NO_THROW(broker->check_part_ranks_for_interpolation());
}

TEST_F(TestInterpolationSendParts, sendPartFaceRank)
{
  if(stk::parallel_machine_size(stk::parallel_machine_world()) != 1) { GTEST_SKIP(); }

  setup_interp_transfer_with_send_recv_parts({"surface_1"}, {"block_1"});

  EXPECT_NO_THROW(broker->check_part_names());
  EXPECT_NO_THROW(broker->check_part_consistency());
  EXPECT_ANY_THROW(broker->check_part_ranks_for_interpolation());
}

TEST_F(TestInterpolationSendParts, sendPartEdgeRank)
{
  if(stk::parallel_machine_size(stk::parallel_machine_world()) != 1) { GTEST_SKIP(); }

  setup_interp_transfer_with_send_recv_parts({"edges_1"}, {"block_1"});

  EXPECT_NO_THROW(broker->check_part_names());
  EXPECT_NO_THROW(broker->check_part_consistency());
  EXPECT_ANY_THROW(broker->check_part_ranks_for_interpolation());
}

TEST_F(TestInterpolationSendParts, sendPartsDifferentRanks)
{
  if(stk::parallel_machine_size(stk::parallel_machine_world()) != 1) { GTEST_SKIP(); }

  setup_interp_transfer_with_send_recv_parts({"block_1", "nodeset_1"}, {"block_1"});

  EXPECT_NO_THROW(broker->check_part_names());
  EXPECT_ANY_THROW(broker->check_part_consistency());
  EXPECT_ANY_THROW(broker->check_part_ranks_for_interpolation());
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//~~~~~~~~~~~~~~~~Recv field rank can be node/edge/face/elem depending on recv type~~~~~~~~~~~~~~~~~~~~
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

TEST(TransferMainInterpolation, recvFieldNodeRankNodeRecvTypeFieldList)
{
  stk::ParallelMachine comm = stk::parallel_machine_world();
  if(stk::parallel_machine_size(comm) != 1) { GTEST_SKIP(); }

  auto sendBulk = stk::unit_test_util::create_mesh_with_field(comm);
  auto recvBulk = stk::unit_test_util::create_mesh_with_field(comm, stk::topology::NODE_RANK);

  stk::transfer_util::TransferMainSettings settings;
  std::vector<std::pair<std::string, std::string>> fieldList = {{"field_1", "field_1"}};
  settings.set_transfer_field(fieldList[0]);
  settings.set_recv_type("NODE");

  stk::transfer_util::TransferMainBroker broker(comm, sendBulk, recvBulk, settings);

  EXPECT_NO_THROW(broker.check_and_create_fields_impl());
  EXPECT_NO_THROW(broker.check_and_create_fields());
}

TEST(TransferMainInterpolation, recvFieldNodeRankElemRecvTypeFieldList)
{
  stk::ParallelMachine comm = stk::parallel_machine_world();
  if(stk::parallel_machine_size(comm) != 1) { GTEST_SKIP(); }

  auto sendBulk = stk::unit_test_util::create_mesh_with_field(comm);
  auto recvBulk = stk::unit_test_util::create_mesh_with_field(comm, stk::topology::NODE_RANK);

  stk::transfer_util::TransferMainSettings settings;
  std::vector<std::pair<std::string, std::string>> fieldList = {{"field_1", "field_1"}};
  settings.set_transfer_field(fieldList[0]);
  settings.set_recv_type("ELEMENT_GAUSS_POINT");

  stk::transfer_util::TransferMainBroker broker(comm, sendBulk, recvBulk, settings);

  EXPECT_ANY_THROW(broker.check_and_create_fields_impl());
  EXPECT_ANY_THROW(broker.check_and_create_fields());
}

TEST(TransferMainInterpolation, recvFieldNodeRankFaceRecvTypeFieldList)
{
  stk::ParallelMachine comm = stk::parallel_machine_world();
  if(stk::parallel_machine_size(comm) != 1) { GTEST_SKIP(); }

  auto sendBulk = stk::unit_test_util::create_mesh_with_field(comm);
  auto recvBulk = stk::unit_test_util::create_mesh_with_field(comm, stk::topology::NODE_RANK);

  stk::transfer_util::TransferMainSettings settings;
  std::vector<std::pair<std::string, std::string>> fieldList = {{"field_1", "field_1"}};
  settings.set_transfer_field(fieldList[0]);
  settings.set_recv_type("FACE_CENTROID");

  stk::transfer_util::TransferMainBroker broker(comm, sendBulk, recvBulk, settings);

  EXPECT_ANY_THROW(broker.check_and_create_fields_impl());
  EXPECT_ANY_THROW(broker.check_and_create_fields());
}

TEST(TransferMainInterpolation, recvFieldNodeRankEdgeRecvTypeFieldList)
{
  stk::ParallelMachine comm = stk::parallel_machine_world();
  if(stk::parallel_machine_size(comm) != 1) { GTEST_SKIP(); }

  auto sendBulk = stk::unit_test_util::create_mesh_with_field(comm);
  auto recvBulk = stk::unit_test_util::create_mesh_with_field(comm, stk::topology::NODE_RANK);

  stk::transfer_util::TransferMainSettings settings;
  std::vector<std::pair<std::string, std::string>> fieldList = {{"field_1", "field_1"}};
  settings.set_transfer_field(fieldList[0]);
  settings.set_recv_type("EDGE_CENTROID");

  stk::transfer_util::TransferMainBroker broker(comm, sendBulk, recvBulk, settings);

  EXPECT_ANY_THROW(broker.check_and_create_fields_impl());
  EXPECT_ANY_THROW(broker.check_and_create_fields());
}

TEST(TransferMainInterpolation, recvFieldNodeRankNodeRecvTypeNoFieldList)
{
  stk::ParallelMachine comm = stk::parallel_machine_world();
  if(stk::parallel_machine_size(comm) != 1) { GTEST_SKIP(); }

  auto sendBulk = stk::unit_test_util::create_mesh_with_field(comm);
  auto recvBulk = stk::unit_test_util::create_mesh_with_field(comm, stk::topology::NODE_RANK);

  stk::transfer_util::TransferMainSettings settings;
  settings.set_recv_type("NODE");

  stk::transfer_util::TransferMainBroker broker(comm, sendBulk, recvBulk, settings);

  EXPECT_NO_THROW(broker.check_and_create_fields_impl());
  EXPECT_NO_THROW(broker.check_and_create_fields());
}

TEST(TransferMainInterpolation, recvFieldElemRankNodeRecvType)
{
  stk::ParallelMachine comm = stk::parallel_machine_world();
  if(stk::parallel_machine_size(comm) != 1) { GTEST_SKIP(); }

  auto sendBulk = stk::unit_test_util::create_mesh_with_field(comm);
  auto recvBulk = stk::unit_test_util::create_mesh_with_field(comm, stk::topology::ELEM_RANK);

  stk::transfer_util::TransferMainSettings settings;
  settings.set_recv_type("NODE");

  stk::transfer_util::TransferMainBroker broker(comm, sendBulk, recvBulk, settings);

  EXPECT_ANY_THROW(broker.check_and_create_fields_impl());
  EXPECT_ANY_THROW(broker.check_and_create_fields());
}

TEST(TransferMainInterpolation, recvFieldElemRankElemRecvType)
{
  stk::ParallelMachine comm = stk::parallel_machine_world();
  if(stk::parallel_machine_size(comm) != 1) { GTEST_SKIP(); }

  auto sendBulk = stk::unit_test_util::create_mesh_with_field(comm);
  auto recvBulk = stk::unit_test_util::create_mesh_with_field(comm, stk::topology::ELEM_RANK);

  stk::transfer_util::TransferMainSettings settings;
  settings.set_recv_type("ELEMENT_CENTROID");

  stk::transfer_util::TransferMainBroker broker(comm, sendBulk, recvBulk, settings);

  EXPECT_NO_THROW(broker.check_and_create_fields_impl());
  EXPECT_NO_THROW(broker.check_and_create_fields());
}

TEST(TransferMainInterpolation, recvFieldElemRankFaceRecvType)
{
  stk::ParallelMachine comm = stk::parallel_machine_world();
  if(stk::parallel_machine_size(comm) != 1) { GTEST_SKIP(); }

  auto sendBulk = stk::unit_test_util::create_mesh_with_field(comm);
  auto recvBulk = stk::unit_test_util::create_mesh_with_field(comm, stk::topology::ELEM_RANK);

  stk::transfer_util::TransferMainSettings settings;
  settings.set_recv_type("FACE_CENTROID");

  stk::transfer_util::TransferMainBroker broker(comm, sendBulk, recvBulk, settings);

  EXPECT_ANY_THROW(broker.check_and_create_fields_impl());
  EXPECT_ANY_THROW(broker.check_and_create_fields());
}

TEST(TransferMainInterpolation, recvFieldElemRankEdgeRecvType)
{
  stk::ParallelMachine comm = stk::parallel_machine_world();
  if(stk::parallel_machine_size(comm) != 1) { GTEST_SKIP(); }

  auto sendBulk = stk::unit_test_util::create_mesh_with_field(comm);
  auto recvBulk = stk::unit_test_util::create_mesh_with_field(comm, stk::topology::ELEM_RANK);

  stk::transfer_util::TransferMainSettings settings;
  settings.set_recv_type("EDGE_GAUSS_POINT");

  stk::transfer_util::TransferMainBroker broker(comm, sendBulk, recvBulk, settings);

  EXPECT_ANY_THROW(broker.check_and_create_fields_impl());
  EXPECT_ANY_THROW(broker.check_and_create_fields());
}

TEST(TransferMainInterpolation, recvFieldFaceRankNodeRecvType)
{
  stk::ParallelMachine comm = stk::parallel_machine_world();
  if(stk::parallel_machine_size(comm) != 1) { GTEST_SKIP(); }

  auto sendBulk = stk::unit_test_util::create_mesh_with_field(comm);
  auto recvBulk = stk::unit_test_util::create_mesh_with_field(comm, stk::topology::FACE_RANK);

  stk::transfer_util::TransferMainSettings settings;
  settings.set_recv_type("NODE");

  stk::transfer_util::TransferMainBroker broker(comm, sendBulk, recvBulk, settings);

  EXPECT_ANY_THROW(broker.check_and_create_fields_impl());
  EXPECT_ANY_THROW(broker.check_and_create_fields());
}

TEST(TransferMainInterpolation, recvFieldFaceRankElemRecvType)
{
  stk::ParallelMachine comm = stk::parallel_machine_world();
  if(stk::parallel_machine_size(comm) != 1) { GTEST_SKIP(); }

  auto sendBulk = stk::unit_test_util::create_mesh_with_field(comm);
  auto recvBulk = stk::unit_test_util::create_mesh_with_field(comm, stk::topology::FACE_RANK);

  stk::transfer_util::TransferMainSettings settings;
  settings.set_recv_type("ELEMENT_GAUSS_POINT");

  stk::transfer_util::TransferMainBroker broker(comm, sendBulk, recvBulk, settings);

  EXPECT_ANY_THROW(broker.check_and_create_fields_impl());
  EXPECT_ANY_THROW(broker.check_and_create_fields());
}

TEST(TransferMainInterpolation, recvFieldFaceRankFaceRecvType)
{
  stk::ParallelMachine comm = stk::parallel_machine_world();
  if(stk::parallel_machine_size(comm) != 1) { GTEST_SKIP(); }

  auto sendBulk = stk::unit_test_util::create_mesh_with_field(comm);
  auto recvBulk = stk::unit_test_util::create_mesh_with_field(comm, stk::topology::FACE_RANK);

  stk::transfer_util::TransferMainSettings settings;
  settings.set_recv_type("FACE_GAUSS_POINT");

  stk::transfer_util::TransferMainBroker broker(comm, sendBulk, recvBulk, settings);

  EXPECT_NO_THROW(broker.check_and_create_fields_impl());
  EXPECT_NO_THROW(broker.check_and_create_fields());
}

TEST(TransferMainInterpolation, recvFieldFaceRankEdgeRecvType)
{
  stk::ParallelMachine comm = stk::parallel_machine_world();
  if(stk::parallel_machine_size(comm) != 1) { GTEST_SKIP(); }

  auto sendBulk = stk::unit_test_util::create_mesh_with_field(comm);
  auto recvBulk = stk::unit_test_util::create_mesh_with_field(comm, stk::topology::FACE_RANK);

  stk::transfer_util::TransferMainSettings settings;
  settings.set_recv_type("EDGE_CENTROID");

  stk::transfer_util::TransferMainBroker broker(comm, sendBulk, recvBulk, settings);

  EXPECT_ANY_THROW(broker.check_and_create_fields_impl());
  EXPECT_ANY_THROW(broker.check_and_create_fields());
}

TEST(TransferMainInterpolation, recvFieldEdgeRankNodeRecvType)
{
  stk::ParallelMachine comm = stk::parallel_machine_world();
  if(stk::parallel_machine_size(comm) != 1) { GTEST_SKIP(); }

  auto sendBulk = stk::unit_test_util::create_mesh_with_field(comm);
  auto recvBulk = stk::unit_test_util::create_mesh_with_field(comm, stk::topology::EDGE_RANK);

  stk::transfer_util::TransferMainSettings settings;
  settings.set_recv_type("NODE");

  stk::transfer_util::TransferMainBroker broker(comm, sendBulk, recvBulk, settings);

  EXPECT_ANY_THROW(broker.check_and_create_fields_impl());
  EXPECT_ANY_THROW(broker.check_and_create_fields());
}

TEST(TransferMainInterpolation, recvFieldEdgeRankElemRecvType)
{
  stk::ParallelMachine comm = stk::parallel_machine_world();
  if(stk::parallel_machine_size(comm) != 1) { GTEST_SKIP(); }

  auto sendBulk = stk::unit_test_util::create_mesh_with_field(comm);
  auto recvBulk = stk::unit_test_util::create_mesh_with_field(comm, stk::topology::EDGE_RANK);

  stk::transfer_util::TransferMainSettings settings;
  settings.set_recv_type("ELEMENT_CENTROID");

  stk::transfer_util::TransferMainBroker broker(comm, sendBulk, recvBulk, settings);

  EXPECT_ANY_THROW(broker.check_and_create_fields_impl());
  EXPECT_ANY_THROW(broker.check_and_create_fields());
}
TEST(TransferMainInterpolation, recvFieldEdgeRankFaceRecvType)
{
  stk::ParallelMachine comm = stk::parallel_machine_world();
  if(stk::parallel_machine_size(comm) != 1) { GTEST_SKIP(); }

  auto sendBulk = stk::unit_test_util::create_mesh_with_field(comm);
  auto recvBulk = stk::unit_test_util::create_mesh_with_field(comm, stk::topology::EDGE_RANK);

  stk::transfer_util::TransferMainSettings settings;
  settings.set_recv_type("FACE_GAUSS_POINT");

  stk::transfer_util::TransferMainBroker broker(comm, sendBulk, recvBulk, settings);

  EXPECT_ANY_THROW(broker.check_and_create_fields_impl());
  EXPECT_ANY_THROW(broker.check_and_create_fields());
}

TEST(TransferMainInterpolation, recvFieldEdgeRankEdgeRecvType)
{
  stk::ParallelMachine comm = stk::parallel_machine_world();
  if(stk::parallel_machine_size(comm) != 1) { GTEST_SKIP(); }

  auto sendBulk = stk::unit_test_util::create_mesh_with_field(comm);
  auto recvBulk = stk::unit_test_util::create_mesh_with_field(comm, stk::topology::EDGE_RANK);

  stk::transfer_util::TransferMainSettings settings;
  settings.set_recv_type("EDGE_GAUSS_POINT");

  stk::transfer_util::TransferMainBroker broker(comm, sendBulk, recvBulk, settings);

  EXPECT_NO_THROW(broker.check_and_create_fields_impl());
  EXPECT_NO_THROW(broker.check_and_create_fields());
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//~~~~~~~~~~~~~~~~~~~~Recv part rank can be edge/face/elem depending on recv type~~~~~~~~~~~~~~~~~~~~~~
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

TEST(TransferMainInterpolation, recvPartNodeRank)
{
  stk::ParallelMachine comm = stk::parallel_machine_world();
  if(stk::parallel_machine_size(comm) != 1) { GTEST_SKIP(); }

  auto sendBulk = stk::unit_test_util::create_mesh_with_field(comm);
  auto recvBulk = stk::unit_test_util::create_mesh_with_field(comm);

  stk::transfer_util::TransferMainSettings settings;
  settings.set_transfer_send_parts({"block_1"});
  settings.set_transfer_recv_parts({"nodeset_1"});

  stk::transfer_util::TransferMainBroker broker(comm, sendBulk, recvBulk, settings);

  auto sendPartRank = sendBulk->mesh_meta_data().get_part("block_1")->primary_entity_rank();
  EXPECT_EQ(sendPartRank, stk::topology::ELEM_RANK);

  auto recvPartRank = recvBulk->mesh_meta_data().get_part("nodeset_1")->primary_entity_rank();
  EXPECT_EQ(recvPartRank, stk::topology::NODE_RANK);

  EXPECT_NO_THROW(broker.check_part_names());
  EXPECT_NO_THROW(broker.check_part_consistency());
  EXPECT_ANY_THROW(broker.check_part_ranks_for_interpolation());
  EXPECT_ANY_THROW(broker.check_parts());
}

TEST(TransferMainInterpolation, recvPartElemRankElemRecvType)
{
  stk::ParallelMachine comm = stk::parallel_machine_world();
  if(stk::parallel_machine_size(comm) != 1) { GTEST_SKIP(); }

  auto sendBulk = stk::unit_test_util::create_mesh_with_field(comm);
  auto recvBulk = stk::unit_test_util::create_mesh_with_field(comm);

  stk::transfer_util::TransferMainSettings settings;
  settings.set_recv_type("ELEMENT_CENTROID");
  settings.set_transfer_send_parts({"block_1"});
  settings.set_transfer_recv_parts({"block_1"});

  stk::transfer_util::TransferMainBroker broker(comm, sendBulk, recvBulk, settings);

  auto sendPartRank = sendBulk->mesh_meta_data().get_part("block_1")->primary_entity_rank();
  EXPECT_EQ(sendPartRank, stk::topology::ELEM_RANK);

  auto recvPartRank = recvBulk->mesh_meta_data().get_part("block_1")->primary_entity_rank();
  EXPECT_EQ(recvPartRank, stk::topology::ELEM_RANK);

  EXPECT_NO_THROW(broker.check_part_names());
  EXPECT_NO_THROW(broker.check_part_consistency());
  EXPECT_NO_THROW(broker.check_part_ranks_for_interpolation());
  EXPECT_NO_THROW(broker.check_parts());
}

TEST(TransferMainInterpolation, recvPartElemRankNodeRecvType)
{
  stk::ParallelMachine comm = stk::parallel_machine_world();
  if(stk::parallel_machine_size(comm) != 1) { GTEST_SKIP(); }

  auto sendBulk = stk::unit_test_util::create_mesh_with_field(comm);
  auto recvBulk = stk::unit_test_util::create_mesh_with_field(comm);

  stk::transfer_util::TransferMainSettings settings;
  settings.set_recv_type("NODE");
  settings.set_transfer_send_parts({"block_1"});
  settings.set_transfer_recv_parts({"block_1"});

  stk::transfer_util::TransferMainBroker broker(comm, sendBulk, recvBulk, settings);

  auto sendPartRank = sendBulk->mesh_meta_data().get_part("block_1")->primary_entity_rank();
  EXPECT_EQ(sendPartRank, stk::topology::ELEM_RANK);

  auto recvPartRank = recvBulk->mesh_meta_data().get_part("block_1")->primary_entity_rank();
  EXPECT_EQ(recvPartRank, stk::topology::ELEM_RANK);

  EXPECT_NO_THROW(broker.check_part_names());
  EXPECT_NO_THROW(broker.check_part_consistency());
  EXPECT_NO_THROW(broker.check_part_ranks_for_interpolation());
  EXPECT_NO_THROW(broker.check_parts());
}

TEST(TransferMainInterpolation, recvPartElemRankFaceRecvType)
{
  stk::ParallelMachine comm = stk::parallel_machine_world();
  if(stk::parallel_machine_size(comm) != 1) { GTEST_SKIP(); }

  auto sendBulk = stk::unit_test_util::create_mesh_with_field(comm);
  auto recvBulk = stk::unit_test_util::create_mesh_with_field(comm);

  stk::transfer_util::TransferMainSettings settings;
  settings.set_recv_type("FACE_GAUSS_POINT");
  settings.set_transfer_send_parts({"block_1"});
  settings.set_transfer_recv_parts({"block_1"});

  stk::transfer_util::TransferMainBroker broker(comm, sendBulk, recvBulk, settings);

  auto sendPartRank = sendBulk->mesh_meta_data().get_part("block_1")->primary_entity_rank();
  EXPECT_EQ(sendPartRank, stk::topology::ELEM_RANK);

  auto recvPartRank = recvBulk->mesh_meta_data().get_part("block_1")->primary_entity_rank();
  EXPECT_EQ(recvPartRank, stk::topology::ELEM_RANK);

  EXPECT_NO_THROW(broker.check_part_names());
  EXPECT_NO_THROW(broker.check_part_consistency());
  EXPECT_ANY_THROW(broker.check_part_ranks_for_interpolation());
  EXPECT_ANY_THROW(broker.check_parts());
}

TEST(TransferMainInterpolation, recvPartFaceRankElemRecvType)
{
  stk::ParallelMachine comm = stk::parallel_machine_world();
  if(stk::parallel_machine_size(comm) != 1) { GTEST_SKIP(); }

  auto sendBulk = stk::unit_test_util::create_mesh_with_field(comm);
  auto recvBulk = stk::unit_test_util::create_mesh_with_field(comm);

  stk::transfer_util::TransferMainSettings settings;
  settings.set_recv_type("ELEMENT_GAUSS_POINT");
  settings.set_transfer_send_parts({"block_1"});
  settings.set_transfer_recv_parts({"surface_1"});

  stk::transfer_util::TransferMainBroker broker(comm, sendBulk, recvBulk, settings);

  auto sendPartRank = sendBulk->mesh_meta_data().get_part("block_1")->primary_entity_rank();
  EXPECT_EQ(sendPartRank, stk::topology::ELEM_RANK);

  auto recvPartRank = recvBulk->mesh_meta_data().get_part("surface_1")->primary_entity_rank();
  EXPECT_EQ(recvPartRank, stk::topology::FACE_RANK);

  EXPECT_NO_THROW(broker.check_part_names());
  EXPECT_NO_THROW(broker.check_part_consistency());
  EXPECT_ANY_THROW(broker.check_part_ranks_for_interpolation());
  EXPECT_ANY_THROW(broker.check_parts());
}

TEST(TransferMainInterpolation, recvPartFaceRankFaceRecvType)
{
  stk::ParallelMachine comm = stk::parallel_machine_world();
  if(stk::parallel_machine_size(comm) != 1) { GTEST_SKIP(); }

  auto sendBulk = stk::unit_test_util::create_mesh_with_field(comm);
  auto recvBulk = stk::unit_test_util::create_mesh_with_field(comm);

  stk::transfer_util::TransferMainSettings settings;
  settings.set_recv_type("FACE_CENTROID");
  settings.set_transfer_send_parts({"block_1"});
  settings.set_transfer_recv_parts({"surface_1"});

  stk::transfer_util::TransferMainBroker broker(comm, sendBulk, recvBulk, settings);

  auto sendPartRank = sendBulk->mesh_meta_data().get_part("block_1")->primary_entity_rank();
  EXPECT_EQ(sendPartRank, stk::topology::ELEM_RANK);

  auto recvPartRank = recvBulk->mesh_meta_data().get_part("surface_1")->primary_entity_rank();
  EXPECT_EQ(recvPartRank, stk::topology::FACE_RANK);

  EXPECT_NO_THROW(broker.check_part_names());
  EXPECT_NO_THROW(broker.check_part_consistency());
  EXPECT_NO_THROW(broker.check_part_ranks_for_interpolation());
  EXPECT_NO_THROW(broker.check_parts());
}

TEST(TransferMainInterpolation, recvPartFaceRankNodeRecvType)
{
  stk::ParallelMachine comm = stk::parallel_machine_world();
  if(stk::parallel_machine_size(comm) != 1) { GTEST_SKIP(); }

  auto sendBulk = stk::unit_test_util::create_mesh_with_field(comm);
  auto recvBulk = stk::unit_test_util::create_mesh_with_field(comm);

  stk::transfer_util::TransferMainSettings settings;
  settings.set_recv_type("NODE");
  settings.set_transfer_send_parts({"block_1"});
  settings.set_transfer_recv_parts({"surface_1"});

  stk::transfer_util::TransferMainBroker broker(comm, sendBulk, recvBulk, settings);

  auto sendPartRank = sendBulk->mesh_meta_data().get_part("block_1")->primary_entity_rank();
  EXPECT_EQ(sendPartRank, stk::topology::ELEM_RANK);

  auto recvPartRank = recvBulk->mesh_meta_data().get_part("surface_1")->primary_entity_rank();
  EXPECT_EQ(recvPartRank, stk::topology::FACE_RANK);

  EXPECT_NO_THROW(broker.check_part_names());
  EXPECT_NO_THROW(broker.check_part_consistency());
  EXPECT_NO_THROW(broker.check_part_ranks_for_interpolation());
  EXPECT_NO_THROW(broker.check_parts());
}

TEST(TransferMainInterpolation, recvPartEdgeRankElemRecvType)
{
  stk::ParallelMachine comm = stk::parallel_machine_world();
  if(stk::parallel_machine_size(comm) != 1) { GTEST_SKIP(); }

  auto sendBulk = stk::unit_test_util::create_mesh_with_field(comm);
  auto recvBulk = stk::unit_test_util::create_mesh_with_field(comm);

  stk::transfer_util::TransferMainSettings settings;
  settings.set_recv_type("ELEMENT_GAUSS_POINT");
  settings.set_transfer_send_parts({"block_1"});
  settings.set_transfer_recv_parts({"edges_1"});

  stk::transfer_util::TransferMainBroker broker(comm, sendBulk, recvBulk, settings);

  auto sendPartRank = sendBulk->mesh_meta_data().get_part("block_1")->primary_entity_rank();
  EXPECT_EQ(sendPartRank, stk::topology::ELEM_RANK);

  auto recvPartRank = recvBulk->mesh_meta_data().get_part("edges_1")->primary_entity_rank();
  EXPECT_EQ(recvPartRank, stk::topology::EDGE_RANK);

  EXPECT_NO_THROW(broker.check_part_names());
  EXPECT_NO_THROW(broker.check_part_consistency());
  EXPECT_ANY_THROW(broker.check_part_ranks_for_interpolation());
  EXPECT_ANY_THROW(broker.check_parts());
}

TEST(TransferMainInterpolation, recvPartEdgeRankFaceRecvType)
{
  stk::ParallelMachine comm = stk::parallel_machine_world();
  if(stk::parallel_machine_size(comm) != 1) { GTEST_SKIP(); }

  auto sendBulk = stk::unit_test_util::create_mesh_with_field(comm);
  auto recvBulk = stk::unit_test_util::create_mesh_with_field(comm);

  stk::transfer_util::TransferMainSettings settings;
  settings.set_recv_type("FACE_CENTROID");
  settings.set_transfer_send_parts({"block_1"});
  settings.set_transfer_recv_parts({"edges_1"});

  stk::transfer_util::TransferMainBroker broker(comm, sendBulk, recvBulk, settings);

  auto sendPartRank = sendBulk->mesh_meta_data().get_part("block_1")->primary_entity_rank();
  EXPECT_EQ(sendPartRank, stk::topology::ELEM_RANK);

  auto recvPartRank = recvBulk->mesh_meta_data().get_part("edges_1")->primary_entity_rank();
  EXPECT_EQ(recvPartRank, stk::topology::EDGE_RANK);

  EXPECT_NO_THROW(broker.check_part_names());
  EXPECT_NO_THROW(broker.check_part_consistency());
  EXPECT_ANY_THROW(broker.check_part_ranks_for_interpolation());
  EXPECT_ANY_THROW(broker.check_parts());
}

TEST(TransferMainInterpolation, recvPartEdgeRankEdgeRecvType)
{
  stk::ParallelMachine comm = stk::parallel_machine_world();
  if(stk::parallel_machine_size(comm) != 1) { GTEST_SKIP(); }

  auto sendBulk = stk::unit_test_util::create_mesh_with_field(comm);
  auto recvBulk = stk::unit_test_util::create_mesh_with_field(comm);

  stk::transfer_util::TransferMainSettings settings;
  settings.set_recv_type("EDGE_CENTROID");
  settings.set_transfer_send_parts({"block_1"});
  settings.set_transfer_recv_parts({"edges_1"});

  stk::transfer_util::TransferMainBroker broker(comm, sendBulk, recvBulk, settings);

  auto sendPartRank = sendBulk->mesh_meta_data().get_part("block_1")->primary_entity_rank();
  EXPECT_EQ(sendPartRank, stk::topology::ELEM_RANK);

  auto recvPartRank = recvBulk->mesh_meta_data().get_part("edges_1")->primary_entity_rank();
  EXPECT_EQ(recvPartRank, stk::topology::EDGE_RANK);

  EXPECT_NO_THROW(broker.check_part_names());
  EXPECT_NO_THROW(broker.check_part_consistency());
  EXPECT_NO_THROW(broker.check_part_ranks_for_interpolation());
  EXPECT_NO_THROW(broker.check_parts());
}

TEST(TransferMainInterpolation, recvPartEdgeRankNodeRecvType)
{
  stk::ParallelMachine comm = stk::parallel_machine_world();
  if(stk::parallel_machine_size(comm) != 1) { GTEST_SKIP(); }

  auto sendBulk = stk::unit_test_util::create_mesh_with_field(comm);
  auto recvBulk = stk::unit_test_util::create_mesh_with_field(comm);

  stk::transfer_util::TransferMainSettings settings;
  settings.set_recv_type("NODE");
  settings.set_transfer_send_parts({"block_1"});
  settings.set_transfer_recv_parts({"edges_1"});

  stk::transfer_util::TransferMainBroker broker(comm, sendBulk, recvBulk, settings);

  auto sendPartRank = sendBulk->mesh_meta_data().get_part("block_1")->primary_entity_rank();
  EXPECT_EQ(sendPartRank, stk::topology::ELEM_RANK);

  auto recvPartRank = recvBulk->mesh_meta_data().get_part("edges_1")->primary_entity_rank();
  EXPECT_EQ(recvPartRank, stk::topology::EDGE_RANK);

  EXPECT_NO_THROW(broker.check_part_names());
  EXPECT_NO_THROW(broker.check_part_consistency());
  EXPECT_NO_THROW(broker.check_part_ranks_for_interpolation());
  EXPECT_NO_THROW(broker.check_parts());
}

TEST(TransferMainInterpolation, recvPartsDifferentRanks)
{
  stk::ParallelMachine comm = stk::parallel_machine_world();
  if(stk::parallel_machine_size(comm) != 1) { GTEST_SKIP(); }

  auto sendBulk = stk::unit_test_util::create_mesh_with_field(comm);
  auto recvBulk = stk::unit_test_util::create_mesh_with_field(comm);

  stk::transfer_util::TransferMainSettings settings;
  settings.set_transfer_send_parts({"block_1"});
  settings.set_transfer_recv_parts({"block_1", "nodeset_1"});

  stk::transfer_util::TransferMainBroker broker(comm, sendBulk, recvBulk, settings);

  auto sendPartRank = sendBulk->mesh_meta_data().get_part("block_1")->primary_entity_rank();
  EXPECT_EQ(sendPartRank, stk::topology::ELEM_RANK);

  auto recvPartRank = recvBulk->mesh_meta_data().get_part("block_1")->primary_entity_rank();
  EXPECT_EQ(recvPartRank, stk::topology::ELEM_RANK);

  auto recvPartRank2 = recvBulk->mesh_meta_data().get_part("nodeset_1")->primary_entity_rank();
  EXPECT_EQ(recvPartRank2, stk::topology::NODE_RANK);

  EXPECT_NO_THROW(broker.check_part_names());
  EXPECT_ANY_THROW(broker.check_part_consistency());
  EXPECT_ANY_THROW(broker.check_part_ranks_for_interpolation());
  EXPECT_ANY_THROW(broker.check_parts());
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//~~~~~~~~~~~~~~~~~~~~~Check that default parts are compatible with recv types~~~~~~~~~~~~~~~~~~~~~~~~~
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

TEST(TransferMainInterpolation, defaultPartsNodeRecvType)
{
  stk::ParallelMachine comm = stk::parallel_machine_world();
  if(stk::parallel_machine_size(comm) != 1) { GTEST_SKIP(); }

  auto sendBulk = stk::unit_test_util::create_mesh_with_field(comm);
  auto recvBulk = stk::unit_test_util::create_mesh_with_field(comm);

  stk::transfer_util::TransferMainSettings settings;
  settings.set_recv_type("NODE");

  stk::transfer_util::TransferMainBroker broker(comm, sendBulk, recvBulk, settings);

  EXPECT_NO_THROW(broker.check_part_names());
  EXPECT_NO_THROW(broker.check_part_consistency());
  EXPECT_NO_THROW(broker.check_part_ranks_for_interpolation());
  EXPECT_NO_THROW(broker.check_parts());
}

TEST(TransferMainInterpolation, defaultPartsElemCentroidRecvType)
{
  stk::ParallelMachine comm = stk::parallel_machine_world();
  if(stk::parallel_machine_size(comm) != 1) { GTEST_SKIP(); }

  auto sendBulk = stk::unit_test_util::create_mesh_with_field(comm);
  auto recvBulk = stk::unit_test_util::create_mesh_with_field(comm);

  stk::transfer_util::TransferMainSettings settings;
  settings.set_recv_type("ELEMENT_CENTROID");

  stk::transfer_util::TransferMainBroker broker(comm, sendBulk, recvBulk, settings);

  EXPECT_NO_THROW(broker.check_part_names());
  EXPECT_NO_THROW(broker.check_part_consistency());
  EXPECT_NO_THROW(broker.check_part_ranks_for_interpolation());
  EXPECT_NO_THROW(broker.check_parts());
}

TEST(TransferMainInterpolation, defaultPartsElemGPRecvType)
{
  stk::ParallelMachine comm = stk::parallel_machine_world();
  if(stk::parallel_machine_size(comm) != 1) { GTEST_SKIP(); }

  auto sendBulk = stk::unit_test_util::create_mesh_with_field(comm);
  auto recvBulk = stk::unit_test_util::create_mesh_with_field(comm);

  stk::transfer_util::TransferMainSettings settings;
  settings.set_recv_type("ELEMENT_GAUSS_POINT");

  stk::transfer_util::TransferMainBroker broker(comm, sendBulk, recvBulk, settings);

  EXPECT_NO_THROW(broker.check_part_names());
  EXPECT_NO_THROW(broker.check_part_consistency());
  EXPECT_NO_THROW(broker.check_part_ranks_for_interpolation());
  EXPECT_NO_THROW(broker.check_parts());
}

TEST(TransferMainInterpolation, defaultPartsFaceCentroidRecvType)
{
  stk::ParallelMachine comm = stk::parallel_machine_world();
  if(stk::parallel_machine_size(comm) != 1) { GTEST_SKIP(); }

  auto sendBulk = stk::unit_test_util::create_mesh_with_field(comm);
  auto recvBulk = stk::unit_test_util::create_mesh_with_field(comm);

  stk::transfer_util::TransferMainSettings settings;
  settings.set_recv_type("FACE_CENTROID");

  stk::transfer_util::TransferMainBroker broker(comm, sendBulk, recvBulk, settings);

  EXPECT_NO_THROW(broker.check_part_names());
  EXPECT_NO_THROW(broker.check_part_consistency());
  EXPECT_NO_THROW(broker.check_part_ranks_for_interpolation());
  EXPECT_NO_THROW(broker.check_parts());
}

TEST(TransferMainInterpolation, defaultPartsFaceGPRecvType)
{
  stk::ParallelMachine comm = stk::parallel_machine_world();
  if(stk::parallel_machine_size(comm) != 1) { GTEST_SKIP(); }

  auto sendBulk = stk::unit_test_util::create_mesh_with_field(comm);
  auto recvBulk = stk::unit_test_util::create_mesh_with_field(comm);

  stk::transfer_util::TransferMainSettings settings;
  settings.set_recv_type("FACE_GAUSS_POINT");

  stk::transfer_util::TransferMainBroker broker(comm, sendBulk, recvBulk, settings);

  EXPECT_NO_THROW(broker.check_part_names());
  EXPECT_NO_THROW(broker.check_part_consistency());
  EXPECT_NO_THROW(broker.check_part_ranks_for_interpolation());
  EXPECT_NO_THROW(broker.check_parts());
}

TEST(TransferMainInterpolation, defaultPartsEdgeCentroidRecvType)
{
  stk::ParallelMachine comm = stk::parallel_machine_world();
  if(stk::parallel_machine_size(comm) != 1) { GTEST_SKIP(); }

  auto sendBulk = stk::unit_test_util::create_mesh_with_field(comm);
  auto recvBulk = stk::unit_test_util::create_mesh_with_field(comm);

  stk::transfer_util::TransferMainSettings settings;
  settings.set_recv_type("EDGE_CENTROID");

  stk::transfer_util::TransferMainBroker broker(comm, sendBulk, recvBulk, settings);

  EXPECT_NO_THROW(broker.check_part_names());
  EXPECT_NO_THROW(broker.check_part_consistency());
  EXPECT_NO_THROW(broker.check_part_ranks_for_interpolation());
  EXPECT_NO_THROW(broker.check_parts());
}

TEST(TransferMainInterpolation, defaultPartsEdgeGPRecvType)
{
  stk::ParallelMachine comm = stk::parallel_machine_world();
  if(stk::parallel_machine_size(comm) != 1) { GTEST_SKIP(); }

  auto sendBulk = stk::unit_test_util::create_mesh_with_field(comm);
  auto recvBulk = stk::unit_test_util::create_mesh_with_field(comm);

  stk::transfer_util::TransferMainSettings settings;
  settings.set_recv_type("EDGE_GAUSS_POINT");

  stk::transfer_util::TransferMainBroker broker(comm, sendBulk, recvBulk, settings);

  EXPECT_NO_THROW(broker.check_part_names());
  EXPECT_NO_THROW(broker.check_part_consistency());
  EXPECT_NO_THROW(broker.check_part_ranks_for_interpolation());
  EXPECT_NO_THROW(broker.check_parts());
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~Correctness testing~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

TEST(TransferMainBroker, elementSendMesh_masterElementInterpolationToNode)
{
  stk::ParallelMachine comm = stk::parallel_machine_world();
  if(stk::parallel_machine_size(comm) != 1) { GTEST_SKIP(); }

  auto sendBulk = stk::unit_test_util::create_mesh_with_field(comm, stk::topology::NODE_RANK);

  auto recvBulk = stk::unit_test_util::create_mesh_with_field(comm, stk::topology::NODE_RANK);

  stk::transfer_util::TransferMainSettings settings;
  settings.set_recv_type("NODE");
  settings.set_transfer_send_parts({"block_1"});
  settings.set_transfer_recv_parts({"block_1"});

  stk::transfer_util::TransferMainBroker broker(comm, sendBulk, recvBulk, settings);

  broker.check_parts();
  broker.check_and_create_fields();

  stk::unit_test_util::LinearFieldEvaluator eval(3u);
  auto sendField = stk::mesh::get_field_by_name("field_1", sendBulk->mesh_meta_data());
  stk::unit_test_util::set_node_field(*sendBulk, *sendField, eval);

  broker.transfer();

  stk::unit_test_util::test_node_transfer(eval, recvBulk, stk::topology::NODE_RANK);
}

TEST(TransferMainBroker, elementSendMesh_masterElementInterpolationToCentroid)
{
  stk::ParallelMachine comm = stk::parallel_machine_world();
  if(stk::parallel_machine_size(comm) != 1) { GTEST_SKIP(); }

  auto sendBulk = stk::unit_test_util::create_mesh_with_field(comm, stk::topology::NODE_RANK);

  auto recvBulk = stk::unit_test_util::create_mesh_with_field(comm, stk::topology::ELEM_RANK);

  stk::transfer_util::TransferMainSettings settings;
  settings.set_recv_type("ELEMENT_CENTROID");
  settings.set_transfer_send_parts({"block_1"});
  settings.set_transfer_recv_parts({"block_1"});

  stk::transfer_util::TransferMainBroker broker(comm, sendBulk, recvBulk, settings);

  broker.check_parts();
  broker.check_and_create_fields();

  stk::unit_test_util::LinearFieldEvaluator eval(3u);
  auto sendField = stk::mesh::get_field_by_name("field_1", sendBulk->mesh_meta_data());
  stk::unit_test_util::set_node_field(*sendBulk, *sendField, eval);

  broker.transfer();

  stk::unit_test_util::test_centroid_transfer(eval, recvBulk, stk::topology::ELEM_RANK);
}

TEST(TransferMainBroker, elementSendMesh_masterElementInterpolationToGaussPoint)
{
  stk::ParallelMachine comm = stk::parallel_machine_world();
  if(stk::parallel_machine_size(comm) != 1) { GTEST_SKIP(); }

  unsigned numSendCopies = 1;
  auto sendBulk = stk::unit_test_util::create_mesh_with_field(comm, stk::topology::NODE_RANK, numSendCopies);
  std::shared_ptr<stk::search::MasterElementProviderInterface> masterElemProvider = std::make_shared<stk::search::MasterElementProviderIntrepid2>();
  stk::search::SearchTopology hex8Topo(stk::topology::HEX_8);
 
  unsigned numRecvCopies = masterElemProvider->num_integration_points(hex8Topo);
  EXPECT_EQ(8u, numRecvCopies);
  auto recvBulk = stk::unit_test_util::create_mesh_with_field(comm, stk::topology::ELEM_RANK, numRecvCopies);

  stk::transfer_util::TransferMainSettings settings;
  settings.set_recv_type("ELEMENT_GAUSS_POINT");
  settings.set_transfer_send_parts({"block_1"});
  settings.set_transfer_recv_parts({"block_1"});

  stk::transfer_util::TransferMainBroker broker(comm, sendBulk, recvBulk, settings);

  broker.check_parts();
  broker.check_and_create_fields();

  stk::unit_test_util::LinearFieldEvaluator eval(3u);
  auto sendField = stk::mesh::get_field_by_name("field_1", sendBulk->mesh_meta_data());
  stk::unit_test_util::set_node_field(*sendBulk, *sendField, eval);

  broker.transfer();

  stk::unit_test_util::test_gauss_point_transfer(eval, recvBulk, stk::topology::ELEM_RANK, "field_1", masterElemProvider);
}

TEST(TransferMainBroker, elementSendMesh_masterElementInterpolationToGaussPointFieldList)
{
  stk::ParallelMachine comm = stk::parallel_machine_world();
  if(stk::parallel_machine_size(comm) != 1) { GTEST_SKIP(); }

  unsigned numSendCopies = 1;
  auto sendBulk = stk::unit_test_util::create_mesh_with_field(comm, stk::topology::NODE_RANK, numSendCopies);
  std::shared_ptr<stk::search::MasterElementProviderInterface> masterElemProvider = std::make_shared<stk::search::MasterElementProviderIntrepid2>();
  stk::search::SearchTopology hex8Topo(stk::topology::HEX_8);
 
  unsigned numRecvCopies = masterElemProvider->num_integration_points(hex8Topo);
  EXPECT_EQ(8u, numRecvCopies);
  auto recvBulk = stk::unit_test_util::create_mesh_with_field(comm, stk::topology::ELEM_RANK, numRecvCopies);

  stk::transfer_util::TransferMainSettings settings;
  std::vector<std::pair<std::string, std::string>> fieldList = {{"field_1", "recv_field_1"}};
  settings.set_transfer_field(fieldList[0]);
  settings.set_recv_type("ELEMENT_GAUSS_POINT");
  settings.set_transfer_send_parts({"block_1"});
  settings.set_transfer_recv_parts({"block_1"});

  stk::transfer_util::TransferMainBroker broker(comm, sendBulk, recvBulk, settings);

  broker.check_parts();
  broker.check_and_create_fields();

  stk::unit_test_util::LinearFieldEvaluator eval(3u);
  auto sendField = stk::mesh::get_field_by_name("field_1", sendBulk->mesh_meta_data());
  stk::unit_test_util::set_node_field(*sendBulk, *sendField, eval);

  broker.transfer();

  stk::unit_test_util::test_gauss_point_transfer(eval, recvBulk, stk::topology::ELEM_RANK, "recv_field_1", masterElemProvider);
}

TEST(TransferMainBroker, twoDimSendMeshQuad4_masterElementInterpolationToCentroid)
{
  stk::ParallelMachine comm = stk::parallel_machine_world();
  if(stk::parallel_machine_size(comm) != 1) { GTEST_SKIP(); }

  std::string mesh = "textmesh:"
                      "0,1,QUAD_4_2D,1,2,5,4,block_1\n"
                      "0,2,QUAD_4_2D,2,3,6,5,block_1\n"
                      "0,3,QUAD_4_2D,4,5,8,7,block_1\n"
                      "0,4,QUAD_4_2D,5,6,9,8,block_1\n"
                      "|coordinates:"
                      "0,0, 1,0, 2,0, 0,1, 1,1, 2,1, 0,2, 1,2, 2,2|dimension:2";

  auto sendBulk = stk::unit_test_util::create_mesh_with_field_2d(comm, stk::topology::NODE_RANK, mesh);
  auto recvBulk = stk::unit_test_util::create_mesh_with_field_2d(comm, stk::topology::ELEM_RANK, mesh);

  EXPECT_EQ(sendBulk->mesh_meta_data().spatial_dimension(), 2u);
  EXPECT_EQ(recvBulk->mesh_meta_data().spatial_dimension(), 2u);

  stk::transfer_util::TransferMainSettings settings;
  settings.set_recv_type("ELEMENT_CENTROID");
  settings.set_transfer_send_parts({"block_1"});
  settings.set_transfer_recv_parts({"block_1"});

  stk::transfer_util::TransferMainBroker broker(comm, sendBulk, recvBulk, settings);

  broker.check_parts();
  broker.check_and_create_fields();

  stk::unit_test_util::LinearFieldEvaluator eval(2u);
  auto sendField = stk::mesh::get_field_by_name("field_1", sendBulk->mesh_meta_data());
  stk::unit_test_util::set_node_field(*sendBulk, *sendField, eval);

  broker.transfer();

  stk::unit_test_util::test_centroid_transfer(eval, recvBulk, stk::topology::ELEM_RANK);
}

TEST(TransferMainBroker, twoDimSendMeshTri3_masterElementInterpolationToCentroid)
{
  stk::ParallelMachine comm = stk::parallel_machine_world();
  if(stk::parallel_machine_size(comm) != 1) { GTEST_SKIP(); }

  std::string mesh = "textmesh:" 
                     "0,1,TRI_3_2D,1,2,6, block_1\n"
                     "0,2,TRI_3_2D,2,7,6, block_1\n"
                     "0,3,TRI_3_2D,2,3,7, block_1\n"
                     "0,4,TRI_3_2D,4,9,8, block_1\n"
                     "0,5,TRI_3_2D,4,5,9, block_1\n"
                     "0,6,TRI_3_2D,5,10,9, block_1\n"
                     "|dimension:2"
                     "|coordinates:"
                     "0,0, 1,0, 2,0, 3,0, 4,0, 0,1, 1,1, 2,1, 3,1, 4,1";

  auto sendBulk = stk::unit_test_util::create_mesh_with_field_2d(comm, stk::topology::NODE_RANK, mesh);
  auto recvBulk = stk::unit_test_util::create_mesh_with_field_2d(comm, stk::topology::ELEM_RANK, mesh);

  EXPECT_EQ(sendBulk->mesh_meta_data().spatial_dimension(), 2u);
  EXPECT_EQ(recvBulk->mesh_meta_data().spatial_dimension(), 2u);

  stk::transfer_util::TransferMainSettings settings;
  settings.set_recv_type("ELEMENT_CENTROID");
  settings.set_transfer_send_parts({"block_1"});
  settings.set_transfer_recv_parts({"block_1"});

  stk::transfer_util::TransferMainBroker broker(comm, sendBulk, recvBulk, settings);

  broker.check_parts();
  broker.check_and_create_fields();

  stk::unit_test_util::LinearFieldEvaluator eval(2u);
  auto sendField = stk::mesh::get_field_by_name("field_1", sendBulk->mesh_meta_data());
  stk::unit_test_util::set_node_field(*sendBulk, *sendField, eval);

  broker.transfer();

  stk::unit_test_util::test_centroid_transfer(eval, recvBulk, stk::topology::ELEM_RANK);
}

}
