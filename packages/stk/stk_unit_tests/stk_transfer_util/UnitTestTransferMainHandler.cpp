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
#include "stk_transfer_util/TransferMainHandler.hpp"
#include "stk_mesh/base/MetaData.hpp"
#include "stk_mesh/base/BulkData.hpp"
#include "stk_mesh/base/MeshBuilder.hpp"
#include "stk_unit_test_utils/ioUtils.hpp"
#include "stk_unit_test_utils/CommandLineArgs.hpp"
#include "stk_unit_test_utils/FieldEvaluator.hpp"
#include "stk_unit_test_utils/TransferTesters.hpp"
#include "stk_io/FillMesh.hpp"
#include "stk_io/IossBridge.hpp"
#include "stk_io/StkMeshIoBroker.hpp"
#include "stk_mesh/base/GetEntities.hpp"
#include "stk_search_util/ObjectCoordinates.hpp"
#include <stk_util/parallel/OutputStreams.hpp>
#include <stk_util/stk_config.h>
#include <fstream>

namespace {

void build_serial_mesh(const std::string fileName, const std::string fieldName = "field_1") {
  stk::unit_test_util::generated_mesh_with_transient_data_to_file_in_serial("4x1x1", fileName,
                                                          fieldName, stk::topology::NODE_RANK,
                                                          "global_var", {1},
                                                          stk::unit_test_util::IdAndTimeFieldValueSetter());
}

TEST(TransferMainHandler, basic)
{
  stk::ParallelMachine comm = stk::parallel_machine_world();
  if(stk::parallel_machine_size(comm) != 1) { GTEST_SKIP(); }

  std::string sendMesh = "source.exo";
  std::string recvMesh = "target.exo";

  build_serial_mesh(sendMesh);
  build_serial_mesh(recvMesh);
  stk::unit_test_util::Args args({"stk_transfer", "--send-mesh", sendMesh, "--recv-mesh", recvMesh, "--use-master-elements", "intrepid2"});

  stk::set_outputP0(&stk::outputNull(), comm);
  stk::transfer_util::TransferMainHandler handlerBasic(comm, args.argc(), args.argv());
  handlerBasic.run();

  EXPECT_EQ(handlerBasic.exit_code(), stk::transfer_util::TransferMainStatus::SUCCESS);

  unlink(sendMesh.c_str());
  unlink(recvMesh.c_str());
}

TEST(TransferMainHandler, outputMeshHasPreExistingRecvMeshFields)
{
  stk::ParallelMachine comm = stk::parallel_machine_world();
  if(stk::parallel_machine_size(comm) != 1) { GTEST_SKIP(); }

  std::string sendMesh = "source.exo";
  std::string recvMesh = "target.exo";
  std::string outputMesh = "output.exo";

  build_serial_mesh(sendMesh);
  std::string field2("field_2");
  build_serial_mesh(recvMesh, field2);
  stk::unit_test_util::Args args({"stk_transfer", "--send-mesh", sendMesh, "--recv-mesh", recvMesh, "--use-master-elements", "intrepid2", "--output-mesh", outputMesh});

  stk::set_outputP0(&stk::outputNull(), comm);
  stk::transfer_util::TransferMainHandler handlerBasic(comm, args.argc(), args.argv());
  handlerBasic.run();

  EXPECT_EQ(handlerBasic.exit_code(), stk::transfer_util::TransferMainStatus::SUCCESS);

  {
    auto output = stk::mesh::MeshBuilder(comm).create();
    stk::io::fill_mesh_with_fields(outputMesh, *output);
    const stk::mesh::MetaData& meta = output->mesh_meta_data();

    const auto* preExistingScalarFieldOnRecvMesh = meta.get_field(stk::topology::NODE_RANK, field2+"_scalar");
    EXPECT_TRUE(preExistingScalarFieldOnRecvMesh != nullptr);

    const auto* preExistingVectorFieldOnRecvMesh = meta.get_field(stk::topology::NODE_RANK, field2+"_vector");
    EXPECT_TRUE(preExistingVectorFieldOnRecvMesh != nullptr);
  }

  unlink(sendMesh.c_str());
  unlink(recvMesh.c_str());
  unlink(outputMesh.c_str());
}

TEST(TransferMainHandler, all)
{
  stk::ParallelMachine comm = stk::parallel_machine_world();
  if(stk::parallel_machine_size(comm) != 1) { GTEST_SKIP(); }

  std::string sendMesh = "source.exo";
  std::string recvMesh = "target.exo";

  build_serial_mesh(sendMesh);
  build_serial_mesh(recvMesh);
  stk::unit_test_util::Args args({"stk_transfer", "--send-mesh", sendMesh, "--recv-mesh", recvMesh,
                                 "--field-list", "field_1_scalar", "--part-list", "block_1:block_1", 
                                 "--extrapolate-option", "PROJECT", "--xfer-type", "INTERP",
                                 "--use-master-elements", "intrepid2" });
  
  stk::set_outputP0(&stk::outputNull(), comm);
  stk::transfer_util::TransferMainHandler handlerAll(comm, args.argc(), args.argv());
  handlerAll.run();

  EXPECT_EQ(handlerAll.exit_code(), stk::transfer_util::TransferMainStatus::SUCCESS);

  unlink(sendMesh.c_str());
  unlink(recvMesh.c_str());
}

#ifdef STK_HAS_YAML_CPP

void write_mesh_with_node_and_elem_fields(stk::ParallelMachine comm, const std::string& meshSpec,
                                          const std::string& fileName, const stk::unit_test_util::FieldEvaluator& eval)
{
  std::shared_ptr<stk::mesh::BulkData> bulk = stk::mesh::MeshBuilder(comm).set_spatial_dimension(3).create();
  stk::mesh::MetaData& meta = bulk->mesh_meta_data();

  stk::mesh::Field<double>& nodeField = meta.declare_field<double>(stk::topology::NODE_RANK, "nodeField", 1);
  stk::mesh::put_field_on_mesh(nodeField, meta.universal_part(), nullptr);
  stk::io::set_field_role(nodeField, Ioss::Field::TRANSIENT);

  stk::mesh::Field<double>& elemField = meta.declare_field<double>(stk::topology::ELEM_RANK, "elemField", 1);
  stk::mesh::put_field_on_mesh(elemField, meta.universal_part(), nullptr);
  stk::io::set_field_role(elemField, Ioss::Field::TRANSIENT);

  stk::io::fill_mesh("generated:" + meshSpec, *bulk);

  stk::unit_test_util::set_node_field(*bulk, nodeField, eval);
  stk::unit_test_util::set_entity_field(*bulk, elemField, eval);

  stk::io::StkMeshIoBroker ioBroker(comm);
  ioBroker.set_bulk_data(bulk);
  size_t outputIndex = ioBroker.create_output_mesh(fileName, stk::io::WRITE_RESULTS);
  ioBroker.add_field(outputIndex, nodeField);
  ioBroker.add_field(outputIndex, elemField);
  ioBroker.process_output_request(outputIndex, 1.0);
}

void expect_field_matches_eval(const stk::mesh::BulkData& bulk, const std::string& fieldName,
                               stk::mesh::EntityRank rank, const stk::unit_test_util::FieldEvaluator& eval)
{
  const stk::mesh::MetaData& meta = bulk.mesh_meta_data();
  const stk::mesh::FieldBase* field = meta.get_field(rank, fieldName);
  ASSERT_TRUE(field != nullptr) << fieldName;

  stk::mesh::EntityVector entities;
  stk::mesh::get_selected_entities(stk::mesh::Selector(*field), bulk.buckets(rank), entities);
  ASSERT_FALSE(entities.empty());

  std::vector<double> location;
  for (stk::mesh::Entity entity : entities) {
    stk::search::determine_centroid(meta.spatial_dimension(), entity, *meta.coordinate_field(), location);
    stk::unit_test_util::test_expected_value_at_location(bulk, field, entity, eval, location.data());
  }
}

TEST(TransferMainHandler, inputFileFieldGroupsWithDifferentTransferTypes)
{
  stk::ParallelMachine comm = stk::parallel_machine_world();
  if(stk::parallel_machine_size(comm) != 1) { GTEST_SKIP(); }

  const std::string sendMesh = "fieldGroupsSource.exo";
  const std::string recvMesh = "fieldGroupsTarget.exo";
  const std::string outputMesh = "fieldGroupsOutput.exo";
  const std::string inputFile = "fieldGroups.yaml";

  stk::unit_test_util::LinearFieldEvaluator eval(3u);
  write_mesh_with_node_and_elem_fields(comm, "2x2x2", sendMesh, eval);
  stk::unit_test_util::generated_mesh_to_file_in_serial("3x3x3|bbox:0,0,0,2,2,2", recvMesh);

  {
    std::ofstream yaml(inputFile);
    yaml << "send-mesh: " << sendMesh << "\n"
         << "recv-mesh: " << recvMesh << "\n"
         << "output-mesh: " << outputMesh << "\n"
         << "use-master-elements: intrepid2\n"
         << "field-list:\n"
         << "  - fields: nodeField\n"
         << "    xfer-type: INTERP\n"
         << "    recv-type: NODE\n"
         << "  - fields: elemField\n"
         << "    xfer-type: PATCH\n"
         << "    recv-type: ELEMENT_CENTROID\n";
  }

  stk::unit_test_util::Args args({"stk_transfer", "--input-file", inputFile});
  stk::set_outputP0(&stk::outputNull(), comm);
  stk::transfer_util::TransferMainHandler handler(comm, args.argc(), args.argv());
  ASSERT_EQ(2u, handler.get_transfer_settings().get_field_groups().size());
  handler.run();
  stk::reset_default_output_streams(comm);

  EXPECT_EQ(handler.exit_code(), stk::transfer_util::TransferMainStatus::SUCCESS);

  {
    auto output = stk::mesh::MeshBuilder(comm).create();
    stk::io::fill_mesh_with_fields(outputMesh, *output);
    expect_field_matches_eval(*output, "nodeField", stk::topology::NODE_RANK, eval);
    expect_field_matches_eval(*output, "elemField", stk::topology::ELEM_RANK, eval);
  }

  unlink(sendMesh.c_str());
  unlink(recvMesh.c_str());
  unlink(outputMesh.c_str());
  unlink(inputFile.c_str());
}

TEST(TransferMainHandler, inputFileFieldGroupsWithSameRecvFieldIsExecutionError)
{
  stk::ParallelMachine comm = stk::parallel_machine_world();
  if(stk::parallel_machine_size(comm) != 1) { GTEST_SKIP(); }

  const std::string sendMesh = "fieldGroupsDupSource.exo";
  const std::string recvMesh = "fieldGroupsDupTarget.exo";
  const std::string inputFile = "fieldGroupsDup.yaml";

  stk::unit_test_util::LinearFieldEvaluator eval(3u);
  write_mesh_with_node_and_elem_fields(comm, "2x2x2", sendMesh, eval);
  stk::unit_test_util::generated_mesh_to_file_in_serial("2x2x2", recvMesh);

  {
    std::ofstream yaml(inputFile);
    yaml << "send-mesh: " << sendMesh << "\n"
         << "recv-mesh: " << recvMesh << "\n"
         << "use-master-elements: intrepid2\n"
         << "field-list:\n"
         << "  - {fields: nodeField, xfer-type: INTERP}\n"
         << "  - {fields: nodeField, xfer-type: COPY}\n";
  }

  stk::unit_test_util::Args args({"stk_transfer", "--input-file", inputFile});
  stk::set_outputP0(&stk::outputNull(), comm);
  stk::transfer_util::TransferMainHandler handler(comm, args.argc(), args.argv());
  testing::internal::CaptureStderr();
  handler.run();
  const std::string errorOutput = testing::internal::GetCapturedStderr();
  stk::reset_default_output_streams(comm);

  EXPECT_EQ(handler.exit_code(), stk::transfer_util::TransferMainStatus::EXECUTION_ERROR);
  EXPECT_NE(errorOutput.find("more than one field group"), std::string::npos) << errorOutput;

  unlink(sendMesh.c_str());
  unlink(recvMesh.c_str());
  unlink(inputFile.c_str());
}

#endif // STK_HAS_YAML_CPP

TEST(TransferMainHandler, noArgs)
{
  stk::ParallelMachine comm = stk::parallel_machine_world();
  if(stk::parallel_machine_size(comm) != 1) { GTEST_SKIP(); }

  stk::unit_test_util::Args args({"stk_transfer"});

  testing::internal::CaptureStderr();
  stk::set_outputP0(&stk::outputNull(), comm);
  stk::transfer_util::TransferMainHandler handlerNoArgs(comm, args.argc(), args.argv());
  handlerNoArgs.run();
  testing::internal::GetCapturedStderr();

  EXPECT_EQ(handlerNoArgs.exit_code(), stk::transfer_util::TransferMainStatus::PARSE_ERROR);
}

TEST(TransferMainHandler, missingRequiredFromMesh)
{
  stk::ParallelMachine comm = stk::parallel_machine_world();
  if(stk::parallel_machine_size(comm) != 1) { GTEST_SKIP(); }

  std::string recvMesh = "target.exo";
  build_serial_mesh(recvMesh);
  stk::unit_test_util::Args args({"stk_transfer", "--recv-mesh", recvMesh});

  testing::internal::CaptureStderr();
  stk::set_outputP0(&stk::outputNull(), comm);
  stk::transfer_util::TransferMainHandler handlerNoFrom(comm, args.argc(), args.argv());
  handlerNoFrom.run();
  testing::internal::GetCapturedStderr();

  EXPECT_EQ(handlerNoFrom.exit_code(), stk::transfer_util::TransferMainStatus::PARSE_ERROR);

  unlink(recvMesh.c_str());
}

TEST(TransferMainHandler, missingRequiredToMesh)
{
  stk::ParallelMachine comm = stk::parallel_machine_world();
  if(stk::parallel_machine_size(comm) != 1) { GTEST_SKIP(); }

  std::string sendMesh = "source.exo";
  build_serial_mesh(sendMesh);
  stk::unit_test_util::Args args({"stk_transfer", "--send-mesh", sendMesh});

  testing::internal::CaptureStderr();
  stk::set_outputP0(&stk::outputNull(), comm);
  stk::transfer_util::TransferMainHandler handlerNoTo(comm, args.argc(), args.argv());
  handlerNoTo.run();
  testing::internal::GetCapturedStderr();

  EXPECT_EQ(handlerNoTo.exit_code(), stk::transfer_util::TransferMainStatus::PARSE_ERROR);

  unlink(sendMesh.c_str());
}

TEST(TransferMainHandler, invalidExtrapolateOption)
{
  stk::ParallelMachine comm = stk::parallel_machine_world();
  if(stk::parallel_machine_size(comm) != 1) { GTEST_SKIP(); }

  std::string sendMesh = "source.exo";
  std::string recvMesh = "target.exo";

  build_serial_mesh(sendMesh);
  build_serial_mesh(recvMesh);
  stk::unit_test_util::Args args({"stk_transfer", "--send-mesh", sendMesh, "--recv-mesh", recvMesh,
                                 "--field-list", "field_1", "--extrapolate-option", "FAKE" });

  testing::internal::CaptureStderr();
  stk::set_outputP0(&stk::outputNull(), comm);
  stk::transfer_util::TransferMainHandler handlerAll(comm, args.argc(), args.argv());
  EXPECT_EQ(stk::transfer_util::TransferMainStatus::PARSE_ERROR, handlerAll.exit_code());
  testing::internal::GetCapturedStderr();

  unlink(sendMesh.c_str());
  unlink(recvMesh.c_str());
}

TEST(TransferMainHandler, invalidTransferType)
{
  stk::ParallelMachine comm = stk::parallel_machine_world();
  if(stk::parallel_machine_size(comm) != 1) { GTEST_SKIP(); }

  std::string sendMesh = "source.exo";
  std::string recvMesh = "target.exo";

  build_serial_mesh(sendMesh);
  build_serial_mesh(recvMesh);
  stk::unit_test_util::Args args({"stk_transfer", "--send-mesh", sendMesh, "--recv-mesh", recvMesh,
                                 "--field-list", "field_1", "--xfer-type", "FAKE" });

  testing::internal::CaptureStderr();
  stk::set_outputP0(&stk::outputNull(), comm);
  stk::transfer_util::TransferMainHandler handlerAll(comm, args.argc(), args.argv());
  EXPECT_EQ(stk::transfer_util::TransferMainStatus::PARSE_ERROR, handlerAll.exit_code());
  testing::internal::GetCapturedStderr();

  unlink(sendMesh.c_str());
  unlink(recvMesh.c_str());
}

TEST(TransferMainHandler, invalidRecvType)
{
  stk::ParallelMachine comm = stk::parallel_machine_world();
  if(stk::parallel_machine_size(comm) != 1) { GTEST_SKIP(); }

  std::string sendMesh = "source.exo";
  std::string recvMesh = "target.exo";

  build_serial_mesh(sendMesh);
  build_serial_mesh(recvMesh);
  stk::unit_test_util::Args args({"stk_transfer", "--send-mesh", sendMesh, "--recv-mesh", recvMesh,
                                 "--field-list", "field_1", "--recv-type", "FAKE" });

  testing::internal::CaptureStderr();
  stk::set_outputP0(&stk::outputNull(), comm);
  stk::transfer_util::TransferMainHandler handlerAll(comm, args.argc(), args.argv());
  EXPECT_EQ(stk::transfer_util::TransferMainStatus::PARSE_ERROR, handlerAll.exit_code());
  testing::internal::GetCapturedStderr();

  unlink(sendMesh.c_str());
  unlink(recvMesh.c_str());
}

TEST(TransferMainHandler, sameFromAndToMesh)
{
  stk::ParallelMachine comm = stk::parallel_machine_world();
  if(stk::parallel_machine_size(comm) != 1) { GTEST_SKIP(); }

  std::string sendMesh = "source.exo";
  build_serial_mesh(sendMesh);
  stk::unit_test_util::Args args({"stk_transfer", "--send-mesh", sendMesh, "--recv-mesh", sendMesh});

  stk::set_outputP0(&stk::outputNull(), comm);
  stk::transfer_util::TransferMainHandler handlerSameFromTo(comm, args.argc(), args.argv());
  handlerSameFromTo.run();

  EXPECT_EQ(handlerSameFromTo.exit_code(), stk::transfer_util::TransferMainStatus::SUCCESS);

  unlink(sendMesh.c_str());
}

TEST(TransferMainHandler, help)
{
  stk::ParallelMachine comm = stk::parallel_machine_world();
  if(stk::parallel_machine_size(comm) != 1) { GTEST_SKIP(); }

  stk::unit_test_util::Args args({"stk_transfer", "--help"});

  stk::set_outputP0(&stk::outputNull(), comm);
  stk::transfer_util::TransferMainHandler handlerHelp(comm, args.argc(), args.argv());
  handlerHelp.run();

  EXPECT_EQ(handlerHelp.exit_code(), stk::transfer_util::TransferMainStatus::PARSE_ONLY);
}

TEST(TransferMainHandler, version)
{
  stk::ParallelMachine comm = stk::parallel_machine_world();
  if(stk::parallel_machine_size(comm) != 1) { GTEST_SKIP(); }

  stk::unit_test_util::Args args({"stk_transfer", "--version"});

  stk::set_outputP0(&stk::outputNull(), comm);
  stk::transfer_util::TransferMainHandler handlerVersion(comm, args.argc(), args.argv());
  handlerVersion.run();

  EXPECT_EQ(handlerVersion.exit_code(), stk::transfer_util::TransferMainStatus::PARSE_ONLY);
}

}
