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
#include "stk_unit_test_utils/getOption.h"            // for get_command_lin...
#include "stk_util/parallel/Parallel.hpp"             // for parallel_machin...
#include "stk_util/util/ReportHandler.hpp"            // for ThrowRequireMsg
#include "stk_util/command_line/CommandLineParserUtils.hpp"
#include "stk_unit_test_utils/CommandLineArgs.hpp"
#include "stk_transfer_util/TransferMainParser.hpp"
#include "stk_transfer_util/TransferMainSettings.hpp"
#include "stk_util/parallel/OutputStreams.hpp"

namespace {

class SimpleTransferMainParser : public stk::transfer_util::TransferMainParser {
  public:

  SimpleTransferMainParser(MPI_Comm comm)
    : TransferMainParser(comm)
  {}

  unsigned get_num_input_processors() const { return 1; }

  unsigned get_num_output_processors() const { return 1; }

  std::string get_sendMesh_filename() const
  {
    return m_cmdLineParser.is_option_parsed(m_optionNames.sendMesh) ?
      m_cmdLineParser.get_option_value<std::string>(m_optionNames.sendMesh) : "not parsed";
  }

  std::string get_recvMesh_filename() const
  {
    return m_cmdLineParser.is_option_parsed(m_optionNames.recvMesh) ?
      m_cmdLineParser.get_option_value<std::string>(m_optionNames.recvMesh) : "not parsed";
  }

  std::string get_outputMesh_filename() const
  {
    return m_cmdLineParser.is_option_parsed(m_optionNames.outputMesh) ?
      m_cmdLineParser.get_option_value<std::string>(m_optionNames.outputMesh) : "not parsed";
  }

  std::string get_transfer_fields() const
  {
    return m_cmdLineParser.is_option_parsed(m_optionNames.fieldList) ?
      m_cmdLineParser.get_option_value<std::string>(m_optionNames.fieldList) : "defaulting to all fields";
  }

  std::string get_transfer_parts() const
  {
    return m_cmdLineParser.is_option_parsed(m_optionNames.partList) ?
      m_cmdLineParser.get_option_value<std::string>(m_optionNames.partList) : "defaulting to all parts";
  }

  std::string get_extrapolate_option() const
  {
    return m_cmdLineParser.is_option_parsed(m_optionNames.ExtrapolateOption) ?
      m_cmdLineParser.get_option_value<std::string>(m_optionNames.ExtrapolateOption) : "defaulting to IGNORE";
  }

  std::string get_master_elements_name() const
  {
    return m_cmdLineParser.is_option_parsed(m_optionNames.UseMasterElements) ?
      m_cmdLineParser.get_option_value<std::string>(m_optionNames.UseMasterElements) : "defaulting to intrepid2";
  }

  std::string get_recv_type() const
  {
    return m_cmdLineParser.is_option_parsed(m_optionNames.recvType) ?
      m_cmdLineParser.get_option_value<std::string>(m_optionNames.recvType) : "defaulting to NODE";
  }

  std::string get_time_steps() const
  {
    return m_cmdLineParser.is_option_parsed(m_optionNames.timeSteps) ?
      m_cmdLineParser.get_option_value<std::string>(m_optionNames.timeSteps) : "defaulting to ALL";
  }

  std::string get_transfer_type() const
  {
    return m_cmdLineParser.is_option_parsed(m_optionNames.transferType) ?
      m_cmdLineParser.get_option_value<std::string>(m_optionNames.transferType) : "defaulting to INTERP";
  }

  std::string get_2d_to_3d_mapping_type() const
  {
    return m_cmdLineParser.is_option_parsed(m_optionNames.mappingType2dTo3d) ?
      m_cmdLineParser.get_option_value<std::string>(m_optionNames.mappingType2dTo3d) : "defaulting to ZPLANE";
  }

  std::string get_coord_transf_z_expr() const
  {
    return m_cmdLineParser.is_option_parsed(m_optionNames.coordTransformZExpr) ?
    m_cmdLineParser.get_option_value<std::string>(m_optionNames.coordTransformZExpr) : "defaulting to z=0";
  }

  std::string get_quickExample() const { return m_quickExample; }

  std::string get_longExamples() const { return m_longExamples; }

};

TEST(SimpleTransferMainParser, empty)
{
  stk::ParallelMachine comm = stk::parallel_machine_world();
  if(stk::parallel_machine_size(comm) != 1) { GTEST_SKIP(); }

  stk::unit_test_util::Args args({"stk_transfer"});

  stk::set_outputP0(&stk::outputNull(), comm);
  SimpleTransferMainParser emptyParser(comm);

  EXPECT_EQ("not parsed", emptyParser.get_sendMesh_filename());
  EXPECT_EQ("not parsed", emptyParser.get_recvMesh_filename());
  EXPECT_EQ("defaulting to all fields", emptyParser.get_transfer_fields());
  EXPECT_EQ("defaulting to all parts", emptyParser.get_transfer_parts());
  EXPECT_EQ("defaulting to IGNORE", emptyParser.get_extrapolate_option());
  EXPECT_EQ("defaulting to intrepid2", emptyParser.get_master_elements_name());
  EXPECT_EQ("defaulting to NODE", emptyParser.get_recv_type());
  EXPECT_EQ("defaulting to ALL", emptyParser.get_time_steps());
  EXPECT_EQ("defaulting to INTERP", emptyParser.get_transfer_type());
  EXPECT_EQ("defaulting to ZPLANE", emptyParser.get_2d_to_3d_mapping_type());
  EXPECT_EQ("defaulting to z=0", emptyParser.get_coord_transf_z_expr());

  testing::internal::CaptureStderr();
  emptyParser.parse_command_line_options(args.argc(), args.argv());
  testing::internal::GetCapturedStderr();

  EXPECT_EQ(emptyParser.get_parser_status(), stk::transfer_util::TransferParserStatus::PARSE_ERROR);
}

TEST(SimpleTransferMainParser, basic)
{
  stk::ParallelMachine comm = stk::parallel_machine_world();
  if(stk::parallel_machine_size(comm) != 1) { GTEST_SKIP(); }

  stk::unit_test_util::Args args({"stk_transfer", "--send-mesh", "source.exo", "--recv-mesh", "target.exo"});

  stk::set_outputP0(&stk::outputNull(), comm);
  SimpleTransferMainParser basicParser(comm);

  EXPECT_EQ("not parsed", basicParser.get_sendMesh_filename());
  EXPECT_EQ("not parsed", basicParser.get_recvMesh_filename());
  EXPECT_EQ("defaulting to all fields", basicParser.get_transfer_fields());
  EXPECT_EQ("defaulting to all parts", basicParser.get_transfer_parts());
  EXPECT_EQ("defaulting to IGNORE", basicParser.get_extrapolate_option());
  EXPECT_EQ("defaulting to intrepid2", basicParser.get_master_elements_name());
  EXPECT_EQ("defaulting to NODE", basicParser.get_recv_type());
  EXPECT_EQ("defaulting to ALL", basicParser.get_time_steps());
  EXPECT_EQ("defaulting to INTERP", basicParser.get_transfer_type());
  EXPECT_EQ("defaulting to ZPLANE", basicParser.get_2d_to_3d_mapping_type());
  EXPECT_EQ("defaulting to z=0", basicParser.get_coord_transf_z_expr());

  basicParser.parse_command_line_options(args.argc(), args.argv());

  EXPECT_EQ(basicParser.get_parser_status(), stk::transfer_util::TransferParserStatus::SUCCESS);
  EXPECT_EQ("source.exo", basicParser.get_sendMesh_filename());
  EXPECT_EQ("target.exo", basicParser.get_recvMesh_filename());
  EXPECT_EQ("defaulting to all fields", basicParser.get_transfer_fields());
  EXPECT_EQ("defaulting to all parts", basicParser.get_transfer_parts());
  EXPECT_EQ("defaulting to IGNORE", basicParser.get_extrapolate_option());
  EXPECT_EQ("defaulting to intrepid2", basicParser.get_master_elements_name());
  EXPECT_EQ("defaulting to ALL", basicParser.get_time_steps());
  EXPECT_EQ("defaulting to INTERP", basicParser.get_transfer_type());
  EXPECT_EQ("defaulting to ZPLANE", basicParser.get_2d_to_3d_mapping_type());
  EXPECT_EQ("defaulting to z=0", basicParser.get_coord_transf_z_expr());
  EXPECT_EQ(1u, basicParser.get_num_input_processors());
  EXPECT_EQ(1u, basicParser.get_num_output_processors());

}

TEST(SimpleTransferMainParser, basicAll)
{
  stk::ParallelMachine comm = stk::parallel_machine_world();
  if(stk::parallel_machine_size(comm) != 1) { GTEST_SKIP(); }

  stk::unit_test_util::Args args({"stk_transfer", "--send-mesh", "source.exo", "--recv-mesh", "target.exo",
                                  "--field-list", "field_1", "--part-list", "send_part1:recv_part1",
                                  "--extrapolate-option", "EXTRAPOLATE", "--use-master-elements", "some_name",
                                  "--recv-type", "RECV_TYPE", "--time-steps", "LAST", "--xfer-type", "PATCH",
                                  "--2d-to-3d-mapping-type", "EXTRUDE", "--2d-to-3d-z-coord", "z=x+y"});

  stk::set_outputP0(&stk::outputNull(), comm);
  SimpleTransferMainParser basicAllParser(comm);

  EXPECT_EQ("not parsed", basicAllParser.get_sendMesh_filename());
  EXPECT_EQ("not parsed", basicAllParser.get_recvMesh_filename());
  EXPECT_EQ("not parsed", basicAllParser.get_outputMesh_filename());
  EXPECT_EQ("defaulting to all fields", basicAllParser.get_transfer_fields());
  EXPECT_EQ("defaulting to all parts", basicAllParser.get_transfer_parts());
  EXPECT_EQ("defaulting to IGNORE", basicAllParser.get_extrapolate_option());
  EXPECT_EQ("defaulting to intrepid2", basicAllParser.get_master_elements_name());
  EXPECT_EQ("defaulting to NODE", basicAllParser.get_recv_type());
  EXPECT_EQ("defaulting to ALL", basicAllParser.get_time_steps());
  EXPECT_EQ("defaulting to INTERP", basicAllParser.get_transfer_type());
  EXPECT_EQ("defaulting to ZPLANE", basicAllParser.get_2d_to_3d_mapping_type());
  EXPECT_EQ("defaulting to z=0", basicAllParser.get_coord_transf_z_expr());

  basicAllParser.parse_command_line_options(args.argc(), args.argv());

  EXPECT_EQ(basicAllParser.get_parser_status(), stk::transfer_util::TransferParserStatus::SUCCESS);
  EXPECT_EQ("source.exo", basicAllParser.get_sendMesh_filename());
  EXPECT_EQ("target.exo", basicAllParser.get_recvMesh_filename());
  EXPECT_EQ("not parsed", basicAllParser.get_outputMesh_filename());
  EXPECT_EQ("field_1", basicAllParser.get_transfer_fields());
  EXPECT_EQ("send_part1:recv_part1", basicAllParser.get_transfer_parts());
  EXPECT_EQ("EXTRAPOLATE", basicAllParser.get_extrapolate_option());
  EXPECT_EQ("some_name", basicAllParser.get_master_elements_name());
  EXPECT_EQ("RECV_TYPE", basicAllParser.get_recv_type());
  EXPECT_EQ("LAST", basicAllParser.get_time_steps());
  EXPECT_EQ("PATCH", basicAllParser.get_transfer_type());
  EXPECT_EQ(1u, basicAllParser.get_num_input_processors());
  EXPECT_EQ(1u, basicAllParser.get_num_output_processors());
  EXPECT_EQ("EXTRUDE", basicAllParser.get_2d_to_3d_mapping_type());
  EXPECT_EQ("z=x+y", basicAllParser.get_coord_transf_z_expr());
}

TEST(SimpleTransferMainParser, basicAllShort)
{
  stk::ParallelMachine comm = stk::parallel_machine_world();
  if(stk::parallel_machine_size(comm) != 1) { GTEST_SKIP(); }

  stk::unit_test_util::Args args({"stk_transfer", "-s", "source.exo", "-r", "target.exo", "-o", "output.exo",
                                  "-f", "field_a,field_b", "-p", "part_a:part_b", "-e", "TRUNCATE", "-u", "some_name",
                                  "-t", "RECV_TYPE", "-x", "COPY"});

  stk::set_outputP0(&stk::outputNull(), comm);
  SimpleTransferMainParser basicAllShortParser(comm);

  EXPECT_EQ("not parsed", basicAllShortParser.get_sendMesh_filename());
  EXPECT_EQ("not parsed", basicAllShortParser.get_recvMesh_filename());
  EXPECT_EQ("not parsed", basicAllShortParser.get_outputMesh_filename());
  EXPECT_EQ("defaulting to all fields", basicAllShortParser.get_transfer_fields());
  EXPECT_EQ("defaulting to all parts", basicAllShortParser.get_transfer_parts());
  EXPECT_EQ("defaulting to IGNORE", basicAllShortParser.get_extrapolate_option());
  EXPECT_EQ("defaulting to intrepid2", basicAllShortParser.get_master_elements_name());
  EXPECT_EQ("defaulting to NODE", basicAllShortParser.get_recv_type());
  EXPECT_EQ("defaulting to INTERP", basicAllShortParser.get_transfer_type());

  basicAllShortParser.parse_command_line_options(args.argc(), args.argv());

  EXPECT_EQ(basicAllShortParser.get_parser_status(), stk::transfer_util::TransferParserStatus::SUCCESS);
  EXPECT_EQ("source.exo", basicAllShortParser.get_sendMesh_filename());
  EXPECT_EQ("target.exo", basicAllShortParser.get_recvMesh_filename());
  EXPECT_EQ("output.exo", basicAllShortParser.get_outputMesh_filename());
  EXPECT_EQ("field_a,field_b", basicAllShortParser.get_transfer_fields());
  EXPECT_EQ("part_a:part_b", basicAllShortParser.get_transfer_parts());
  EXPECT_EQ("TRUNCATE", basicAllShortParser.get_extrapolate_option());
  EXPECT_EQ("some_name", basicAllShortParser.get_master_elements_name());
  EXPECT_EQ("RECV_TYPE", basicAllShortParser.get_recv_type());
  EXPECT_EQ("COPY", basicAllShortParser.get_transfer_type());
  EXPECT_EQ(1u, basicAllShortParser.get_num_input_processors());
  EXPECT_EQ(1u, basicAllShortParser.get_num_output_processors());
}

TEST(SimpleTransferMainParser, invalidArg)
{
  stk::ParallelMachine comm = stk::parallel_machine_world();
  if(stk::parallel_machine_size(comm) != 1) { GTEST_SKIP(); }

  stk::unit_test_util::Args args({"stk_transfer", "--send-mesh", "source.exo", "--recv-mesh", "target.exo",
                                  "--field-vector", "field_1", "--extrapolate-option", "EXTRAPOLATE"});

  stk::set_outputP0(&stk::outputNull(), comm);
  SimpleTransferMainParser invalidArgParser(comm);

  EXPECT_EQ("not parsed", invalidArgParser.get_sendMesh_filename());
  EXPECT_EQ("not parsed", invalidArgParser.get_recvMesh_filename());
  EXPECT_EQ("defaulting to all fields", invalidArgParser.get_transfer_fields());
  EXPECT_EQ("defaulting to IGNORE", invalidArgParser.get_extrapolate_option());

  testing::internal::CaptureStderr();
  invalidArgParser.parse_command_line_options(args.argc(), args.argv());
  testing::internal::GetCapturedStderr();

  EXPECT_EQ(invalidArgParser.get_parser_status(), stk::transfer_util::TransferParserStatus::PARSE_ERROR);
}

TEST(SimpleTransferMainParser, invalidArgShort)
{
  stk::ParallelMachine comm = stk::parallel_machine_world();
  if(stk::parallel_machine_size(comm) != 1) { GTEST_SKIP(); }

  stk::unit_test_util::Args args({"stk_transfer", "-s", "source.exo", "-r", "target.exo",
                                  "-f", "field_a", "-z", "TRUNCATE"});

  stk::set_outputP0(&stk::outputNull(), comm);
  SimpleTransferMainParser invalidArgShortParser(comm);

  EXPECT_EQ("not parsed", invalidArgShortParser.get_sendMesh_filename());
  EXPECT_EQ("not parsed", invalidArgShortParser.get_recvMesh_filename());
  EXPECT_EQ("defaulting to all fields", invalidArgShortParser.get_transfer_fields());
  EXPECT_EQ("defaulting to IGNORE", invalidArgShortParser.get_extrapolate_option());

  testing::internal::CaptureStderr();
  invalidArgShortParser.parse_command_line_options(args.argc(), args.argv());
  testing::internal::GetCapturedStderr();

  EXPECT_EQ(invalidArgShortParser.get_parser_status(), stk::transfer_util::TransferParserStatus::PARSE_ERROR);
}

TEST(SimpleTransferMainParser, checkExamples)
{
  stk::ParallelMachine comm = stk::parallel_machine_world();
  if(stk::parallel_machine_size(comm) != 1) { GTEST_SKIP(); }

  stk::unit_test_util::Args args({"stk_transfer", "-f", "source.exo", "-t", "target.exo"});

  stk::set_outputP0(&stk::outputNull(), comm);
  SimpleTransferMainParser examplesParser(comm);

  EXPECT_TRUE(examplesParser.get_quickExample().starts_with("Basic example"));
  EXPECT_TRUE(examplesParser.get_longExamples().starts_with("Examples"));
}

TEST(SimpleTransferMainParser, help)
{
  stk::ParallelMachine comm = stk::parallel_machine_world();
  if(stk::parallel_machine_size(comm) != 1) { GTEST_SKIP(); }

  stk::unit_test_util::Args args({"stk_transfer", "--help"});

  stk::set_outputP0(&stk::outputNull(), comm);
  SimpleTransferMainParser helpParser(comm);

  helpParser.parse_command_line_options(args.argc(), args.argv());

  EXPECT_EQ(helpParser.get_parser_status(), stk::transfer_util::TransferParserStatus::PARSE_ONLY);
}

TEST(SimpleTransferMainParser, version)
{
  stk::ParallelMachine comm = stk::parallel_machine_world();
  if(stk::parallel_machine_size(comm) != 1) { GTEST_SKIP(); }

  stk::unit_test_util::Args args({"stk_transfer", "--version"});

  stk::set_outputP0(&stk::outputNull(), comm);
  SimpleTransferMainParser versionParser(comm);

  versionParser.parse_command_line_options(args.argc(), args.argv());

  EXPECT_EQ(versionParser.get_parser_status(), stk::transfer_util::TransferParserStatus::PARSE_ONLY);
}

TEST(TransferMainParser, empty)
{
  stk::ParallelMachine comm = stk::parallel_machine_world();
  if(stk::parallel_machine_size(comm) != 1) { GTEST_SKIP(); }

  stk::unit_test_util::Args args({"stk_transfer"});

  stk::set_outputP0(&stk::outputNull(), comm);
  stk::transfer_util::TransferMainParser emptyParser(comm);


  testing::internal::CaptureStderr();
  emptyParser.parse_command_line_options(args.argc(), args.argv());
  testing::internal::GetCapturedStderr();

  EXPECT_EQ(emptyParser.get_parser_status(), stk::transfer_util::TransferParserStatus::PARSE_ERROR);

  stk::transfer_util::TransferMainSettings emptySettings = emptyParser.generate_transfer_settings();
  EXPECT_EQ(emptySettings.get_num_input_processors(), 1u);
  EXPECT_EQ(emptySettings.get_num_output_processors(), 1u);
  EXPECT_EQ(emptySettings.get_sendMesh_filename(), "");
  EXPECT_EQ(emptySettings.get_recvMesh_filename(), "");
  EXPECT_EQ(emptySettings.get_transfer_fields().size(), 0u);
  EXPECT_EQ(emptySettings.get_extrapolate_option_string(), "IGNORE");
  EXPECT_EQ(emptySettings.get_master_elements_name(), "");
  EXPECT_EQ(emptySettings.get_recv_type_string(), "NODE");
  EXPECT_EQ(emptySettings.get_transfer_type(), "INTERP");
}

TEST(TransferMainParser, basic)
{
  stk::ParallelMachine comm = stk::parallel_machine_world();
  if(stk::parallel_machine_size(comm) != 1) { GTEST_SKIP(); }

  stk::unit_test_util::Args args({"stk_transfer", "--send-mesh", "source.exo", "--recv-mesh", "target.exo"});

  stk::set_outputP0(&stk::outputNull(), comm);
  stk::transfer_util::TransferMainParser basicParser(comm);

  basicParser.parse_command_line_options(args.argc(), args.argv());
  EXPECT_EQ(basicParser.get_parser_status(), stk::transfer_util::TransferParserStatus::SUCCESS);

  stk::transfer_util::TransferMainSettings basicSettings = basicParser.generate_transfer_settings();
  EXPECT_EQ(basicSettings.get_num_input_processors(), 1u);
  EXPECT_EQ(basicSettings.get_num_output_processors(), 1u);
  EXPECT_EQ(basicSettings.get_sendMesh_filename(), "source.exo");
  EXPECT_EQ(basicSettings.get_recvMesh_filename(), "target.exo");
  EXPECT_EQ(basicSettings.get_outputMesh_filename(), "transferred_target.exo");
  EXPECT_EQ(basicSettings.get_transfer_fields().size(), 0u);
  EXPECT_EQ(basicSettings.get_extrapolate_option_string(), "IGNORE");
  EXPECT_EQ(basicSettings.get_master_elements_name(), stk::transfer_util::TransferMainParser::default_master_element_name);
  EXPECT_EQ(basicSettings.get_recv_type_string(), "NODE");
  EXPECT_EQ(basicSettings.get_transfer_type(), "INTERP");
}

TEST(TransferMainParser, basicAll)
{
  stk::ParallelMachine comm = stk::parallel_machine_world();
  if(stk::parallel_machine_size(comm) != 1) { GTEST_SKIP(); }

  stk::unit_test_util::Args args({"stk_transfer", "--send-mesh", "source.exo", "--recv-mesh", "target.exo",
                                  "--output-mesh", "output.exo",
                                  "--field-list", "field_1",
                                  "--part-list", "send_part1,send_part2:recv_part1",
                                  "--extrapolate-option", "EXTRAPOLATE",
                                  "--use-master-elements", "SOME_NAME",
                                  "--recv-type", "EDGE_CENTROID",
                                  "--xfer-type", "PATCH", "--time-steps", "last", "--2d-to-3d-z-coord", "z=2*y"});

  stk::set_outputP0(&stk::outputNull(), comm);
  stk::transfer_util::TransferMainParser basicAllParser(comm);
  basicAllParser.parse_command_line_options(args.argc(), args.argv());
  EXPECT_EQ(basicAllParser.get_parser_status(), stk::transfer_util::TransferParserStatus::SUCCESS);

  stk::transfer_util::TransferMainSettings basicAllSettings = basicAllParser.generate_transfer_settings();
  EXPECT_EQ(basicAllSettings.get_num_input_processors(), 1u);
  EXPECT_EQ(basicAllSettings.get_num_output_processors(), 1u);
  EXPECT_EQ(basicAllSettings.get_sendMesh_filename(), "source.exo");
  EXPECT_EQ(basicAllSettings.get_recvMesh_filename(), "target.exo");
  EXPECT_EQ(basicAllSettings.get_outputMesh_filename(), "output.exo");
  EXPECT_EQ(basicAllSettings.get_transfer_fields().size(), 1u);
  EXPECT_EQ(basicAllSettings.get_field_list_string(), "field_1:field_1");
  EXPECT_EQ(basicAllSettings.get_part_list_string(), "send_part1, send_part2:recv_part1");
  EXPECT_EQ(basicAllSettings.get_extrapolate_option_string(), "EXTRAPOLATE");
  EXPECT_EQ(basicAllSettings.get_master_elements_name(), "SOME_NAME");
  EXPECT_EQ(basicAllSettings.get_recv_type_string(), "EDGE_CENTROID");
  EXPECT_EQ(basicAllSettings.get_transfer_type(), "PATCH");
  EXPECT_EQ(basicAllSettings.get_time_steps_spec(), "LAST");
  EXPECT_EQ(basicAllSettings.get_coord_transf_z_expr(), "z=2*y");

}

TEST(TransferMainParser, basicSendRecvFields)
{
  stk::ParallelMachine comm = stk::parallel_machine_world();
  if(stk::parallel_machine_size(comm) != 1) { GTEST_SKIP(); }

  stk::unit_test_util::Args args({"stk_transfer", "--send-mesh", "source.exo", "--recv-mesh", "target.exo",
                                  "--field-list", "send_field_1:recv_field_1,send_field_2"});

  stk::set_outputP0(&stk::outputNull(), comm);
  stk::transfer_util::TransferMainParser basicAllParser(comm);
  basicAllParser.parse_command_line_options(args.argc(), args.argv());
  EXPECT_EQ(basicAllParser.get_parser_status(), stk::transfer_util::TransferParserStatus::SUCCESS);

  stk::transfer_util::TransferMainSettings basicAllSettings = basicAllParser.generate_transfer_settings();
  EXPECT_EQ(basicAllSettings.get_num_input_processors(), 1u);
  EXPECT_EQ(basicAllSettings.get_num_output_processors(), 1u);
  EXPECT_EQ(basicAllSettings.get_sendMesh_filename(), "source.exo");
  EXPECT_EQ(basicAllSettings.get_recvMesh_filename(), "target.exo");
  EXPECT_EQ(basicAllSettings.get_transfer_fields().size(), 2u);
  EXPECT_EQ(basicAllSettings.get_field_list_string(), "send_field_1:recv_field_1, send_field_2:send_field_2");
  EXPECT_EQ(basicAllSettings.get_extrapolate_option_string(), "IGNORE");

}

TEST(TransferMainParser, basicAllShort)
{
  stk::ParallelMachine comm = stk::parallel_machine_world();
  if(stk::parallel_machine_size(comm) != 1) { GTEST_SKIP(); }

  stk::unit_test_util::Args args({"stk_transfer", "-s", "source.exo", "-r", "target.exo", "-o", "output.exo",
                                  "-f", "field_a,field_b", "-e", "TRUNCATE", "-t", "FACE_GAUSS_POINT"});

  stk::set_outputP0(&stk::outputNull(), comm);
  stk::transfer_util::TransferMainParser basicAllShortParser(comm);
  basicAllShortParser.parse_command_line_options(args.argc(), args.argv());
  EXPECT_EQ(basicAllShortParser.get_parser_status(), stk::transfer_util::TransferParserStatus::SUCCESS);

  stk::transfer_util::TransferMainSettings basicAllShortSettings = basicAllShortParser.generate_transfer_settings();
  EXPECT_EQ(basicAllShortSettings.get_num_input_processors(), 1u);
  EXPECT_EQ(basicAllShortSettings.get_num_output_processors(), 1u);
  EXPECT_EQ(basicAllShortSettings.get_sendMesh_filename(), "source.exo");
  EXPECT_EQ(basicAllShortSettings.get_recvMesh_filename(), "target.exo");
  EXPECT_EQ(basicAllShortSettings.get_outputMesh_filename(), "output.exo");
  EXPECT_EQ(basicAllShortSettings.get_transfer_fields().size(), 2u);
  EXPECT_EQ(basicAllShortSettings.get_field_list_string(), "field_a:field_a, field_b:field_b");
  EXPECT_EQ(basicAllShortSettings.get_extrapolate_option_string(), "TRUNCATE");
  EXPECT_EQ(basicAllShortSettings.get_recv_type_string(), "FACE_GAUSS_POINT");

}

TEST(TransferMainParser, basicAllShortCaseInsensitivity)
{
  stk::ParallelMachine comm = stk::parallel_machine_world();
  if(stk::parallel_machine_size(comm) != 1) { GTEST_SKIP(); }

  stk::unit_test_util::Args args({"stk_transfer", "-s", "Source.exo", "-r", "targeT.exo", "-o", "outPut.exo",
                                  "-f", "FIELD_A,field_b", "-x", "CoPy", "-e", "truncate", "-t", "FACE_gauss_POINT",
                                  "-u", "InTrePid2"});

  stk::set_outputP0(&stk::outputNull(), comm);
  stk::transfer_util::TransferMainParser basicAllShortParser(comm);
  basicAllShortParser.parse_command_line_options(args.argc(), args.argv());
  EXPECT_EQ(basicAllShortParser.get_parser_status(), stk::transfer_util::TransferParserStatus::SUCCESS);

  stk::transfer_util::TransferMainSettings basicAllShortSettings = basicAllShortParser.generate_transfer_settings();
  EXPECT_EQ(basicAllShortSettings.get_num_input_processors(), 1u);
  EXPECT_EQ(basicAllShortSettings.get_num_output_processors(), 1u);
  EXPECT_EQ(basicAllShortSettings.get_sendMesh_filename(), "Source.exo");
  EXPECT_EQ(basicAllShortSettings.get_recvMesh_filename(), "targeT.exo");
  EXPECT_EQ(basicAllShortSettings.get_outputMesh_filename(), "outPut.exo");
  EXPECT_EQ(basicAllShortSettings.get_transfer_fields().size(), 2u);
  EXPECT_EQ(basicAllShortSettings.get_field_list_string(), "FIELD_A:FIELD_A, field_b:field_b");
  EXPECT_EQ(basicAllShortSettings.get_transfer_type(), "COPY");
  EXPECT_EQ(basicAllShortSettings.get_extrapolate_option_string(), "TRUNCATE");
  EXPECT_EQ(basicAllShortSettings.get_recv_type_string(), "FACE_GAUSS_POINT");
  EXPECT_EQ(basicAllShortSettings.get_master_elements_name(), "INTREPID2");

}

TEST(TransferMainParser, invalidFieldList)
{
  stk::ParallelMachine comm = stk::parallel_machine_world();
  if(stk::parallel_machine_size(comm) != 1) { GTEST_SKIP(); }

  stk::unit_test_util::Args args({"stk_transfer", "--send-mesh", "source.exo", "--recv-mesh", "target.exo",
                                  "--field-list", "send_field_1:recv_field_1:send_field_2"});

  stk::set_outputP0(&stk::outputNull(), comm);
  stk::transfer_util::TransferMainParser basicAllParser(comm);
  basicAllParser.parse_command_line_options(args.argc(), args.argv());
  EXPECT_EQ(basicAllParser.get_parser_status(), stk::transfer_util::TransferParserStatus::SUCCESS);

  EXPECT_ANY_THROW(basicAllParser.generate_transfer_settings());

}

TEST(TransferMainParser, invalidRecvFieldOnly)
{
  stk::ParallelMachine comm = stk::parallel_machine_world();
  if(stk::parallel_machine_size(comm) != 1) { GTEST_SKIP(); }

  stk::unit_test_util::Args args({"stk_transfer", "--send-mesh", "source.exo", "--recv-mesh", "target.exo",
                                  "--field-list", "send_field_1:recv_field_1,:field_2"});

  stk::set_outputP0(&stk::outputNull(), comm);
  stk::transfer_util::TransferMainParser basicAllParser(comm);
  basicAllParser.parse_command_line_options(args.argc(), args.argv());
  EXPECT_EQ(basicAllParser.get_parser_status(), stk::transfer_util::TransferParserStatus::SUCCESS);

  EXPECT_ANY_THROW(basicAllParser.generate_transfer_settings());
}

TEST(TransferMainParser, emptySendPartList)
{
  stk::ParallelMachine comm = stk::parallel_machine_world();
  if(stk::parallel_machine_size(comm) != 1) { GTEST_SKIP(); }

  stk::unit_test_util::Args args({"stk_transfer", "--send-mesh", "source.exo", "--recv-mesh", "target.exo",
                                  "--part-list", ":recv_part_1"});

  stk::set_outputP0(&stk::outputNull(), comm);
  stk::transfer_util::TransferMainParser basicAllParser(comm);
  basicAllParser.parse_command_line_options(args.argc(), args.argv());
  EXPECT_EQ(basicAllParser.get_parser_status(), stk::transfer_util::TransferParserStatus::SUCCESS);

  stk::transfer_util::TransferMainSettings basicAllSettings = basicAllParser.generate_transfer_settings();
  EXPECT_EQ(basicAllSettings.get_transfer_send_parts().size(), 0u);
  EXPECT_EQ(basicAllSettings.get_transfer_recv_parts().size(), 1u);
  EXPECT_EQ(basicAllSettings.get_transfer_recv_parts()[0], "recv_part_1");

}

TEST(TransferMainParser, emptySendPartMultiRecvPartsList)
{
  stk::ParallelMachine comm = stk::parallel_machine_world();
  if(stk::parallel_machine_size(comm) != 1) { GTEST_SKIP(); }

  stk::unit_test_util::Args args({"stk_transfer", "--send-mesh", "source.exo", "--recv-mesh", "target.exo",
                                  "--part-list", ":recv_part_1,recv_part_2"});

  stk::set_outputP0(&stk::outputNull(), comm);
  stk::transfer_util::TransferMainParser basicAllParser(comm);
  basicAllParser.parse_command_line_options(args.argc(), args.argv());
  EXPECT_EQ(basicAllParser.get_parser_status(), stk::transfer_util::TransferParserStatus::SUCCESS);

  stk::transfer_util::TransferMainSettings basicAllSettings = basicAllParser.generate_transfer_settings();
  EXPECT_EQ(basicAllSettings.get_transfer_send_parts().size(), 0u);
  EXPECT_EQ(basicAllSettings.get_transfer_recv_parts().size(), 2u);
  EXPECT_EQ(basicAllSettings.get_transfer_recv_parts()[0], "recv_part_1");
  EXPECT_EQ(basicAllSettings.get_transfer_recv_parts()[1], "recv_part_2");
  EXPECT_EQ(basicAllSettings.get_recv_part_list_string(), "recv_part_1, recv_part_2");

}

TEST(TransferMainParser, emptyRecvPartList)
{
  stk::ParallelMachine comm = stk::parallel_machine_world();
  if(stk::parallel_machine_size(comm) != 1) { GTEST_SKIP(); }

  stk::unit_test_util::Args args({"stk_transfer", "--send-mesh", "source.exo", "--recv-mesh", "target.exo",
                                  "--part-list", "send_part_1"});

  stk::set_outputP0(&stk::outputNull(), comm);
  stk::transfer_util::TransferMainParser basicAllParser(comm);
  basicAllParser.parse_command_line_options(args.argc(), args.argv());
  EXPECT_EQ(basicAllParser.get_parser_status(), stk::transfer_util::TransferParserStatus::SUCCESS);

  stk::transfer_util::TransferMainSettings basicAllSettings = basicAllParser.generate_transfer_settings();
  EXPECT_EQ(basicAllSettings.get_transfer_send_parts().size(), 1u);
  EXPECT_EQ(basicAllSettings.get_transfer_send_parts()[0], "send_part_1");
  EXPECT_EQ(basicAllSettings.get_transfer_recv_parts().size(), 0u);

}

TEST(TransferMainParser, emptyRecvPartMultiSendPartsList)
{
  stk::ParallelMachine comm = stk::parallel_machine_world();
  if(stk::parallel_machine_size(comm) != 1) { GTEST_SKIP(); }

  stk::unit_test_util::Args args({"stk_transfer", "--send-mesh", "source.exo", "--recv-mesh", "target.exo",
                                  "--part-list", "send_part_1,send_part_2"});

  stk::set_outputP0(&stk::outputNull(), comm);
  stk::transfer_util::TransferMainParser basicAllParser(comm);
  basicAllParser.parse_command_line_options(args.argc(), args.argv());
  EXPECT_EQ(basicAllParser.get_parser_status(), stk::transfer_util::TransferParserStatus::SUCCESS);

  stk::transfer_util::TransferMainSettings basicAllSettings = basicAllParser.generate_transfer_settings();
  EXPECT_EQ(basicAllSettings.get_transfer_send_parts().size(), 2u);
  EXPECT_EQ(basicAllSettings.get_transfer_send_parts()[0], "send_part_1");
  EXPECT_EQ(basicAllSettings.get_transfer_send_parts()[1], "send_part_2");
  EXPECT_EQ(basicAllSettings.get_transfer_recv_parts().size(), 0u);
  EXPECT_EQ(basicAllSettings.get_send_part_list_string(), "send_part_1, send_part_2");

}

TEST(TransferMainParser, invalidPartList)
{
  stk::ParallelMachine comm = stk::parallel_machine_world();
  if(stk::parallel_machine_size(comm) != 1) { GTEST_SKIP(); }

  stk::unit_test_util::Args args({"stk_transfer", "--send-mesh", "source.exo", "--recv-mesh", "target.exo",
                                  "--part-list", "send_part_1:recv_part_1:send_part_2"});

  stk::set_outputP0(&stk::outputNull(), comm);
  stk::transfer_util::TransferMainParser basicAllParser(comm);
  basicAllParser.parse_command_line_options(args.argc(), args.argv());
  EXPECT_EQ(basicAllParser.get_parser_status(), stk::transfer_util::TransferParserStatus::SUCCESS);

  EXPECT_ANY_THROW(basicAllParser.generate_transfer_settings());

}

TEST(TransferMainParser, invalidArg)
{
  stk::ParallelMachine comm = stk::parallel_machine_world();
  if(stk::parallel_machine_size(comm) != 1) { GTEST_SKIP(); }

  stk::unit_test_util::Args args({"stk_transfer", "--send-mesh", "source.exo", "--recv-mesh", "target.exo",
                                  "--field-vector", "field_1", "--extrapolate-option", "EXTRAPOLATE"});

  stk::set_outputP0(&stk::outputNull(), comm);
  stk::transfer_util::TransferMainParser invalidArgParser(comm);

  testing::internal::CaptureStderr();
  invalidArgParser.parse_command_line_options(args.argc(), args.argv());
  testing::internal::GetCapturedStderr();

  EXPECT_EQ(invalidArgParser.get_parser_status(), stk::transfer_util::TransferParserStatus::PARSE_ERROR);

  stk::transfer_util::TransferMainSettings invalidArgSettings = invalidArgParser.generate_transfer_settings();
  EXPECT_EQ(invalidArgSettings.get_num_input_processors(), 1u);
  EXPECT_EQ(invalidArgSettings.get_num_output_processors(), 1u);
  EXPECT_EQ(invalidArgSettings.get_sendMesh_filename(), "");
  EXPECT_EQ(invalidArgSettings.get_recvMesh_filename(), "");
  EXPECT_EQ(invalidArgSettings.get_transfer_fields().size(), 0u);
  EXPECT_EQ(invalidArgSettings.get_extrapolate_option_string(), "IGNORE");
}

TEST(TransferMainParser, invalidArgShort)
{
  stk::ParallelMachine comm = stk::parallel_machine_world();
  if(stk::parallel_machine_size(comm) != 1) { GTEST_SKIP(); }

  stk::unit_test_util::Args args({"stk_transfer", "-f", "source.exo", "-t", "target.exo",
                                  "-l", "field_a", "-z", "TRUNCATE"});

  stk::set_outputP0(&stk::outputNull(), comm);
  stk::transfer_util::TransferMainParser invalidArgShortParser(comm);

  testing::internal::CaptureStderr();
  invalidArgShortParser.parse_command_line_options(args.argc(), args.argv());
  testing::internal::GetCapturedStderr();

  EXPECT_EQ(invalidArgShortParser.get_parser_status(), stk::transfer_util::TransferParserStatus::PARSE_ERROR);

  stk::transfer_util::TransferMainSettings invalidArgShortSettings = invalidArgShortParser.generate_transfer_settings();
  EXPECT_EQ(invalidArgShortSettings.get_num_input_processors(), 1u);
  EXPECT_EQ(invalidArgShortSettings.get_num_output_processors(), 1u);
  EXPECT_EQ(invalidArgShortSettings.get_sendMesh_filename(), "");
  EXPECT_EQ(invalidArgShortSettings.get_recvMesh_filename(), "");
  EXPECT_EQ(invalidArgShortSettings.get_transfer_fields().size(), 0u);
  EXPECT_EQ(invalidArgShortSettings.get_extrapolate_option_string(), "IGNORE");
}

namespace {
class TransferParserTest : public ::testing::Test
{
  protected:
    TransferParserTest() :
      basicParser(stk::parallel_machine_world())
    {}

    void SetUp()
    {
      stk::ParallelMachine comm = stk::parallel_machine_world();
      if(stk::parallel_machine_size(comm) != 1)
      {
        GTEST_SKIP();
      }
      stk::set_outputP0(&stk::outputNull(), comm);
    }

    virtual ~TransferParserTest()
    {
      stk::reset_default_output_streams(stk::parallel_machine_world());
    }

    void setArgs(stk::unit_test_util::Args args)
    {
      basicParser.parse_command_line_options(args.argc(), args.argv());
      EXPECT_EQ(basicParser.get_parser_status(), stk::transfer_util::TransferParserStatus::SUCCESS);
      basicSettings =  basicParser.generate_transfer_settings();
    }

    stk::transfer_util::TransferMainParser basicParser;
    stk::transfer_util::TransferMainSettings basicSettings;
};

}

std::array<double, 3> make_array(double x, double y, double z)
{
  return {x, y, z};
}

double pi()
{
  return 2*std::atan2(1, 0);
}

}


TEST_F(TransferParserTest, TwoDTo3DAxisymmetricNotSpecified)
{
  setArgs({{"stk_transfer", "--send-mesh", "source.exo", "--recv-mesh", "target.exo"}});

  EXPECT_FALSE(basicSettings.get_enable_2d_to_3d_axisymmetric_transfer());
}


TEST_F(TransferParserTest, TwoDTo3DAxisymmetricDefaults)
{
  setArgs({{"stk_transfer", "--send-mesh", "source.exo", "--recv-mesh", "target.exo", "--2d-to-3d-axisymmetric"}});

  EXPECT_TRUE(basicSettings.get_enable_2d_to_3d_axisymmetric_transfer());
  stk::transfer_util::TwoDTo3DAxisymmetricParams params = basicSettings.get_2d_to_3d_axisymmetric_transfer_params();
  EXPECT_EQ(params.axis, make_array(0, 0, 1));
  EXPECT_EQ(params.trans, make_array(0, 0, 0));
}

TEST_F(TransferParserTest, TwoDTo3DAxisymmetricShortOption)
{
  setArgs({{"stk_transfer", "--send-mesh", "source.exo", "--recv-mesh", "target.exo", "-2d3daxi"}});

  EXPECT_TRUE(basicSettings.get_enable_2d_to_3d_axisymmetric_transfer());
  stk::transfer_util::TwoDTo3DAxisymmetricParams params = basicSettings.get_2d_to_3d_axisymmetric_transfer_params();
  EXPECT_EQ(params.axis, make_array(0, 0, 1));
  EXPECT_EQ(params.trans, make_array(0, 0, 0));
}

TEST_F(TransferParserTest, TwoDTo3DAxisymmetricAxisSpecified)
{
  setArgs({{"stk_transfer", "--send-mesh", "source.exo", "--recv-mesh", "target.exo", "--2d-to-3d-axisymmetric", "1,2,3"}});

  EXPECT_TRUE(basicSettings.get_enable_2d_to_3d_axisymmetric_transfer());
  stk::transfer_util::TwoDTo3DAxisymmetricParams params = basicSettings.get_2d_to_3d_axisymmetric_transfer_params();
  EXPECT_EQ(params.axis, make_array(1, 2, 3));
  EXPECT_EQ(params.trans, make_array(0, 0, 0));
}

TEST_F(TransferParserTest, TwoDTo3DAxisymmetricAxisAndTranslationSpecified)
{
  setArgs({{"stk_transfer", "--send-mesh", "source.exo", "--recv-mesh", "target.exo", "--2d-to-3d-axisymmetric", "1,2,3,4,5,6"}});

  EXPECT_TRUE(basicSettings.get_enable_2d_to_3d_axisymmetric_transfer());
  stk::transfer_util::TwoDTo3DAxisymmetricParams params = basicSettings.get_2d_to_3d_axisymmetric_transfer_params();
  EXPECT_EQ(params.axis, make_array(1, 2, 3));
  EXPECT_EQ(params.trans, make_array(4, 5, 6));
}

TEST_F(TransferParserTest, IncorrectNumParams)
{
  EXPECT_ANY_THROW(setArgs({{"stk_transfer", "--send-mesh", "source.exo", "--recv-mesh", "target.exo", "--2d-to-3d-axisymmetric", "1,2"}}));
}

TEST_F(TransferParserTest, IncorrectNumParamsWithTrailingComma)
{
  EXPECT_ANY_THROW(setArgs({{"stk_transfer", "--send-mesh", "source.exo", "--recv-mesh", "target.exo", "--2d-to-3d-axisymmetric", "1,2,"}}));
}

TEST_F(TransferParserTest, ThreeDTo3DAxisymmetricNotSpecified)
{
  setArgs({{"stk_transfer", "--send-mesh", "source.exo", "--recv-mesh", "target.exo"}});

  EXPECT_FALSE(basicSettings.get_enable_3d_to_3d_axisymmetric_transfer());
}

TEST_F(TransferParserTest, ThreeDTo3DAxisymmetricThetaOnly)
{
  setArgs({{"stk_transfer", "--send-mesh", "source.exo", "--recv-mesh", "target.exo", "--3d-to-3d-axisymmetric", "30,60"}});

  EXPECT_TRUE(basicSettings.get_enable_3d_to_3d_axisymmetric_transfer());
  stk::transfer_util::ThreeDTo3DAxisymmetricParams params = basicSettings.get_3d_to_3d_axisymmetric_transfer_params();
  EXPECT_NEAR(params.theta_min, pi()/6, 1e-13);
  EXPECT_NEAR(params.theta_max, pi()/3, 1e-13);
  EXPECT_EQ(params.axis, make_array(0, 0, 1));
  EXPECT_EQ(params.trans, make_array(0, 0, 0));
}

TEST_F(TransferParserTest, ThreeDTo3DAxisymmetricShortArg)
{
  setArgs({{"stk_transfer", "--send-mesh", "source.exo", "--recv-mesh", "target.exo", "-3d3daxi", "30,60"}});

  EXPECT_TRUE(basicSettings.get_enable_3d_to_3d_axisymmetric_transfer());
  stk::transfer_util::ThreeDTo3DAxisymmetricParams params = basicSettings.get_3d_to_3d_axisymmetric_transfer_params();
  EXPECT_NEAR(params.theta_min, pi()/6, 1e-13);
  EXPECT_NEAR(params.theta_max, pi()/3, 1e-13);
  EXPECT_EQ(params.axis, make_array(0, 0, 1));
  EXPECT_EQ(params.trans, make_array(0, 0, 0));
}

TEST_F(TransferParserTest, ThreeDTo3DAxisymmetricThetaAndAxis)
{
  setArgs({{"stk_transfer", "--send-mesh", "source.exo", "--recv-mesh", "target.exo", "--3d-to-3d-axisymmetric", "30,60,1,2,3"}});

  EXPECT_TRUE(basicSettings.get_enable_3d_to_3d_axisymmetric_transfer());
  stk::transfer_util::ThreeDTo3DAxisymmetricParams params = basicSettings.get_3d_to_3d_axisymmetric_transfer_params();
  EXPECT_NEAR(params.theta_min, pi()/6, 1e-13);
  EXPECT_NEAR(params.theta_max, pi()/3, 1e-13);
  EXPECT_EQ(params.axis, make_array(1, 2, 3));
  EXPECT_EQ(params.trans, make_array(0, 0, 0));
}

TEST_F(TransferParserTest, ThreeDTo3DAxisymmetricThetaAndAxisAndTranslation)
{
  setArgs({{"stk_transfer", "--send-mesh", "source.exo", "--recv-mesh", "target.exo", "--3d-to-3d-axisymmetric", "30,60,1,2,3,4,5,6"}});

  EXPECT_TRUE(basicSettings.get_enable_3d_to_3d_axisymmetric_transfer());
  stk::transfer_util::ThreeDTo3DAxisymmetricParams params = basicSettings.get_3d_to_3d_axisymmetric_transfer_params();
  EXPECT_NEAR(params.theta_min, pi()/6, 1e-13);
  EXPECT_NEAR(params.theta_max, pi()/3, 1e-13);
  EXPECT_EQ(params.axis, make_array(1, 2, 3));
  EXPECT_EQ(params.trans, make_array(4, 5, 6));
}