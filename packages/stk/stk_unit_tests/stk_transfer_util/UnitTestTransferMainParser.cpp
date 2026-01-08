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

  std::string get_transfer_fields() const
  {
    return m_cmdLineParser.is_option_parsed(m_optionNames.fieldList) ?
      m_cmdLineParser.get_option_value<std::string>(m_optionNames.fieldList) : "defaulting to all fields";
  }

  std::string get_extrapolate_option() const
  {
    return m_cmdLineParser.is_option_parsed(m_optionNames.ExtrapolateOption) ?
      m_cmdLineParser.get_option_value<std::string>(m_optionNames.ExtrapolateOption) : "defaulting to IGNORE";
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
  EXPECT_EQ("defaulting to IGNORE", emptyParser.get_extrapolate_option());

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
  EXPECT_EQ("defaulting to IGNORE", basicParser.get_extrapolate_option());

  basicParser.parse_command_line_options(args.argc(), args.argv());

  EXPECT_EQ(basicParser.get_parser_status(), stk::transfer_util::TransferParserStatus::SUCCESS);
  EXPECT_EQ("source.exo", basicParser.get_sendMesh_filename());
  EXPECT_EQ("target.exo", basicParser.get_recvMesh_filename());
  EXPECT_EQ("defaulting to all fields", basicParser.get_transfer_fields());
  EXPECT_EQ("defaulting to IGNORE", basicParser.get_extrapolate_option());
  EXPECT_EQ(1u, basicParser.get_num_input_processors());
  EXPECT_EQ(1u, basicParser.get_num_output_processors());

}

TEST(SimpleTransferMainParser, basicAll)
{
  stk::ParallelMachine comm = stk::parallel_machine_world();
  if(stk::parallel_machine_size(comm) != 1) { GTEST_SKIP(); }

  stk::unit_test_util::Args args({"stk_transfer", "--send-mesh", "source.exo", "--recv-mesh", "target.exo",
                                  "--field-list", "field_1", "--extrapolate-option", "EXTRAPOLATE"});

  stk::set_outputP0(&stk::outputNull(), comm);
  SimpleTransferMainParser basicAllParser(comm);

  EXPECT_EQ("not parsed", basicAllParser.get_sendMesh_filename());
  EXPECT_EQ("not parsed", basicAllParser.get_recvMesh_filename());
  EXPECT_EQ("defaulting to all fields", basicAllParser.get_transfer_fields());
  EXPECT_EQ("defaulting to IGNORE", basicAllParser.get_extrapolate_option());

  basicAllParser.parse_command_line_options(args.argc(), args.argv());

  EXPECT_EQ(basicAllParser.get_parser_status(), stk::transfer_util::TransferParserStatus::SUCCESS);
  EXPECT_EQ("source.exo", basicAllParser.get_sendMesh_filename());
  EXPECT_EQ("target.exo", basicAllParser.get_recvMesh_filename());
  EXPECT_EQ("field_1", basicAllParser.get_transfer_fields());
  EXPECT_EQ("EXTRAPOLATE", basicAllParser.get_extrapolate_option());
  EXPECT_EQ(1u, basicAllParser.get_num_input_processors());
  EXPECT_EQ(1u, basicAllParser.get_num_output_processors());

}

TEST(SimpleTransferMainParser, basicAllShort)
{
  stk::ParallelMachine comm = stk::parallel_machine_world();
  if(stk::parallel_machine_size(comm) != 1) { GTEST_SKIP(); }

  stk::unit_test_util::Args args({"stk_transfer", "-s", "source.exo", "-r", "target.exo",
                                  "-f", "field_a,field_b", "-e", "TRUNCATE"});

  stk::set_outputP0(&stk::outputNull(), comm);
  SimpleTransferMainParser basicAllShortParser(comm);

  EXPECT_EQ("not parsed", basicAllShortParser.get_sendMesh_filename());
  EXPECT_EQ("not parsed", basicAllShortParser.get_recvMesh_filename());
  EXPECT_EQ("defaulting to all fields", basicAllShortParser.get_transfer_fields());
  EXPECT_EQ("defaulting to IGNORE", basicAllShortParser.get_extrapolate_option());

  basicAllShortParser.parse_command_line_options(args.argc(), args.argv());

  EXPECT_EQ(basicAllShortParser.get_parser_status(), stk::transfer_util::TransferParserStatus::SUCCESS);
  EXPECT_EQ("source.exo", basicAllShortParser.get_sendMesh_filename());
  EXPECT_EQ("target.exo", basicAllShortParser.get_recvMesh_filename());
  EXPECT_EQ("field_a,field_b", basicAllShortParser.get_transfer_fields());
  EXPECT_EQ("TRUNCATE", basicAllShortParser.get_extrapolate_option());
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

  EXPECT_EQ(examplesParser.get_quickExample(), "quick example not provided yet");
  EXPECT_EQ(examplesParser.get_longExamples(), "long examples not provided yet");
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
  EXPECT_EQ(basicSettings.get_transfer_fields().size(), 0u);
  EXPECT_EQ(basicSettings.get_extrapolate_option_string(), "IGNORE");
}

TEST(TransferMainParser, basicAll)
{
  stk::ParallelMachine comm = stk::parallel_machine_world();
  if(stk::parallel_machine_size(comm) != 1) { GTEST_SKIP(); }

  stk::unit_test_util::Args args({"stk_transfer", "--send-mesh", "source.exo", "--recv-mesh", "target.exo",
                                  "--field-list", "field_1", "--extrapolate-option", "EXTRAPOLATE"});

  stk::set_outputP0(&stk::outputNull(), comm);
  stk::transfer_util::TransferMainParser basicAllParser(comm);
  basicAllParser.parse_command_line_options(args.argc(), args.argv());
  EXPECT_EQ(basicAllParser.get_parser_status(), stk::transfer_util::TransferParserStatus::SUCCESS);

  stk::transfer_util::TransferMainSettings basicAllSettings = basicAllParser.generate_transfer_settings();
  EXPECT_EQ(basicAllSettings.get_num_input_processors(), 1u); 
  EXPECT_EQ(basicAllSettings.get_num_output_processors(), 1u); 
  EXPECT_EQ(basicAllSettings.get_sendMesh_filename(), "source.exo");
  EXPECT_EQ(basicAllSettings.get_recvMesh_filename(), "target.exo");
  EXPECT_EQ(basicAllSettings.get_transfer_fields().size(), 1u);
  EXPECT_EQ(basicAllSettings.get_field_list_string(), "field_1:field_1");
  EXPECT_EQ(basicAllSettings.get_extrapolate_option_string(), "EXTRAPOLATE");

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

  stk::unit_test_util::Args args({"stk_transfer", "-s", "source.exo", "-r", "target.exo",
                                  "-f", "field_a,field_b", "-e", "TRUNCATE"});

  stk::set_outputP0(&stk::outputNull(), comm);
  stk::transfer_util::TransferMainParser basicAllShortParser(comm);
  basicAllShortParser.parse_command_line_options(args.argc(), args.argv());
  EXPECT_EQ(basicAllShortParser.get_parser_status(), stk::transfer_util::TransferParserStatus::SUCCESS);
  
  stk::transfer_util::TransferMainSettings basicAllShortSettings = basicAllShortParser.generate_transfer_settings();
  EXPECT_EQ(basicAllShortSettings.get_num_input_processors(), 1u); 
  EXPECT_EQ(basicAllShortSettings.get_num_output_processors(), 1u); 
  EXPECT_EQ(basicAllShortSettings.get_sendMesh_filename(), "source.exo");
  EXPECT_EQ(basicAllShortSettings.get_recvMesh_filename(), "target.exo");
  EXPECT_EQ(basicAllShortSettings.get_transfer_fields().size(), 2u);
  EXPECT_EQ(basicAllShortSettings.get_field_list_string(), "field_a:field_a, field_b:field_b");
  EXPECT_EQ(basicAllShortSettings.get_extrapolate_option_string(), "TRUNCATE");

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

}
