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
#include <algorithm>
#include <array>
#include <cmath>
#include <string>
#include <vector>
#include "stk_transfer_util/TransferMainInputFile.hpp"
#include "stk_transfer_util/TransferMainParser.hpp"
#include "stk_transfer_util/TransferMainSettings.hpp"
#include "stk_unit_test_utils/CommandLineArgs.hpp"
#include "stk_util/parallel/OutputStreams.hpp"
#include "stk_util/parallel/Parallel.hpp"
#include "stk_util/stk_config.h"

namespace {

using stk::transfer_util::TransferParserStatus;

// All tests parse YAML from in-memory strings, via either
// TransferMainParser::set_input_file_contents or parse_transfer_input_file_yaml.
class TransferMainInputFile : public ::testing::Test
{
protected:
  TransferMainInputFile()
    : m_comm(stk::parallel_machine_world()),
      m_parser(m_comm)
  {}

  void SetUp() override
  {
    if (stk::parallel_machine_size(m_comm) != 1) { GTEST_SKIP(); }
    stk::set_outputP0(&stk::outputNull(), m_comm);
  }

  void TearDown() override
  {
    stk::reset_default_output_streams(m_comm);
  }

  // Loads the yaml text as if it were read from --input-file, then parses the given
  // command-line arguments (which take precedence over the yaml options).
  TransferParserStatus parse(const std::string& yamlText, const std::vector<std::string>& extraArgs = {})
  {
    m_parser.set_input_file_contents(yamlText);

    std::vector<std::string> argStrings = {"stk_transfer"};
    argStrings.insert(argStrings.end(), extraArgs.begin(), extraArgs.end());
    stk::unit_test_util::Args args(argStrings);

    testing::internal::CaptureStderr();
    m_parser.parse_command_line_options(args.argc(), args.argv());
    m_stderr = testing::internal::GetCapturedStderr();

    return m_parser.get_parser_status();
  }

  stk::transfer_util::TransferMainSettings parse_and_generate(const std::string& yamlText,
                                                              const std::vector<std::string>& extraArgs = {})
  {
    EXPECT_EQ(TransferParserStatus::SUCCESS, parse(yamlText, extraArgs)) << m_stderr;
    return m_parser.generate_transfer_settings();
  }

  void expect_yaml_error(const std::string& yamlText, const std::string& expectedMessageSubstring)
  {
    try {
      stk::transfer_util::parse_transfer_input_file_yaml(yamlText);
      FAIL() << "Expected an exception for input:\n" << yamlText;
    } catch (const std::exception& e) {
      EXPECT_NE(std::string(e.what()).find(expectedMessageSubstring), std::string::npos)
        << "Message '" << e.what() << "' does not contain '" << expectedMessageSubstring << "'";
    }
  }

  static std::array<double, 3> make_array(double x, double y, double z) { return {x, y, z}; }
  static double pi() { return 2*std::atan2(1, 0); }

  stk::ParallelMachine m_comm;
  stk::transfer_util::TransferMainParser m_parser;
  std::string m_stderr;
};

#ifdef STK_HAS_YAML_CPP

TEST_F(TransferMainInputFile, isSupported)
{
  EXPECT_TRUE(stk::transfer_util::transfer_input_file_supported());
}

TEST_F(TransferMainInputFile, sendRecvOutputMesh)
{
  const std::string yaml = R"(
send-mesh: sendMesh.exo
recv-mesh: receiveMesh.exo
output-mesh: output.exo
)";
  stk::transfer_util::TransferMainSettings settings = parse_and_generate(yaml);

  EXPECT_EQ("sendMesh.exo", settings.get_sendMesh_filename());
  EXPECT_EQ("receiveMesh.exo", settings.get_recvMesh_filename());
  EXPECT_EQ("output.exo", settings.get_outputMesh_filename());
  EXPECT_EQ(1u, settings.get_num_input_processors());
  EXPECT_EQ(1u, settings.get_num_output_processors());
}

TEST_F(TransferMainInputFile, defaultsMatchCommandLineDefaults)
{
  const std::string yaml = R"(
send-mesh: source.exo
recv-mesh: target.exo
)";
  stk::transfer_util::TransferMainSettings settings = parse_and_generate(yaml);

  EXPECT_EQ("source.exo", settings.get_sendMesh_filename());
  EXPECT_EQ("target.exo", settings.get_recvMesh_filename());
  EXPECT_EQ("transferred_target.exo", settings.get_outputMesh_filename());
  EXPECT_EQ(0u, settings.get_transfer_fields().size());
  EXPECT_EQ(0u, settings.get_transfer_send_parts().size());
  EXPECT_EQ(0u, settings.get_transfer_recv_parts().size());
  EXPECT_EQ("IGNORE", settings.get_extrapolate_option_string());
  EXPECT_EQ(stk::transfer_util::TransferMainParser::default_master_element_name, settings.get_master_elements_name());
  EXPECT_EQ("NODE", settings.get_recv_type_string());
  EXPECT_EQ("INTERP", settings.get_transfer_type());
  EXPECT_EQ("ALL", settings.get_time_steps_spec());
  EXPECT_EQ("ZPLANE", settings.get_2d_to_3d_mapping_type());
  EXPECT_EQ("z=0", settings.get_coord_transf_z_expr());
  EXPECT_FALSE(settings.get_enable_2d_to_3d_axisymmetric_transfer());
  EXPECT_FALSE(settings.get_enable_3d_to_3d_axisymmetric_transfer());
}

TEST_F(TransferMainInputFile, allScalarOptions)
{
  const std::string yaml = R"(
send-mesh: source.exo
recv-mesh: target.exo
output-mesh: output.exo
field-list:
  - field_1
part-list:
  send: [send_part1, send_part2]
  recv: [recv_part1]
extrapolate-option: EXTRAPOLATE
use-master-elements: SOME_NAME
recv-type: EDGE_CENTROID
xfer-type: PATCH
time-steps: last
2d-to-3d-mapping-type: EXTRUDE
2d-to-3d-z-coord: z=2*y
)";
  stk::transfer_util::TransferMainSettings settings = parse_and_generate(yaml);

  EXPECT_EQ("source.exo", settings.get_sendMesh_filename());
  EXPECT_EQ("target.exo", settings.get_recvMesh_filename());
  EXPECT_EQ("output.exo", settings.get_outputMesh_filename());
  EXPECT_EQ(1u, settings.get_transfer_fields().size());
  EXPECT_EQ("field_1:field_1", settings.get_field_list_string());
  EXPECT_EQ("send_part1, send_part2:recv_part1", settings.get_part_list_string());
  EXPECT_EQ("EXTRAPOLATE", settings.get_extrapolate_option_string());
  EXPECT_EQ("SOME_NAME", settings.get_master_elements_name());
  EXPECT_EQ("EDGE_CENTROID", settings.get_recv_type_string());
  EXPECT_EQ("PATCH", settings.get_transfer_type());
  EXPECT_EQ("LAST", settings.get_time_steps_spec());
  EXPECT_EQ("EXTRUDE", settings.get_2d_to_3d_mapping_type());
  EXPECT_EQ("z=2*y", settings.get_coord_transf_z_expr());
}

TEST_F(TransferMainInputFile, caseInsensitivity)
{
  const std::string yaml = R"(
send-mesh: Source.exo
recv-mesh: targeT.exo
output-mesh: outPut.exo
field-list: [FIELD_A, field_b]
xfer-type: CoPy
extrapolate-option: truncate
recv-type: FACE_gauss_POINT
use-master-elements: InTrePid2
)";
  stk::transfer_util::TransferMainSettings settings = parse_and_generate(yaml);

  EXPECT_EQ("Source.exo", settings.get_sendMesh_filename());
  EXPECT_EQ("targeT.exo", settings.get_recvMesh_filename());
  EXPECT_EQ("outPut.exo", settings.get_outputMesh_filename());
  EXPECT_EQ("FIELD_A:FIELD_A, field_b:field_b", settings.get_field_list_string());
  EXPECT_EQ("COPY", settings.get_transfer_type());
  EXPECT_EQ("TRUNCATE", settings.get_extrapolate_option_string());
  EXPECT_EQ("FACE_GAUSS_POINT", settings.get_recv_type_string());
  EXPECT_EQ("INTREPID2", settings.get_master_elements_name());
}

TEST_F(TransferMainInputFile, quotedValuesAndComments)
{
  const std::string yaml = R"(
# stk_transfer input file
send-mesh: "send mesh.exo"   # names may contain spaces when quoted
recv-mesh: 'recv.exo'
2d-to-3d-z-coord: "z = x + y"
)";
  stk::transfer_util::TransferMainSettings settings = parse_and_generate(yaml);

  EXPECT_EQ("send mesh.exo", settings.get_sendMesh_filename());
  EXPECT_EQ("recv.exo", settings.get_recvMesh_filename());
  EXPECT_EQ("z = x + y", settings.get_coord_transf_z_expr());
}

TEST_F(TransferMainInputFile, fieldListSendRecvMaps)
{
  const std::string yaml = R"(
send-mesh: source.exo
recv-mesh: target.exo
field-list:
  - {send: send_field_1, recv: recv_field_1}
  - send_field_2
  - send: send_field_3
    recv: recv_field_3
  - {send: send_field_4}
)";
  stk::transfer_util::TransferMainSettings settings = parse_and_generate(yaml);

  ASSERT_EQ(4u, settings.get_transfer_fields().size());
  EXPECT_EQ("send_field_1:recv_field_1, send_field_2:send_field_2, send_field_3:recv_field_3, send_field_4:send_field_4",
            settings.get_field_list_string());
}

TEST_F(TransferMainInputFile, fieldListSingleEntryWithoutSequence)
{
  {
    stk::transfer_util::TransferInputFileOptions options =
      stk::transfer_util::parse_transfer_input_file_yaml("field-list: field_1");
    ASSERT_TRUE(options.fields.has_value());
    ASSERT_EQ(1u, options.fields->size());
    EXPECT_EQ("field_1", (*options.fields)[0].first);
    EXPECT_EQ("field_1", (*options.fields)[0].second);
  }
  {
    stk::transfer_util::TransferInputFileOptions options =
      stk::transfer_util::parse_transfer_input_file_yaml("field-list: {send: a, recv: b}");
    ASSERT_TRUE(options.fields.has_value());
    ASSERT_EQ(1u, options.fields->size());
    EXPECT_EQ("a", (*options.fields)[0].first);
    EXPECT_EQ("b", (*options.fields)[0].second);
  }
}

TEST_F(TransferMainInputFile, fieldListErrors)
{
  expect_yaml_error("field-list:\n  - send_field_1:recv_field_1\n", "contains ':'");
  expect_yaml_error("field-list: field_1, field_2\n", "contains ','");
  expect_yaml_error("field-list:\n  - {recv: recv_field_1}\n", "requires a 'send' field name");
  expect_yaml_error("field-list:\n  - {send: a, receive: b}\n", "unrecognized key 'receive'");
  expect_yaml_error("field-list:\n  - [a, b]\n", "must be a field name or a {send: name, recv: name} map");
  expect_yaml_error("field-list:\n", "must be a sequence");
  expect_yaml_error("field-list:\n  - {send: [a, b]}\n", "'field-list[0].send' must be a single value");
}

TEST_F(TransferMainInputFile, partListSendOnly)
{
  const std::string yaml = R"(
send-mesh: source.exo
recv-mesh: target.exo
part-list:
  send: [send_part_1, send_part_2]
)";
  stk::transfer_util::TransferMainSettings settings = parse_and_generate(yaml);

  ASSERT_EQ(2u, settings.get_transfer_send_parts().size());
  EXPECT_EQ("send_part_1", settings.get_transfer_send_parts()[0]);
  EXPECT_EQ("send_part_2", settings.get_transfer_send_parts()[1]);
  EXPECT_EQ(0u, settings.get_transfer_recv_parts().size());
}

TEST_F(TransferMainInputFile, partListRecvOnly)
{
  const std::string yaml = R"(
send-mesh: source.exo
recv-mesh: target.exo
part-list:
  recv:
    - recv_part_1
    - recv_part_2
)";
  stk::transfer_util::TransferMainSettings settings = parse_and_generate(yaml);

  EXPECT_EQ(0u, settings.get_transfer_send_parts().size());
  ASSERT_EQ(2u, settings.get_transfer_recv_parts().size());
  EXPECT_EQ("recv_part_1, recv_part_2", settings.get_recv_part_list_string());
}

TEST_F(TransferMainInputFile, partListSingleNamesAndEmptyLists)
{
  const std::string yaml = R"(
send-mesh: source.exo
recv-mesh: target.exo
part-list:
  send: block_1
  recv: []
)";
  stk::transfer_util::TransferMainSettings settings = parse_and_generate(yaml);

  ASSERT_EQ(1u, settings.get_transfer_send_parts().size());
  EXPECT_EQ("block_1", settings.get_transfer_send_parts()[0]);
  EXPECT_EQ(0u, settings.get_transfer_recv_parts().size());
}

TEST_F(TransferMainInputFile, partListErrors)
{
  expect_yaml_error("part-list: block_1:block_2\n", "'part-list' must be a map");
  expect_yaml_error("part-list: [block_1]\n", "'part-list' must be a map");
  expect_yaml_error("part-list:\n  sender: [block_1]\n", "unrecognized key 'sender'");
  expect_yaml_error("part-list:\n  send: {a: b}\n", "'part-list.send' must be a name or a sequence of names");
  expect_yaml_error("part-list:\n  send: [[a]]\n", "'part-list.send[0]' must be a single value");
}

TEST_F(TransferMainInputFile, timeStepsScalar)
{
  stk::transfer_util::TransferMainSettings settings =
    parse_and_generate("send-mesh: s.exo\nrecv-mesh: r.exo\ntime-steps: \"1:5:2\"\n");
  EXPECT_EQ("1:5:2", settings.get_time_steps_spec());
}

TEST_F(TransferMainInputFile, timeStepsMap)
{
  stk::transfer_util::TransferMainSettings settings =
    parse_and_generate("send-mesh: s.exo\nrecv-mesh: r.exo\ntime-steps:\n  begin: 1\n  end: 5\n  stride: 2\n");
  EXPECT_EQ("1:5:2", settings.get_time_steps_spec());
}

TEST_F(TransferMainInputFile, timeStepsMapDefaultStride)
{
  stk::transfer_util::TransferMainSettings settings =
    parse_and_generate("send-mesh: s.exo\nrecv-mesh: r.exo\ntime-steps: {begin: 2, end: 4}\n");
  EXPECT_EQ("2:4:1", settings.get_time_steps_spec());
}

TEST_F(TransferMainInputFile, timeStepsErrors)
{
  expect_yaml_error("time-steps: {begin: 1}\n", "requires both 'begin' and 'end'");
  expect_yaml_error("time-steps: {begin: 5, end: 1}\n", "begin <= end");
  expect_yaml_error("time-steps: {begin: 1, end: 5, stride: 0}\n", "positive stride");
  expect_yaml_error("time-steps: {begin: one, end: 5}\n", "'time-steps.begin' must be an integer");
  expect_yaml_error("time-steps: {first: 1, end: 5}\n", "unrecognized key 'first'");
  expect_yaml_error("time-steps: [1, 5]\n", "'time-steps' must be ALL, FIRST, LAST");
}

TEST_F(TransferMainInputFile, twoDTo3DAxisymmetricNotSpecified)
{
  stk::transfer_util::TransferMainSettings settings = parse_and_generate("send-mesh: s.exo\nrecv-mesh: r.exo\n");
  EXPECT_FALSE(settings.get_enable_2d_to_3d_axisymmetric_transfer());
}

TEST_F(TransferMainInputFile, twoDTo3DAxisymmetricDefaults)
{
  stk::transfer_util::TransferMainSettings settings =
    parse_and_generate("send-mesh: s.exo\nrecv-mesh: r.exo\n2d-to-3d-axisymmetric: {}\n");

  EXPECT_TRUE(settings.get_enable_2d_to_3d_axisymmetric_transfer());
  const stk::transfer_util::TwoDTo3DAxisymmetricParams& params = settings.get_2d_to_3d_axisymmetric_transfer_params();
  EXPECT_EQ(make_array(0, 0, 1), params.axis);
  EXPECT_EQ(make_array(0, 0, 0), params.trans);
}

TEST_F(TransferMainInputFile, twoDTo3DAxisymmetricAxisSpecified)
{
  stk::transfer_util::TransferMainSettings settings =
    parse_and_generate("send-mesh: s.exo\nrecv-mesh: r.exo\n2d-to-3d-axisymmetric:\n  axis: [1, 2, 3]\n");

  EXPECT_TRUE(settings.get_enable_2d_to_3d_axisymmetric_transfer());
  const stk::transfer_util::TwoDTo3DAxisymmetricParams& params = settings.get_2d_to_3d_axisymmetric_transfer_params();
  EXPECT_EQ(make_array(1, 2, 3), params.axis);
  EXPECT_EQ(make_array(0, 0, 0), params.trans);
}

TEST_F(TransferMainInputFile, twoDTo3DAxisymmetricAxisAndTranslationSpecified)
{
  const std::string yaml = R"(
send-mesh: s.exo
recv-mesh: r.exo
2d-to-3d-axisymmetric:
  axis: [1, 2, 3]
  translation: [4.5, 5, -6]
)";
  stk::transfer_util::TransferMainSettings settings = parse_and_generate(yaml);

  EXPECT_TRUE(settings.get_enable_2d_to_3d_axisymmetric_transfer());
  const stk::transfer_util::TwoDTo3DAxisymmetricParams& params = settings.get_2d_to_3d_axisymmetric_transfer_params();
  EXPECT_EQ(make_array(1, 2, 3), params.axis);
  EXPECT_EQ(make_array(4.5, 5, -6), params.trans);
}

TEST_F(TransferMainInputFile, twoDTo3DAxisymmetricErrors)
{
  expect_yaml_error("2d-to-3d-axisymmetric:\n", "'2d-to-3d-axisymmetric' must be a map");
  expect_yaml_error("2d-to-3d-axisymmetric: true\n", "'2d-to-3d-axisymmetric' must be a map");
  expect_yaml_error("2d-to-3d-axisymmetric: 0,1,0\n", "'2d-to-3d-axisymmetric' must be a map");
  expect_yaml_error("2d-to-3d-axisymmetric:\n  axis: [1, 2]\n", "'2d-to-3d-axisymmetric.axis' must be a sequence of 3 numbers");
  expect_yaml_error("2d-to-3d-axisymmetric:\n  translation: [1, 2, 3, 4]\n", "'2d-to-3d-axisymmetric.translation' must be a sequence of 3 numbers");
  expect_yaml_error("2d-to-3d-axisymmetric:\n  axis: [1, y, 3]\n", "'2d-to-3d-axisymmetric.axis[1]' must be a number");
  expect_yaml_error("2d-to-3d-axisymmetric:\n  trans: [1, 2, 3]\n", "unrecognized key 'trans'");
}

TEST_F(TransferMainInputFile, threeDTo3DAxisymmetricThetaOnly)
{
  stk::transfer_util::TransferMainSettings settings =
    parse_and_generate("send-mesh: s.exo\nrecv-mesh: r.exo\n3d-to-3d-axisymmetric: {theta-min: 30, theta-max: 60}\n");

  EXPECT_TRUE(settings.get_enable_3d_to_3d_axisymmetric_transfer());
  const stk::transfer_util::ThreeDTo3DAxisymmetricParams& params = settings.get_3d_to_3d_axisymmetric_transfer_params();
  EXPECT_NEAR(pi()/6, params.theta_min, 1e-13);
  EXPECT_NEAR(pi()/3, params.theta_max, 1e-13);
  EXPECT_EQ(make_array(0, 0, 1), params.axis);
  EXPECT_EQ(make_array(0, 0, 0), params.trans);
}

TEST_F(TransferMainInputFile, threeDTo3DAxisymmetricAllParams)
{
  const std::string yaml = R"(
send-mesh: s.exo
recv-mesh: r.exo
3d-to-3d-axisymmetric:
  theta-min: 30
  theta-max: 60
  axis: [1, 2, 3]
  translation: [4, 5, 6]
)";
  stk::transfer_util::TransferMainSettings settings = parse_and_generate(yaml);

  EXPECT_TRUE(settings.get_enable_3d_to_3d_axisymmetric_transfer());
  const stk::transfer_util::ThreeDTo3DAxisymmetricParams& params = settings.get_3d_to_3d_axisymmetric_transfer_params();
  EXPECT_NEAR(pi()/6, params.theta_min, 1e-13);
  EXPECT_NEAR(pi()/3, params.theta_max, 1e-13);
  EXPECT_EQ(make_array(1, 2, 3), params.axis);
  EXPECT_EQ(make_array(4, 5, 6), params.trans);
}

TEST_F(TransferMainInputFile, threeDTo3DAxisymmetricErrors)
{
  expect_yaml_error("3d-to-3d-axisymmetric: {theta-min: 30}\n", "requires both 'theta-min' and 'theta-max'");
  expect_yaml_error("3d-to-3d-axisymmetric: {}\n", "requires both 'theta-min' and 'theta-max'");
  expect_yaml_error("3d-to-3d-axisymmetric: 30,60\n", "'3d-to-3d-axisymmetric' must be a map");
  expect_yaml_error("3d-to-3d-axisymmetric: {theta-min: 30, theta-max: big}\n", "'3d-to-3d-axisymmetric.theta-max' must be a number");
  expect_yaml_error("3d-to-3d-axisymmetric: {theta-min: 30, theta-max: 60, axis: [0, 1]}\n", "must be a sequence of 3 numbers");
}

TEST_F(TransferMainInputFile, topLevelErrors)
{
  expect_yaml_error("send-mesh: s.exo\nfield-vector: [f]\n", "unrecognized key 'field-vector'");
  expect_yaml_error("s: s.exo\n", "unrecognized key 's'");
  expect_yaml_error("send-mesh: s.exo\nsend-mesh: t.exo\n", "more than once");
  expect_yaml_error("send-mesh: [s.exo, t.exo]\n", "'send-mesh' must be a single value");
  expect_yaml_error("send-mesh:\n", "'send-mesh' must be a single value");
  expect_yaml_error("send-mesh: ''\n", "'send-mesh' must not be empty");
  expect_yaml_error("send-mesh: s.exo\nrecv-mesh: [unclosed\n", "invalid YAML syntax");
  expect_yaml_error("- send-mesh: s.exo\n", "expected a map");
  expect_yaml_error("just a string\n", "expected a map");
  expect_yaml_error("", "expected a map");
}

TEST_F(TransferMainInputFile, invalidYamlThrowsFromParser)
{
  EXPECT_ANY_THROW(m_parser.set_input_file_contents("send-mesh: s.exo\nbogus-option: 1\n"));
}

TEST_F(TransferMainInputFile, missingRequiredMeshesIsParseError)
{
  EXPECT_EQ(TransferParserStatus::PARSE_ERROR, parse("output-mesh: out.exo\n"));
  EXPECT_NE(std::string::npos, m_stderr.find("send-mesh"));
}

TEST_F(TransferMainInputFile, missingRecvMeshIsParseError)
{
  EXPECT_EQ(TransferParserStatus::PARSE_ERROR, parse("send-mesh: s.exo\n"));
  EXPECT_NE(std::string::npos, m_stderr.find("recv-mesh"));
}

TEST_F(TransferMainInputFile, missingMeshMessageMentionsInputFile)
{
  stk::unit_test_util::Args args({"stk_transfer", "-r", "r.exo"});

  testing::internal::CaptureStderr();
  m_parser.parse_command_line_options(args.argc(), args.argv());
  const std::string errMsg = testing::internal::GetCapturedStderr();

  EXPECT_EQ(TransferParserStatus::PARSE_ERROR, m_parser.get_parser_status());
  EXPECT_NE(std::string::npos, errMsg.find("'--send-mesh' not found"));
  EXPECT_NE(std::string::npos, errMsg.find("--input-file"));
}

TEST_F(TransferMainInputFile, inputErrorsPrintOnlyTheErrorMessage)
{
  stk::unit_test_util::Args args({"stk_transfer", "-r", "r.exo"});

  testing::internal::CaptureStderr();
  m_parser.parse_command_line_options(args.argc(), args.argv());
  const std::string errMsg = testing::internal::GetCapturedStderr();

  EXPECT_EQ(TransferParserStatus::PARSE_ERROR, m_parser.get_parser_status());
  EXPECT_EQ(0u, errMsg.find("Error: required option '--send-mesh' not found.")) << errMsg;
  EXPECT_EQ(std::string::npos, errMsg.find("Requirement(")) << errMsg;
  EXPECT_EQ(std::string::npos, errMsg.find("Error occurred at")) << errMsg;
  EXPECT_EQ(1, std::count(errMsg.begin(), errMsg.end(), '\n')) << errMsg;
}

TEST_F(TransferMainInputFile, setInputFileContentsAfterParseThrows)
{
  stk::unit_test_util::Args args({"stk_transfer", "-s", "s.exo", "-r", "r.exo"});
  m_parser.parse_command_line_options(args.argc(), args.argv());
  EXPECT_EQ(TransferParserStatus::SUCCESS, m_parser.get_parser_status());

  EXPECT_ANY_THROW(m_parser.set_input_file_contents("output-mesh: out.exo\n"));
}

TEST_F(TransferMainInputFile, commandLineOverridesInputFile)
{
  const std::string yaml = R"(
send-mesh: file_source.exo
recv-mesh: file_target.exo
time-steps: FIRST
xfer-type: COPY
field-list: [file_field]
part-list: {send: [file_send_part], recv: [file_recv_part]}
2d-to-3d-axisymmetric: {axis: [1, 0, 0]}
3d-to-3d-axisymmetric: {theta-min: 10, theta-max: 20}
)";
  stk::transfer_util::TransferMainSettings settings =
    parse_and_generate(yaml, {"--send-mesh", "cli_source.exo", "--time-steps", "LAST", "--field-list", "cli_field",
                              "--part-list", "cli_send_part", "--2d-to-3d-axisymmetric", "0,1,0",
                              "--3d-to-3d-axisymmetric", "30,60"});

  EXPECT_EQ("cli_source.exo", settings.get_sendMesh_filename());
  EXPECT_EQ("file_target.exo", settings.get_recvMesh_filename());
  EXPECT_EQ("transferred_file_target.exo", settings.get_outputMesh_filename());
  EXPECT_EQ("LAST", settings.get_time_steps_spec());
  EXPECT_EQ("COPY", settings.get_transfer_type());
  EXPECT_EQ("cli_field:cli_field", settings.get_field_list_string());
  ASSERT_EQ(1u, settings.get_transfer_send_parts().size());
  EXPECT_EQ("cli_send_part", settings.get_transfer_send_parts()[0]);
  EXPECT_EQ(0u, settings.get_transfer_recv_parts().size());
  EXPECT_EQ(make_array(0, 1, 0), settings.get_2d_to_3d_axisymmetric_transfer_params().axis);
  EXPECT_NEAR(pi()/6, settings.get_3d_to_3d_axisymmetric_transfer_params().theta_min, 1e-13);
}

TEST_F(TransferMainInputFile, commandLineOutputMeshUsesRecvMeshFromFile)
{
  stk::transfer_util::TransferMainSettings settings =
    parse_and_generate("send-mesh: s.exo\nrecv-mesh: r.exo\noutput-mesh: file_out.exo\n", {"-o", "cli_out.exo"});
  EXPECT_EQ("cli_out.exo", settings.get_outputMesh_filename());
}

TEST_F(TransferMainInputFile, mixSendMeshFromFileRecvMeshFromCommandLine)
{
  stk::transfer_util::TransferMainSettings settings =
    parse_and_generate("send-mesh: source.exo\nrecv-type: element_centroid\n", {"-r", "target.exo"});

  EXPECT_EQ("source.exo", settings.get_sendMesh_filename());
  EXPECT_EQ("target.exo", settings.get_recvMesh_filename());
  EXPECT_EQ("transferred_target.exo", settings.get_outputMesh_filename());
  EXPECT_EQ("ELEMENT_CENTROID", settings.get_recv_type_string());
}

TEST_F(TransferMainInputFile, invalidValuesInFileThrowLikeCommandLine)
{
  {
    stk::transfer_util::TransferMainParser parser(m_comm);
    parser.set_input_file_contents("send-mesh: s.exo\nrecv-mesh: r.exo\nxfer-type: FAKE\n");
    stk::unit_test_util::Args args({"stk_transfer"});
    parser.parse_command_line_options(args.argc(), args.argv());
    EXPECT_EQ(TransferParserStatus::SUCCESS, parser.get_parser_status());
    EXPECT_ANY_THROW(parser.generate_transfer_settings());
  }
  {
    stk::transfer_util::TransferMainParser parser(m_comm);
    parser.set_input_file_contents("send-mesh: s.exo\nrecv-mesh: r.exo\nrecv-type: FAKE\n");
    stk::unit_test_util::Args args({"stk_transfer"});
    parser.parse_command_line_options(args.argc(), args.argv());
    EXPECT_ANY_THROW(parser.generate_transfer_settings());
  }
  {
    stk::transfer_util::TransferMainParser parser(m_comm);
    parser.set_input_file_contents("send-mesh: s.exo\nrecv-mesh: r.exo\nextrapolate-option: FAKE\n");
    stk::unit_test_util::Args args({"stk_transfer"});
    parser.parse_command_line_options(args.argc(), args.argv());
    EXPECT_ANY_THROW(parser.generate_transfer_settings());
  }
}

TEST_F(TransferMainInputFile, nonexistentInputFileIsParseError)
{
  stk::unit_test_util::Args args({"stk_transfer", "--input-file", "this_input_file_does_not_exist.yaml"});

  testing::internal::CaptureStderr();
  m_parser.parse_command_line_options(args.argc(), args.argv());
  const std::string errMsg = testing::internal::GetCapturedStderr();

  EXPECT_EQ(TransferParserStatus::PARSE_ERROR, m_parser.get_parser_status());
  EXPECT_NE(std::string::npos, errMsg.find("could not open input file"));
}

TEST_F(TransferMainInputFile, inputFileOptionIsRecognized)
{
  // -i is the short form of --input-file. Unrecognized options would be a parse error,
  // so a missing file (rather than an unrecognized option) proves the option exists.
  stk::unit_test_util::Args args({"stk_transfer", "-i", "this_input_file_does_not_exist.yaml", "-s", "s.exo", "-r", "r.exo"});

  testing::internal::CaptureStderr();
  m_parser.parse_command_line_options(args.argc(), args.argv());
  const std::string errMsg = testing::internal::GetCapturedStderr();

  EXPECT_EQ(TransferParserStatus::PARSE_ERROR, m_parser.get_parser_status());
  EXPECT_NE(std::string::npos, errMsg.find("could not open input file"));
}

#else // STK_HAS_YAML_CPP

TEST_F(TransferMainInputFile, notSupportedWithoutYamlCpp)
{
  EXPECT_FALSE(stk::transfer_util::transfer_input_file_supported());
  EXPECT_ANY_THROW(stk::transfer_util::parse_transfer_input_file_yaml("send-mesh: s.exo\n"));

  stk::unit_test_util::Args args({"stk_transfer", "--input-file", "input.yaml"});
  testing::internal::CaptureStderr();
  m_parser.parse_command_line_options(args.argc(), args.argv());
  const std::string errMsg = testing::internal::GetCapturedStderr();

  EXPECT_EQ(TransferParserStatus::PARSE_ERROR, m_parser.get_parser_status());
  EXPECT_NE(std::string::npos, errMsg.find("yaml-cpp"));
}

TEST_F(TransferMainInputFile, missingMeshMessageDoesNotMentionInputFileWithoutYamlCpp)
{
  stk::unit_test_util::Args args({"stk_transfer", "-r", "r.exo"});

  testing::internal::CaptureStderr();
  m_parser.parse_command_line_options(args.argc(), args.argv());
  const std::string errMsg = testing::internal::GetCapturedStderr();

  EXPECT_EQ(TransferParserStatus::PARSE_ERROR, m_parser.get_parser_status());
  EXPECT_NE(std::string::npos, errMsg.find("'--send-mesh' not found"));
  EXPECT_EQ(std::string::npos, errMsg.find("--input-file"));
}

TEST_F(TransferMainInputFile, commandLineStillWorksWithoutYamlCpp)
{
  stk::unit_test_util::Args args({"stk_transfer", "-s", "s.exo", "-r", "r.exo"});
  m_parser.parse_command_line_options(args.argc(), args.argv());
  EXPECT_EQ(TransferParserStatus::SUCCESS, m_parser.get_parser_status());
  EXPECT_EQ("s.exo", m_parser.generate_transfer_settings().get_sendMesh_filename());
}

#endif // STK_HAS_YAML_CPP

}
