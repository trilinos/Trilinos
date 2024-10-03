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
#include "stk_unit_test_utils/MeshFixture.hpp"
#include "stk_balance/setup/DefaultSettings.hpp"
#include <stk_balance/setup/M2NParser.hpp>
#include "stk_util/environment/Env.hpp"
#include <stk_unit_test_utils/TextMeshToFile.hpp>
#include <stk_balance/balanceUtils.hpp>
#include <stk_balance/m2n/m2nRebalance.hpp>
#include <vector>
#include <fstream>

namespace {

class TestM2NLogFile : public stk::unit_test_util::MeshFixture
{
protected:
  void clean_up_file(const std::string & fileName)
  {
    if (get_parallel_rank() == 0) {
      unlink(fileName.c_str());
    }
    MPI_Barrier(get_comm());
  }

  std::vector<const char*> assemble_args(const std::vector<const char*> & customArgs)
  {
    std::vector<const char*> args {"stk_balance_m2n", "dummy_mesh.g", "16"};

    for (const char * customArg : customArgs) {
      args.push_back(customArg);
    }

    return args;
  }

  void set_up_output_streams(const std::vector<const char*> & customArgs)
  {
    stk::balance::M2NBalanceSettings balanceSettings;
    std::vector<const char*> args = assemble_args(customArgs);
    stk::balance::M2NParser parser(get_comm());
    parser.parse_command_line_options(args.size(), args.data(), balanceSettings);
    stk::balance::m2n::set_output_streams(get_comm(), balanceSettings);
  }
};

bool test_file_exists(const std::string & fileName) {
  std::ifstream stream(fileName);
  return stream.good();
}

std::string get_file_contents(const std::string & fileName) {
  std::ifstream stream(fileName);
  std::stringstream buffer;
  buffer << stream.rdbuf();
  return buffer.str();
}

TEST_F(TestM2NLogFile, defaultLogFile)
{
  const int initialNumProcs = stk::parallel_machine_size(MPI_COMM_WORLD); 
  const std::string logFileName{"dummy_mesh." + std::to_string(initialNumProcs) + "_to_16.log"};
  clean_up_file(logFileName);
  set_up_output_streams({});

  const std::string expectedOutput{"This is a test\n"};
  sierra::Env::outputP0() << expectedOutput;
  sierra::Env::outputP0().flush();

  if (get_parallel_rank() == 0) {
    ASSERT_TRUE(test_file_exists(logFileName));
    EXPECT_EQ(get_file_contents(logFileName), expectedOutput);
  }

  clean_up_file(logFileName);
}

TEST_F(TestM2NLogFile, customLogFile)
{
  const std::string logFileName{"custom.log"};
  clean_up_file(logFileName);
  set_up_output_streams({"-l", logFileName.c_str()});

  const std::string expectedOutput{"This is a test\n"};
  sierra::Env::outputP0() << expectedOutput;
  sierra::Env::outputP0().flush();

  if (get_parallel_rank() == 0) {
    ASSERT_TRUE(test_file_exists(logFileName));
    EXPECT_EQ(get_file_contents(logFileName), expectedOutput);
  }

  clean_up_file(logFileName);
}

TEST_F(TestM2NLogFile, standardOutLog)
{
  const std::string logFileName{"cout"};
  set_up_output_streams({"-l", logFileName.c_str()});

  const std::string expectedOutput{"This is a test\n"};

  testing::internal::CaptureStdout();
  sierra::Env::outputP0() << expectedOutput;
  sierra::Env::outputP0().flush();
  std::string stdoutString = testing::internal::GetCapturedStdout();

  if (get_parallel_rank() == 0) {
    EXPECT_EQ(stdoutString, expectedOutput);
  }
}

TEST_F(TestM2NLogFile, standardErrLog)
{
  const std::string logFileName{"cerr"};
  set_up_output_streams({"-l", logFileName.c_str()});

  const std::string expectedOutput{"This is a test\n"};

  testing::internal::CaptureStderr();
  sierra::Env::outputP0() << expectedOutput;
  sierra::Env::outputP0().flush();
  std::string stderrString = testing::internal::GetCapturedStderr();

  if (get_parallel_rank() == 0) {
    EXPECT_EQ(stderrString, expectedOutput);
  }
}

}
