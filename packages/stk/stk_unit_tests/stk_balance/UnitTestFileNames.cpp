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
#include "stk_balance/setup/FileValidator.hpp"
#include "stk_balance/setup/Parser.hpp"

class InputSanity : public stk::unit_test_util::MeshFixture
{
public:
  InputSanity()
    : validator(get_comm())
  { }

protected:
  stk::balance::FileValidator validator;
};

TEST_F(InputSanity, verifyThrowIfInputFileEqualsOutputFile)
{
  if (get_parallel_size() == 1) {
    std::string serialMeshName = "sixteen_hex_transient.e";
    std::string parallelOutputMeshName = "sixteen_hex_transient.e";

    EXPECT_TRUE(validator.serial_input_equals_output(serialMeshName, parallelOutputMeshName));
  }
}

TEST_F(InputSanity, verifyThrowIfInputFileEqualsDotSlashOutputFile)
{
  if (get_parallel_size() == 1) {
    std::string serialMeshName = "sixteen_hex_transient.e";
    std::string parallelOutputMeshName = "./sixteen_hex_transient.e";

    EXPECT_TRUE(validator.serial_input_equals_output(serialMeshName, parallelOutputMeshName));
  }
}

TEST(BalanceOutputFile, checkOutputFileFromInputFileWithoutPath)
{
  EXPECT_EQ(stk::balance::construct_output_file_name("/home/code/results", "input.e"), "/home/code/results/input.e");
}

TEST(BalanceOutputFile, checkOutputFileFromInputFileAndOutputFileWithExtraBackslash)
{
  EXPECT_EQ(stk::balance::construct_output_file_name("/home/code/results/", "input.e"), "/home/code/results//input.e");
}

TEST(BalanceOutputFile, checkOutputFileFromInputFileWithPath)
{
  EXPECT_EQ(stk::balance::construct_output_file_name("/home/code/results", "/another/directory/input.e"), "/home/code/results/input.e");
}

