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

#include <gtest/gtest.h>                // for AssertHelper, ASSERT_TRUE, etc
#include <stk_util/stk_config.h>
#include <stk_util/util/ParameterList.hpp>
#include <cstdint>
#include <string>
#include <vector>

namespace {

TEST(UnitTestParameterList, set_and_get_int)
{
  stk::util::ParameterList params;
  EXPECT_EQ(0, std::distance(params.begin(), params.end()));

  int intValue = 3;

  params.set_param("intValue", intValue, true, true);

  EXPECT_EQ(1u, std::distance(params.begin(), params.end()));

  stk::util::Parameter& param = params.get_param("intValue");
  EXPECT_TRUE(stk::util::ParameterType::INTEGER == param.type);
  EXPECT_TRUE(param.toResultsFile);
  EXPECT_TRUE(param.toRestartFile);

  EXPECT_EQ(intValue, params.get_value<int>("intValue"));
}

TEST(UnitTestParameterList, set_and_get_vectorint)
{
  stk::util::ParameterList params;
  EXPECT_EQ(0, std::distance(params.begin(), params.end()));

  std::vector<int> vectorIntValues = {1, 2, 3};

  params.set_value("vectorIntValue", vectorIntValues);

  EXPECT_EQ(1u, std::distance(params.begin(), params.end()));

  stk::util::Parameter& param = params.get_param("vectorIntValue");
  EXPECT_TRUE(stk::util::ParameterType::INTEGERVECTOR == param.type);
  EXPECT_FALSE(param.toResultsFile);
  EXPECT_FALSE(param.toRestartFile);

  EXPECT_EQ(vectorIntValues, params.get_value<std::vector<int>>("vectorIntValue"));
}

TEST(UnitTestParameterList, set_and_find)
{
  stk::util::ParameterList params;

  int intValue = 3;
  std::vector<int> vectorIntValues = {1, 2, 3};

  params.set_param("intValue", intValue, true, true);
  params.set_value("vectorIntValue", vectorIntValues);

  EXPECT_TRUE(params.find("intValue") != params.end());
  EXPECT_TRUE(params.find("vectorIntValue") != params.end());
  EXPECT_TRUE(params.find("badName") == params.end());
}

}

