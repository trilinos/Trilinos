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

#include "gtest/gtest.h"
#include "stk_util/util/GetEnv.hpp"  // for get_env_var_as_bool, get_env_var_as_int
#include <cstdlib>                   // for unsetenv, setenv
#include <string>                    // for string

void check_env_var(bool expectedValue, const std::string& varName, const std::string& varVal)
{
  unsetenv(varName.c_str());
  setenv(varName.c_str(), varVal.c_str(), 1);

  bool boolValue = stk::get_env_var_as_bool(varName, !expectedValue);
  EXPECT_EQ(expectedValue, boolValue);

  unsetenv(varName.c_str());
}

void check_env_var(int expectedValue, const std::string& varName, const std::string& varVal,
                   int defaultValue)
{
  unsetenv(varName.c_str());
  setenv(varName.c_str(), varVal.c_str(), 1);

  int intValue = stk::get_env_var_as_int(varName, defaultValue);
  EXPECT_EQ(expectedValue, intValue);

  unsetenv(varName.c_str());
}

TEST(GetEnv, asBool)
{
  EXPECT_TRUE(stk::get_env_var_as_bool("VAR_DOESNT_EXIST", true));
  EXPECT_FALSE(stk::get_env_var_as_bool("VAR_DOESNT_EXIST", false));

  std::string varName("STK_TEST_BOOL");

  check_env_var(true, varName, "true");
  check_env_var(true, varName, "TRUE");
  check_env_var(true, varName, "on");
  check_env_var(true, varName, "On");

  check_env_var(false, varName, "false");
  check_env_var(false, varName, "FOO");
  check_env_var(false, varName, "ONLY");
  check_env_var(false, varName, "truly");
}

TEST(GetEnv, asInt)
{
  int defaultValue = -99;
  EXPECT_EQ(defaultValue, stk::get_env_var_as_int("VAR_DOESNT_EXIST", defaultValue));

  std::string varName("STK_TEST_INT");

  check_env_var(defaultValue, varName, "X100", defaultValue);
  check_env_var(defaultValue, varName, "foo", defaultValue);
  check_env_var(100, varName, "100", defaultValue);
}

