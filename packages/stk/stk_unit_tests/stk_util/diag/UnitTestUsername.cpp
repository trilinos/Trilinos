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
#include "stk_util/diag/Platform.hpp"  // for domainname, hostname
#include <stdio.h>                     // for sprintf
#include <stdlib.h>                     // for putenv
#include <cstdlib>                     // for getenv
#include <string>                      // for string, operator==, basic_string


TEST(StkEnv, getenv)
{
  std::string userKey("MYUSER=");
  putenv(userKey.data());

  const char* env_user = std::getenv("MYUSER");

  //getenv can return non-null but empty strings...
  if (env_user) {
    std::string str_env_user(env_user);
    if (str_env_user.empty()) {
      env_user = NULL;
    }
  }

  char*expected_env_user = NULL;

  EXPECT_EQ(expected_env_user, env_user);

  unsetenv("MYUSER");
  std::string userKeyVal("MYUSER=the_man");
  putenv(userKeyVal.data());

  env_user = std::getenv("MYUSER");
  std::string str_env_user(env_user);

  std::string expected_user("the_man");

  bool they_match = expected_user == str_env_user;
  EXPECT_TRUE(they_match);

  unsetenv("MYUSER");
}

TEST(StkEnv, HostName)
{
  std::string host_name = sierra::Env::hostname();
  EXPECT_TRUE( !host_name.empty() );
}

TEST(StkEnv, DomainName)
{
  std::string host_name = sierra::Env::domainname();
  EXPECT_TRUE( !host_name.empty() );
}

