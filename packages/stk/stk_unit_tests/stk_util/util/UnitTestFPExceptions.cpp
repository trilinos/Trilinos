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
#include "stk_util/stk_config.h"
#include "stk_util/util/FPExceptions.hpp"
#include <sstream>

TEST(FPExceptions, SimpleAdditionNoError)
{
  if (!stk::util::have_errno() && !stk::util::have_errexcept()) GTEST_SKIP();

  stk::util::clear_fp_errors();
  double x = 1.0 + 2.0;
  EXPECT_NO_THROW(stk::util::throw_on_fp_error());
  EXPECT_EQ(x, 3.0);  // appease the compiler
}

TEST(FPExceptions, Log0Error)
{
  if (!stk::util::have_errno() && !stk::util::have_errexcept()) GTEST_SKIP();

  stk::util::clear_fp_errors();
  std::log(0.0);
  EXPECT_ANY_THROW(stk::util::throw_on_fp_error());
}

TEST(FPExceptions, FlagsAreClearedAfterThrow)
{
  if (!stk::util::have_errno() && !stk::util::have_errexcept()) GTEST_SKIP();

  stk::util::clear_fp_errors();
  std::log(0.0);
  EXPECT_ANY_THROW(stk::util::throw_on_fp_error());
  EXPECT_NO_THROW(stk::util::throw_on_fp_error());
}

TEST(FPExceptions, ErrorMessageContainsName)
{
  if (!stk::util::have_errno() && !stk::util::have_errexcept()) GTEST_SKIP();

  stk::util::clear_fp_errors();
  std::log(0.0);

  std::string fname = "my_very_specific_and_clear_function_name";
  try {
    stk::util::throw_on_fp_error(fname.c_str());
  } catch (std::exception& except)
  {
    std::string msg(except.what());
    size_t pos = msg.find(fname);
    EXPECT_NE(pos, std::string::npos);
  }
}

TEST(FPExceptions, NoWarning)
{
  if (!stk::util::have_errno() && !stk::util::have_errexcept()) GTEST_SKIP();

  stk::util::clear_fp_errors();
  double x = 1.0 + 2.0;
  std::stringstream ss;
  EXPECT_FALSE(stk::util::warn_on_fp_error(nullptr, ss));
  EXPECT_EQ(ss.str().size(), 0U);
  EXPECT_EQ(x, 3.0);  // appease the compiler
}

TEST(FPExceptions, Warning)
{
  if (!stk::util::have_errno() && !stk::util::have_errexcept()) GTEST_SKIP();

  stk::util::clear_fp_errors();
  std::log(0.0);
  std::stringstream ss;
  EXPECT_TRUE(stk::util::warn_on_fp_error(nullptr, ss));
  EXPECT_GT(ss.str().size(), 0U);
}
