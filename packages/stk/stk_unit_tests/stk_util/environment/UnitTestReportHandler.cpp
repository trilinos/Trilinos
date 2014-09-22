// Copyright (c) 2013, Sandia Corporation.
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Governement retains certain rights in this software.
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
//     * Neither the name of Sandia Corporation nor the names of its
//       contributors may be used to endorse or promote products derived
//       from this software without specific prior written permission.
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

#include <iostream>                     // for operator<<, ostringstream, etc
#include <stdexcept>                    // for runtime_error
#include <stk_util/environment/ReportHandler.hpp>
#include <gtest/gtest.h>
#include <string>                       // for operator==, basic_string, etc



namespace {

std::ostringstream &
test_ostringstream() 
{
  static std::ostringstream     s_testStringStream;

  return s_testStringStream;
}


std::ostream &
test_stream() 
{
  return test_ostringstream();
}


void
test_report_handler(
  const char *          message,
  int                   type)
{
  test_stream() << "Message type " << type << ": " << message << std::endl;
}

} // namespace <empty>

TEST(UnitTestReportHandler, UnitTest)
{
  // Set and restore report handler.
  stk::report("This is a test", 0);

  stk::REH original_reh = stk::set_report_handler(test_report_handler);

  stk::report("This is a test", 0);
  
  stk::set_report_handler(original_reh);

  ASSERT_THROW(stk::set_report_handler(0), std::runtime_error);
  
  ASSERT_EQ((std::string("Message type 0: This is a test\n") == test_ostringstream().str()), true);
  ASSERT_EQ((std::string("Test.cpp") == stk::source_relative_path("/src/Test.cpp")), true);
  ASSERT_EQ((std::string("Test.hpp") == stk::source_relative_path("/include/Test.hpp")), true);
  ASSERT_EQ((std::string("Apps_Test.cpp") == stk::source_relative_path("/Apps_Test.cpp")), true);
  ASSERT_EQ((std::string("stk_Test.cpp") == stk::source_relative_path("/stk_Test.cpp")), true);
  ASSERT_EQ((std::string("/smile/Test.cpp") == stk::source_relative_path("/smile/Test.cpp")), true);
}

