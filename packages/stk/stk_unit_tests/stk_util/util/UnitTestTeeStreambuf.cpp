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
#include "stk_util/util/TeeStreambuf.hpp"  // for tee_streambuf
#include <sstream>                         // for ostringstream, basic_ostringstream<>::__string...
#include <string>                          // for operator==, allocator, operator<<, operator+


TEST(UnitTestTeeStreambuf, UnitTest)
{
  stk::tee_streambuf    out_tee_streambuf;

  std::ostream          my_out(&out_tee_streambuf);
  
  std::ostringstream    dest1;
  std::ostringstream    dest2;

  out_tee_streambuf.add(&dest1);
  out_tee_streambuf.add(&dest2);

  std::string message1("This is a test");
  std::string message2("This is a test");

  std::string message3 = message1 + message2;
  
  my_out << message1;

  ASSERT_EQ((dest1.str() == message1), true);
  ASSERT_EQ((dest2.str() == message1), true);

  out_tee_streambuf.remove(&dest2);

  my_out << message2;
  
  ASSERT_EQ((dest1.str() == message3), true);
  ASSERT_EQ((dest2.str() == message1), true);

  out_tee_streambuf.remove(&dest1);

  my_out << message2;

  ASSERT_EQ((dest1.str() == message3), true);
  ASSERT_EQ((dest2.str() == message1), true);
}

