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
#include "stk_util/environment/FormatTime.hpp"  // for formatTime, TIMEFORMAT_HMS, TIMEFORMAT_MI...
#include <string>                               // for operator==, string



TEST(UnitTestFormatTime, UnitTest)
{
  double time = 41.399684906;
  
  ASSERT_EQ(std::string("00:00:41"), stk::formatTime(time, stk::TIMEFORMAT_HMS));
  ASSERT_EQ(std::string("00:00:41.400"), stk::formatTime(time, stk::TIMEFORMAT_HMS | stk::TIMEFORMAT_MILLIS));
  ASSERT_EQ(std::string("41"), stk::formatTime(time, stk::TIMEFORMAT_SECONDS));
  ASSERT_EQ(std::string("41.400"), stk::formatTime(time, stk::TIMEFORMAT_SECONDS | stk::TIMEFORMAT_MILLIS));
  ASSERT_EQ(std::string("41.3997"), stk::formatTime(time, stk::TIMEFORMAT_NONE));

  time = 441.399684906;
  
  ASSERT_EQ(std::string("00:07:21"), stk::formatTime(time, stk::TIMEFORMAT_HMS));
  ASSERT_EQ(std::string("00:07:21.400"), stk::formatTime(time, stk::TIMEFORMAT_HMS | stk::TIMEFORMAT_MILLIS));
  ASSERT_EQ(std::string("441"), stk::formatTime(time, stk::TIMEFORMAT_SECONDS));
  ASSERT_EQ(std::string("441.400"), stk::formatTime(time, stk::TIMEFORMAT_SECONDS | stk::TIMEFORMAT_MILLIS));
  ASSERT_EQ(std::string("441.4"), stk::formatTime(time, stk::TIMEFORMAT_NONE));

  time = 5441.399684906;
  
  ASSERT_EQ(std::string("01:30:41"), stk::formatTime(time, stk::TIMEFORMAT_HMS));
  ASSERT_EQ(std::string("01:30:41.400"), stk::formatTime(time, stk::TIMEFORMAT_HMS | stk::TIMEFORMAT_MILLIS));
  ASSERT_EQ(std::string("5441"), stk::formatTime(time, stk::TIMEFORMAT_SECONDS));
  ASSERT_EQ(std::string("5441.400"), stk::formatTime(time, stk::TIMEFORMAT_SECONDS | stk::TIMEFORMAT_MILLIS));
  ASSERT_EQ(std::string("5441.4"), stk::formatTime(time, stk::TIMEFORMAT_NONE));

  time = 1251305441.399684906;
  
  ASSERT_EQ(std::string("347584:50:41"), stk::formatTime(time, stk::TIMEFORMAT_HMS));
  ASSERT_EQ(std::string("347584:50:41.400"), stk::formatTime(time, stk::TIMEFORMAT_HMS | stk::TIMEFORMAT_MILLIS));
  ASSERT_EQ(std::string("1251305441"), stk::formatTime(time, stk::TIMEFORMAT_SECONDS));
  ASSERT_EQ(std::string("1251305441.400"), stk::formatTime(time, stk::TIMEFORMAT_SECONDS | stk::TIMEFORMAT_MILLIS));
  ASSERT_EQ(std::string("1.25131e+09"), stk::formatTime(time, stk::TIMEFORMAT_NONE));

  time = -41.399684906;
  
  ASSERT_EQ(std::string("-00:00:41"), stk::formatTime(time, stk::TIMEFORMAT_HMS));
  ASSERT_EQ(std::string("-00:00:41.400"), stk::formatTime(time, stk::TIMEFORMAT_HMS | stk::TIMEFORMAT_MILLIS));
  ASSERT_EQ(std::string("-41"), stk::formatTime(time, stk::TIMEFORMAT_SECONDS));
  ASSERT_EQ(std::string("-41.400"), stk::formatTime(time, stk::TIMEFORMAT_SECONDS | stk::TIMEFORMAT_MILLIS));
  ASSERT_EQ(std::string("-41.3997"), stk::formatTime(time, stk::TIMEFORMAT_NONE));
}



