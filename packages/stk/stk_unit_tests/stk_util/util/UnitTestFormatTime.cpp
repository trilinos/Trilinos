/*------------------------------------------------------------------------*/
/*                 Copyright (c) 2013, Sandia Corporation.
/*                 Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
/*                 the U.S. Governement retains certain rights in this software.
/*                 
/*                 Redistribution and use in source and binary forms, with or without
/*                 modification, are permitted provided that the following conditions are
/*                 met:
/*                 
/*                     * Redistributions of source code must retain the above copyright
/*                       notice, this list of conditions and the following disclaimer.
/*                 
/*                     * Redistributions in binary form must reproduce the above
/*                       copyright notice, this list of conditions and the following
/*                       disclaimer in the documentation and/or other materials provided
/*                       with the distribution.
/*                 
/*                     * Neither the name of Sandia Corporation nor the names of its
/*                       contributors may be used to endorse or promote products derived
/*                       from this software without specific prior written permission.
/*                 
/*                 THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
/*                 "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
/*                 LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
/*                 A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
/*                 OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
/*                 SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
/*                 LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
/*                 DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
/*                 THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
/*                 (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
/*                 OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
/*                 
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/

#include <stk_util/environment/FormatTime.hpp>  // for formatTime, etc
#include <gtest/gtest.h>
#include <string>                       // for operator==, basic_string, etc



TEST(UnitTestFormatTime, UnitTest)
{
  double time = 41.399684906;
  
  ASSERT_EQ((std::string("41") == stk::formatTime(time, stk::TIMEFORMAT_HMS)), true);
  ASSERT_EQ((std::string("41.400") == stk::formatTime(time, stk::TIMEFORMAT_HMS | stk::TIMEFORMAT_MILLIS)), true);
  ASSERT_EQ((std::string("41") == stk::formatTime(time, stk::TIMEFORMAT_SECONDS)), true);
  ASSERT_EQ((std::string("41.400") == stk::formatTime(time, stk::TIMEFORMAT_SECONDS | stk::TIMEFORMAT_MILLIS)), true);
  ASSERT_EQ((std::string("41.3997") == stk::formatTime(time, stk::TIMEFORMAT_NONE)), true);

  time = 441.399684906;
  
  ASSERT_EQ((std::string("7:21") == stk::formatTime(time, stk::TIMEFORMAT_HMS)), true);
  ASSERT_EQ((std::string("7:21.400") == stk::formatTime(time, stk::TIMEFORMAT_HMS | stk::TIMEFORMAT_MILLIS)), true);
  ASSERT_EQ((std::string("441") == stk::formatTime(time, stk::TIMEFORMAT_SECONDS)), true);
  ASSERT_EQ((std::string("441.400") == stk::formatTime(time, stk::TIMEFORMAT_SECONDS | stk::TIMEFORMAT_MILLIS)), true);
  ASSERT_EQ((std::string("441.4") == stk::formatTime(time, stk::TIMEFORMAT_NONE)), true);

  time = 5441.399684906;
  
  ASSERT_EQ((std::string("1:30:41") == stk::formatTime(time, stk::TIMEFORMAT_HMS)), true);
  ASSERT_EQ((std::string("1:30:41.400") == stk::formatTime(time, stk::TIMEFORMAT_HMS | stk::TIMEFORMAT_MILLIS)), true);
  ASSERT_EQ((std::string("5441") == stk::formatTime(time, stk::TIMEFORMAT_SECONDS)), true);
  ASSERT_EQ((std::string("5441.400") == stk::formatTime(time, stk::TIMEFORMAT_SECONDS | stk::TIMEFORMAT_MILLIS)), true);
  ASSERT_EQ((std::string("5441.4") == stk::formatTime(time, stk::TIMEFORMAT_NONE)), true);

  time = 1251305441.399684906;
  
  ASSERT_EQ((std::string("347584:50:41") == stk::formatTime(time, stk::TIMEFORMAT_HMS)), true);
  ASSERT_EQ((std::string("347584:50:41.400") == stk::formatTime(time, stk::TIMEFORMAT_HMS | stk::TIMEFORMAT_MILLIS)), true);
  ASSERT_EQ((std::string("1251305441") == stk::formatTime(time, stk::TIMEFORMAT_SECONDS)), true);
  ASSERT_EQ((std::string("1251305441.400") == stk::formatTime(time, stk::TIMEFORMAT_SECONDS | stk::TIMEFORMAT_MILLIS)), true);
  ASSERT_EQ((std::string("1.25131e+09") == stk::formatTime(time, stk::TIMEFORMAT_NONE)), true);

  time = -41.399684906;
  
  ASSERT_EQ((std::string("-41") == stk::formatTime(time, stk::TIMEFORMAT_HMS)), true);
  ASSERT_EQ((std::string("-41.400") == stk::formatTime(time, stk::TIMEFORMAT_HMS | stk::TIMEFORMAT_MILLIS)), true);
  ASSERT_EQ((std::string("-41") == stk::formatTime(time, stk::TIMEFORMAT_SECONDS)), true);
  ASSERT_EQ((std::string("-41.400") == stk::formatTime(time, stk::TIMEFORMAT_SECONDS | stk::TIMEFORMAT_MILLIS)), true);
  ASSERT_EQ((std::string("-41.3997") == stk::formatTime(time, stk::TIMEFORMAT_NONE)), true);
}



