/*------------------------------------------------------------------------*/
/*                 Copyright 2010 Sandia Corporation.                     */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/

#include <iostream>
#include <stdexcept>

#include <stk_util/environment/FormatTime.hpp>

#include <stk_util/unit_test_support/stk_utest_macros.hpp>

STKUNIT_UNIT_TEST(UnitTestFormatTime, UnitTest)
{
  double time = 41.399684906;
  
  STKUNIT_ASSERT_EQUAL((std::string("41") == stk::formatTime(time, stk::TIMEFORMAT_HMS)), true);
  STKUNIT_ASSERT_EQUAL((std::string("41.400") == stk::formatTime(time, stk::TIMEFORMAT_HMS | stk::TIMEFORMAT_MILLIS)), true);
  STKUNIT_ASSERT_EQUAL((std::string("41") == stk::formatTime(time, stk::TIMEFORMAT_SECONDS)), true);
  STKUNIT_ASSERT_EQUAL((std::string("41.400") == stk::formatTime(time, stk::TIMEFORMAT_SECONDS | stk::TIMEFORMAT_MILLIS)), true);
  STKUNIT_ASSERT_EQUAL((std::string("41.3997") == stk::formatTime(time, stk::TIMEFORMAT_NONE)), true);

  time = 441.399684906;
  
  STKUNIT_ASSERT_EQUAL((std::string("7:21") == stk::formatTime(time, stk::TIMEFORMAT_HMS)), true);
  STKUNIT_ASSERT_EQUAL((std::string("7:21.400") == stk::formatTime(time, stk::TIMEFORMAT_HMS | stk::TIMEFORMAT_MILLIS)), true);
  STKUNIT_ASSERT_EQUAL((std::string("441") == stk::formatTime(time, stk::TIMEFORMAT_SECONDS)), true);
  STKUNIT_ASSERT_EQUAL((std::string("441.400") == stk::formatTime(time, stk::TIMEFORMAT_SECONDS | stk::TIMEFORMAT_MILLIS)), true);
  STKUNIT_ASSERT_EQUAL((std::string("441.4") == stk::formatTime(time, stk::TIMEFORMAT_NONE)), true);

  time = 5441.399684906;
  
  STKUNIT_ASSERT_EQUAL((std::string("1:30:41") == stk::formatTime(time, stk::TIMEFORMAT_HMS)), true);
  STKUNIT_ASSERT_EQUAL((std::string("1:30:41.400") == stk::formatTime(time, stk::TIMEFORMAT_HMS | stk::TIMEFORMAT_MILLIS)), true);
  STKUNIT_ASSERT_EQUAL((std::string("5441") == stk::formatTime(time, stk::TIMEFORMAT_SECONDS)), true);
  STKUNIT_ASSERT_EQUAL((std::string("5441.400") == stk::formatTime(time, stk::TIMEFORMAT_SECONDS | stk::TIMEFORMAT_MILLIS)), true);
  STKUNIT_ASSERT_EQUAL((std::string("5441.4") == stk::formatTime(time, stk::TIMEFORMAT_NONE)), true);

  time = 1251305441.399684906;
  
  STKUNIT_ASSERT_EQUAL((std::string("347584:50:41") == stk::formatTime(time, stk::TIMEFORMAT_HMS)), true);
  STKUNIT_ASSERT_EQUAL((std::string("347584:50:41.400") == stk::formatTime(time, stk::TIMEFORMAT_HMS | stk::TIMEFORMAT_MILLIS)), true);
  STKUNIT_ASSERT_EQUAL((std::string("1251305441") == stk::formatTime(time, stk::TIMEFORMAT_SECONDS)), true);
  STKUNIT_ASSERT_EQUAL((std::string("1251305441.400") == stk::formatTime(time, stk::TIMEFORMAT_SECONDS | stk::TIMEFORMAT_MILLIS)), true);
  STKUNIT_ASSERT_EQUAL((std::string("1.25131e+09") == stk::formatTime(time, stk::TIMEFORMAT_NONE)), true);

  time = -41.399684906;
  
  STKUNIT_ASSERT_EQUAL((std::string("-41") == stk::formatTime(time, stk::TIMEFORMAT_HMS)), true);
  STKUNIT_ASSERT_EQUAL((std::string("-41.400") == stk::formatTime(time, stk::TIMEFORMAT_HMS | stk::TIMEFORMAT_MILLIS)), true);
  STKUNIT_ASSERT_EQUAL((std::string("-41") == stk::formatTime(time, stk::TIMEFORMAT_SECONDS)), true);
  STKUNIT_ASSERT_EQUAL((std::string("-41.400") == stk::formatTime(time, stk::TIMEFORMAT_SECONDS | stk::TIMEFORMAT_MILLIS)), true);
  STKUNIT_ASSERT_EQUAL((std::string("-41.3997") == stk::formatTime(time, stk::TIMEFORMAT_NONE)), true);
}



