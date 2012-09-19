/*------------------------------------------------------------------------*/
/*                 Copyright 2010 Sandia Corporation.                     */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/

#include <iostream>
#include <stdexcept>

#include <stk_util/environment/ReportHandler.hpp>

#include <stk_util/unit_test_support/stk_utest_macros.hpp>

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

STKUNIT_UNIT_TEST(UnitTestReportHandler, UnitTest)
{
  // Set and restore report handler.
  stk::report("This is a test", 0);

  stk::REH original_reh = stk::set_report_handler(test_report_handler);

  stk::report("This is a test", 0);
  
  stk::set_report_handler(original_reh);

  STKUNIT_ASSERT_THROW(stk::set_report_handler(0), std::runtime_error);
  
  STKUNIT_ASSERT_EQUAL((std::string("Message type 0: This is a test\n") == test_ostringstream().str()), true);
  STKUNIT_ASSERT_EQUAL((std::string("Test.cpp") == stk::source_relative_path("/src/Test.cpp")), true);
  STKUNIT_ASSERT_EQUAL((std::string("Test.hpp") == stk::source_relative_path("/include/Test.hpp")), true);
  STKUNIT_ASSERT_EQUAL((std::string("Apps_Test.cpp") == stk::source_relative_path("/Apps_Test.cpp")), true);
  STKUNIT_ASSERT_EQUAL((std::string("stk_Test.cpp") == stk::source_relative_path("/stk_Test.cpp")), true);
  STKUNIT_ASSERT_EQUAL((std::string("/smile/Test.cpp") == stk::source_relative_path("/smile/Test.cpp")), true);
}

