/*------------------------------------------------------------------------*/
/*                 Copyright 2010 Sandia Corporation.                     */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/

#include <sstream>
#include <string>
#include <stk_util/util/TeeStreambuf.hpp>

#include <stk_util/unit_test_support/stk_utest_macros.hpp>

STKUNIT_UNIT_TEST(UnitTestTeeStreambuf, UnitTest)
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

  STKUNIT_ASSERT_EQUAL(dest1.str(), message1);
  STKUNIT_ASSERT_EQUAL(dest2.str(), message1);

  out_tee_streambuf.remove(&dest2);

  my_out << message2;
  
  STKUNIT_ASSERT_EQUAL(dest1.str(), message3);
  STKUNIT_ASSERT_EQUAL(dest2.str(), message1);

  out_tee_streambuf.remove(&dest1);

  my_out << message2;

  STKUNIT_ASSERT_EQUAL(dest1.str(), message3);
  STKUNIT_ASSERT_EQUAL(dest2.str(), message1);
}

