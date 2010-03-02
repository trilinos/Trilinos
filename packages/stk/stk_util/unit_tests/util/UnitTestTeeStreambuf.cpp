/*------------------------------------------------------------------------*/
/*                 Copyright 2010 Sandia Corporation.                     */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/

#include <cppunit/TestCase.h>
#include <cppunit/extensions/HelperMacros.h>
#include <cppunit/extensions/TestFactoryRegistry.h>
#include <cppunit/ui/text/TestRunner.h>

#include <sstream>
#include <string>
#include <stk_util/util/TeeStreambuf.hpp>

class UnitTestTeeStreambuf : public CppUnit::TestCase {
private:
  CPPUNIT_TEST_SUITE(UnitTestTeeStreambuf);
  CPPUNIT_TEST(testUnit);
  CPPUNIT_TEST_SUITE_END();

public:
  void setUp()
  {}

  void tearDown()
  {}

  void testUnit();
};


void UnitTestTeeStreambuf::testUnit()
{
  stk::tee_streambuf    out_tee_streambuf;

  std::ostream          out(&out_tee_streambuf);
  
  std::ostringstream    dest1;
  std::ostringstream    dest2;

  out_tee_streambuf.add(&dest1);
  out_tee_streambuf.add(&dest2);

  std::string message1("This is a test");
  std::string message2("This is a test");

  std::string message3 = message1 + message2;
  
  out << message1;

  CPPUNIT_ASSERT_EQUAL(dest1.str(), message1);
  CPPUNIT_ASSERT_EQUAL(dest2.str(), message1);

  out_tee_streambuf.remove(&dest2);

  out << message2;
  
  CPPUNIT_ASSERT_EQUAL(dest1.str(), message3);
  CPPUNIT_ASSERT_EQUAL(dest2.str(), message1);

  out_tee_streambuf.remove(&dest1);

  out << message2;

  CPPUNIT_ASSERT_EQUAL(dest1.str(), message3);
  CPPUNIT_ASSERT_EQUAL(dest2.str(), message1);
}

CPPUNIT_TEST_SUITE_REGISTRATION(UnitTestTeeStreambuf);

