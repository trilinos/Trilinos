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

#include <iostream>

#include <stk_util/environment/CPUTime.hpp>

class UnitTestCPUTime : public CppUnit::TestCase {
private:
  CPPUNIT_TEST_SUITE( UnitTestCPUTime );
  CPPUNIT_TEST( testUnit );
  CPPUNIT_TEST_SUITE_END();

public:
  void setUp() {}
  void tearDown() {}
  void testUnit();
};

CPPUNIT_TEST_SUITE_REGISTRATION(UnitTestCPUTime);

void
UnitTestCPUTime::testUnit()
{
  std::ostringstream oss;
  
  double cpu_now = stk::cpu_time();
  
  double x = 0.0;
  for (int i = 0; i < 10000; ++i)
    x += 1.0;

  // This makes sure that the loop isn't optimized away (hopefully)
  if (x > 100000.0)
    oss << x << std::endl;
  
  double cpu_delta = stk::cpu_time() - cpu_now;
  
  CPPUNIT_ASSERT(cpu_delta >= 0.0 && cpu_delta <= 1.0);
}
