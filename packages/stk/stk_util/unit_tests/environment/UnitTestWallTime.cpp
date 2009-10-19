#include <cppunit/TestCase.h>
#include <cppunit/extensions/HelperMacros.h>
#include <cppunit/extensions/TestFactoryRegistry.h>
#include <cppunit/ui/text/TestRunner.h>

#include <unistd.h>

#include <iostream>

#include <stk_util/environment/WallTime.hpp>

class UnitTestWallTime : public CppUnit::TestCase {
private:
  CPPUNIT_TEST_SUITE( UnitTestWallTime );
  CPPUNIT_TEST( testUnit );
  CPPUNIT_TEST_SUITE_END();

public:
  void setUp() {}
  void tearDown() {}
  void testUnit();
};

CPPUNIT_TEST_SUITE_REGISTRATION(UnitTestWallTime);

void UnitTestWallTime::testUnit()
{
  double wall_now = stk::wall_time();
  
  ::sleep(1);
  
  double wall_delta = stk::wall_time() - wall_now;
  
  CPPUNIT_ASSERT(wall_delta >= 1.0 && wall_delta <= 2.0);

  double wall_delta2 = stk::wall_dtime(wall_now);

  CPPUNIT_ASSERT(wall_delta2 >= 1.0 && wall_delta2 <= 2.0);
}
