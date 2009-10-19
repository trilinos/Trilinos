#include <cppunit/TestCase.h>
#include <cppunit/extensions/HelperMacros.h>
#include <cppunit/extensions/TestFactoryRegistry.h>
#include <cppunit/ui/text/TestRunner.h>

#include <iostream>
#include <stdexcept>

#include <mpi.h>

// #include <stk_util/use_cases/UseCaseEnvironment.hpp>
#include <stk_util/util/FormatTime.hpp>

// using namespace use_case;


class UnitTestFormatTime : public CppUnit::TestCase {
private:
  CPPUNIT_TEST_SUITE(UnitTestFormatTime);
  CPPUNIT_TEST(testUnit);
  CPPUNIT_TEST_SUITE_END();

public:
  void setUp()
  {}
  
  void tearDown()
  {}

  void testUnit();
};

CPPUNIT_TEST_SUITE_REGISTRATION(UnitTestFormatTime);

void
UnitTestFormatTime::testUnit()
{
  double time = 41.399684906;
  
  CPPUNIT_ASSERT_EQUAL(std::string("41"), stk::formatTime(time, stk::TIMEFORMAT_HMS));
  CPPUNIT_ASSERT_EQUAL(std::string("41.400"), stk::formatTime(time, stk::TIMEFORMAT_HMS | stk::TIMEFORMAT_MILLIS));
  CPPUNIT_ASSERT_EQUAL(std::string("41"), stk::formatTime(time, stk::TIMEFORMAT_SECONDS));
  CPPUNIT_ASSERT_EQUAL(std::string("41.400"), stk::formatTime(time, stk::TIMEFORMAT_SECONDS | stk::TIMEFORMAT_MILLIS));
  CPPUNIT_ASSERT_EQUAL(std::string("41.3997"), stk::formatTime(time, stk::TIMEFORMAT_NONE));

  time = 441.399684906;
  
  CPPUNIT_ASSERT_EQUAL(std::string("7:21"), stk::formatTime(time, stk::TIMEFORMAT_HMS));
  CPPUNIT_ASSERT_EQUAL(std::string("7:21.400"), stk::formatTime(time, stk::TIMEFORMAT_HMS | stk::TIMEFORMAT_MILLIS));
  CPPUNIT_ASSERT_EQUAL(std::string("441"), stk::formatTime(time, stk::TIMEFORMAT_SECONDS));
  CPPUNIT_ASSERT_EQUAL(std::string("441.400"), stk::formatTime(time, stk::TIMEFORMAT_SECONDS | stk::TIMEFORMAT_MILLIS));
  CPPUNIT_ASSERT_EQUAL(std::string("441.4"), stk::formatTime(time, stk::TIMEFORMAT_NONE));

  time = 5441.399684906;
  
  CPPUNIT_ASSERT_EQUAL(std::string("1:30:41"), stk::formatTime(time, stk::TIMEFORMAT_HMS));
  CPPUNIT_ASSERT_EQUAL(std::string("1:30:41.400"), stk::formatTime(time, stk::TIMEFORMAT_HMS | stk::TIMEFORMAT_MILLIS));
  CPPUNIT_ASSERT_EQUAL(std::string("5441"), stk::formatTime(time, stk::TIMEFORMAT_SECONDS));
  CPPUNIT_ASSERT_EQUAL(std::string("5441.400"), stk::formatTime(time, stk::TIMEFORMAT_SECONDS | stk::TIMEFORMAT_MILLIS));
  CPPUNIT_ASSERT_EQUAL(std::string("5441.4"), stk::formatTime(time, stk::TIMEFORMAT_NONE));

  time = 1251305441.399684906;
  
  CPPUNIT_ASSERT_EQUAL(std::string("347584:50:41"), stk::formatTime(time, stk::TIMEFORMAT_HMS));
  CPPUNIT_ASSERT_EQUAL(std::string("347584:50:41.400"), stk::formatTime(time, stk::TIMEFORMAT_HMS | stk::TIMEFORMAT_MILLIS));
  CPPUNIT_ASSERT_EQUAL(std::string("1251305441"), stk::formatTime(time, stk::TIMEFORMAT_SECONDS));
  CPPUNIT_ASSERT_EQUAL(std::string("1251305441.400"), stk::formatTime(time, stk::TIMEFORMAT_SECONDS | stk::TIMEFORMAT_MILLIS));
  CPPUNIT_ASSERT_EQUAL(std::string("1.25131e+09"), stk::formatTime(time, stk::TIMEFORMAT_NONE));

  time = -41.399684906;
  
  CPPUNIT_ASSERT_EQUAL(std::string("-41"), stk::formatTime(time, stk::TIMEFORMAT_HMS));
  CPPUNIT_ASSERT_EQUAL(std::string("-41.400"), stk::formatTime(time, stk::TIMEFORMAT_HMS | stk::TIMEFORMAT_MILLIS));
  CPPUNIT_ASSERT_EQUAL(std::string("-41"), stk::formatTime(time, stk::TIMEFORMAT_SECONDS));
  CPPUNIT_ASSERT_EQUAL(std::string("-41.400"), stk::formatTime(time, stk::TIMEFORMAT_SECONDS | stk::TIMEFORMAT_MILLIS));
  CPPUNIT_ASSERT_EQUAL(std::string("-41.3997"), stk::formatTime(time, stk::TIMEFORMAT_NONE));
}



