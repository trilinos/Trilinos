#include <cppunit/TestCase.h>
#include <cppunit/extensions/HelperMacros.h>
#include <cppunit/extensions/TestFactoryRegistry.h>
#include <cppunit/ui/text/TestRunner.h>

#include <iostream>
#include <stdexcept>

#include <stk_util/environment/ReportHandler.hpp>

class UnitTestReportHandler : public CppUnit::TestCase {
private:
  CPPUNIT_TEST_SUITE(UnitTestReportHandler);
  CPPUNIT_TEST(testUnit);
  CPPUNIT_TEST_SUITE_END();

public:
  void setUp()
  {}

  void tearDown()
  {}

  void testUnit();
};


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



void UnitTestReportHandler::testUnit()
{
// Set and restore report handler.
  stk::report("This is a test", 0);

  stk::REH original_reh = stk::set_report_handler(test_report_handler);

  stk::report("This is a test", 0);
  
  stk::set_report_handler(original_reh);

  CPPUNIT_ASSERT_THROW(stk::set_report_handler(0), std::runtime_error);
  
  CPPUNIT_ASSERT_EQUAL(std::string("Message type 0: This is a test\n"), test_ostringstream().str());

  CPPUNIT_ASSERT_EQUAL(std::string("Test.cpp"), stk::source_relative_path("/src/Test.cpp"));
  CPPUNIT_ASSERT_EQUAL(std::string("Test.hpp"), stk::source_relative_path("/include/Test.hpp"));
  CPPUNIT_ASSERT_EQUAL(std::string("Apps_Test.cpp"), stk::source_relative_path("/Apps_Test.cpp"));
  CPPUNIT_ASSERT_EQUAL(std::string("stk_Test.cpp"), stk::source_relative_path("/stk_Test.cpp"));
  CPPUNIT_ASSERT_EQUAL(std::string("/smile/Test.cpp"), stk::source_relative_path("/smile/Test.cpp"));
}

CPPUNIT_TEST_SUITE_REGISTRATION(UnitTestReportHandler);
