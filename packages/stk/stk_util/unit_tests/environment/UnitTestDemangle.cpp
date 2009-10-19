#include <cppunit/TestCase.h>
#include <cppunit/extensions/HelperMacros.h>
#include <cppunit/extensions/TestFactoryRegistry.h>
#include <cppunit/ui/text/TestRunner.h>

#include <iostream>
#include <unistd.h>

#include <stk_util/environment/Demangle.hpp>

class UnitTestDemangle : public CppUnit::TestCase {
private:
  CPPUNIT_TEST_SUITE( UnitTestDemangle );
  CPPUNIT_TEST( testUnit );
  CPPUNIT_TEST_SUITE_END();

public:
  void setUp() {}
  void tearDown() {}
  void testUnit();
};

CPPUNIT_TEST_SUITE_REGISTRATION(UnitTestDemangle);

bool
utest_demangle()
{
  return true;
}

void
UnitTestDemangle::testUnit()
{
#if defined(__PGIC)
  {
    std::string linux_name("ThisIsJunk");
    std::string demangled_name = stk::demangle(linux_name.c_str());
    CPPUNIT_ASSERT_EQUAL(linux_name, demangled_name);
  }

  {
    std::string linux_name("bool()");
    std::string demangled_name = stk::demangle(typeid(utest_demangle).name());
    CPPUNIT_ASSERT_EQUAL(linux_name, demangled_name);
  }

  {
    typedef std::vector<double> DoubleVector;

    DoubleVector double_vector;

    std::string linux_name("std::vector<double, std::allocator<double>>");
    std::string demangled_name = stk::demangle(typeid(double_vector).name());
    CPPUNIT_ASSERT_EQUAL(linux_name, demangled_name);
  }

#elif defined(__sun)
  {
    std::string linux_name("ThisIsJunk");
    std::string demangled_name = stk::demangle(linux_name.c_str());
    CPPUNIT_ASSERT_EQUAL(linux_name, demangled_name);
  }

  {
    std::string linux_name("bool()");
    std::string demangled_name = stk::demangle(typeid(utest_demangle).name());
    CPPUNIT_ASSERT_EQUAL(linux_name, demangled_name);
  }

  {
    typedef std::vector<double> DoubleVector;

    DoubleVector double_vector;

    std::string linux_name("std::vector<double, std::allocator<double>>");
    std::string demangled_name = stk::demangle(typeid(double_vector).name());
    CPPUNIT_ASSERT_EQUAL(linux_name, demangled_name);
  }

#elif defined(__linux__)
  {
    std::string linux_name("ThisIsJunk");
    std::string demangled_name = stk::demangle(linux_name.c_str());
    CPPUNIT_ASSERT_EQUAL(linux_name, demangled_name);
  }

  {
    std::string linux_name("bool ()()");
    std::string demangled_name = stk::demangle(typeid(utest_demangle).name());
    CPPUNIT_ASSERT_EQUAL(linux_name, demangled_name);
  }

  {
    typedef std::vector<double> DoubleVector;

    DoubleVector double_vector;

    std::string linux_name("std::vector<double, std::allocator<double> >");
    std::string demangled_name = stk::demangle(typeid(double_vector).name());
    CPPUNIT_ASSERT_EQUAL(linux_name, demangled_name);
  }
#endif
}
