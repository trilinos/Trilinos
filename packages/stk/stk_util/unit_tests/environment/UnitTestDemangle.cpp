/*------------------------------------------------------------------------*/
/*                 Copyright 2010, 2011 Sandia Corporation.                     */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/

#include <iostream>
#include <unistd.h>

#include <stk_util/environment/Demangle.hpp>

#include <stk_util/unit_test_support/stk_utest_macros.hpp>

bool
utest_demangle()
{
  return true;
}

STKUNIT_UNIT_TEST(UnitTestDemangle, UnitTest)
{
#if defined(__PGI)
  {
    std::string linux_name("ThisIsJunk");
    std::string demangled_name = stk::demangle(linux_name.c_str());
    STKUNIT_ASSERT_EQUAL(linux_name, demangled_name);
  }

  {
    std::string linux_name("bool ()");
    std::string demangled_name = stk::demangle(typeid(utest_demangle).name());
    STKUNIT_ASSERT_EQUAL(linux_name, demangled_name);
  }

  {
    typedef std::vector<double> DoubleVector;

    DoubleVector double_vector;

    std::string linux_name("std::vector<double, std::allocator<double>>");
    std::string demangled_name = stk::demangle(typeid(double_vector).name());
    STKUNIT_ASSERT_EQUAL(linux_name, demangled_name);
  }

#elif defined(__sun)
  {
    std::string linux_name("ThisIsJunk");
    std::string demangled_name = stk::demangle(linux_name.c_str());
    STKUNIT_ASSERT_EQUAL(linux_name, demangled_name);
  }

  {
    std::string linux_name("bool()");
    std::string demangled_name = stk::demangle(typeid(utest_demangle).name());
    STKUNIT_ASSERT_EQUAL(linux_name, demangled_name);
  }

  {
    typedef std::vector<double> DoubleVector;

    DoubleVector double_vector;

    std::string linux_name("std::vector<double, std::allocator<double>>");
    std::string demangled_name = stk::demangle(typeid(double_vector).name());
    STKUNIT_ASSERT_EQUAL(linux_name, demangled_name);
  }

#elif defined(__xlC__)
  {
    std::string linux_name("ThisIsJunk");
    std::string demangled_name = stk::demangle(linux_name.c_str());
//    STKUNIT_ASSERT_EQUAL((linux_name == demangled_name), true);
  }

  {
    std::string linux_name("bool ()()");
    std::string demangled_name = stk::demangle(typeid(utest_demangle).name());
//    STKUNIT_ASSERT_EQUAL(linux_name, demangled_name);
  }

  {
    typedef std::vector<double> DoubleVector;

    DoubleVector double_vector;

    std::string linux_name("std::vector<double, std::allocator<double> >");
    std::string demangled_name = stk::demangle(typeid(double_vector).name());
//    STKUNIT_ASSERT_EQUAL((linux_name == demangled_name), true);
  }
#elif defined(__linux__)
  {
    std::string linux_name("ThisIsJunk");
    std::string demangled_name = stk::demangle(linux_name.c_str());
    STKUNIT_ASSERT_EQUAL((linux_name == demangled_name), true);
  }

  {
    std::string linux_name("bool ()()");
    std::string demangled_name = stk::demangle(typeid(utest_demangle).name());
    STKUNIT_ASSERT_EQUAL(linux_name, demangled_name);
  }

  {
    typedef std::vector<double> DoubleVector;

    DoubleVector double_vector;

#ifdef _GLIBCXX_DEBUG
    std::string linux_name("__gnu_debug_def::vector<double, std::allocator<double> >");
#else
    std::string linux_name("std::vector<double, std::allocator<double> >");
#endif
    std::string demangled_name = stk::demangle(typeid(double_vector).name());
    STKUNIT_ASSERT_EQUAL((linux_name == demangled_name), true);
  }
#endif
}
