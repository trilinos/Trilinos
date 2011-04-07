/*------------------------------------------------------------------------*/
/*                 Copyright 2010 Sandia Corporation.                     */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/

#include <utility>
#include <iostream>
#include <iomanip>
#include <cstring>
#include <string>

//boost tr1 headers...
//On the sun, couldn't get '#include <memory>' to work, so we're using the boost
//form instead...
#include <boost/tr1/memory.hpp>
#include <boost/unordered_set.hpp>
#include <boost/shared_array.hpp>

#include <stk_util/util/ci_string.hpp>

#include <stk_util/unit_test_support/stk_utest_macros.hpp>

namespace {
  
char * my_strdup(const char *s) {
  return std::strcpy(new char[std::strlen(s) + 1], s);
}

}

namespace boost {

template <>
struct hash<ci_string>
{
  std::size_t operator()(const ci_string &s) const {
    std::size_t seed = 0;
    
    for(ci_string::const_iterator first = s.begin(); first != s.end(); ++first) {
      boost::hash<char> hasher;
      seed ^= hasher(std::tolower(*first)) + 0x9e3779b9 + (seed<<6) + (seed>>2);
    }
    
    return seed;
  }
};
  
} // namespace boost

STKUNIT_UNIT_TEST(UnitTestBoost, testUnit)
{
  {  
    double* d = new double;
    boost::shared_ptr<double> dptr(d);

    STKUNIT_ASSERT_EQUAL( dptr.get(), d);

    double* d2 = new double[1];
    boost::shared_array<double> dptr2(d2);

    STKUNIT_ASSERT_EQUAL( dptr2.get(), d2);
  }
  
  // Had to comment this out because boost/tr1/array.hpp is incompatible with
  // stk_utest_macros.hpp
  /*
  boost::array<double,5> my_array;

  my_array[0] = 5.0;

  STKUNIT_ASSERT_EQUAL( my_array[0], 5.0 );
  STKUNIT_ASSERT_EQUAL( my_array.size(), (boost::array<double,5>::size_type)5 );
  */

  boost::unordered_set<int> int_set;

  int_set.insert(5);

  STKUNIT_ASSERT_EQUAL( int_set.size(), (boost::unordered_set<int>::size_type)1 );

  boost::unordered_set<ci_string> ci_string_set;

  ci_string_set.insert("Test");
  std::pair<boost::unordered_set<ci_string>::iterator, bool> res = ci_string_set.insert("test");

  STKUNIT_ASSERT_EQUAL( ci_string_set.size(), (boost::unordered_set<ci_string>::size_type)1 );
  STKUNIT_ASSERT_EQUAL( res.second, false );

  ci_string s("This is a test");

  STKUNIT_ASSERT( s == "this is a test" );
  
  std::cout << s << std::endl;
  
}

