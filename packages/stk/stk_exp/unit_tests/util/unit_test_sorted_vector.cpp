#include <gtest/gtest.h>

#include <sierra/util/sorted_vector.hpp>
#include <sierra/util/algorithms.hpp>

#include <boost/foreach.hpp>

#include <iostream>
#include <functional>
#include <algorithm>
#include <iterator>

using namespace sierra::util;

TEST( sorted_vector, not_unique)
{
  sorted_vector<int> a;

  EXPECT_EQ( a.size(), 0u);
  EXPECT_EQ( a.capacity(), 0u);

  a.reserve(10);

  EXPECT_EQ( a.capacity(), 10u);

  // insert in reverse order
  for (int i=10; i > 0; --i) {
    a.insert(i);
  }

  EXPECT_TRUE( sorted( a.begin(), a.end()) );

  sorted_vector<int> b(a);
  sorted_vector<int> c;
  c = a;
  EXPECT_TRUE( std::equal( a.begin(), a.end(), b.begin()));
  EXPECT_TRUE( std::equal( a.begin(), a.end(), c.begin()));

  a.clear();
  b.clear();
  c.clear();

  for( int i=100; i>0; i-=3) {
    a.insert(i);
    b.insert(b.end(),i);
  }
  EXPECT_TRUE( sorted( a.begin(), a.end()) );
  EXPECT_TRUE( sorted( b.begin(), b.end()) );

}
