// Copyright (c) 2013, Sandia Corporation.
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
// 
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
// 
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
// 
//     * Redistributions in binary form must reproduce the above
//       copyright notice, this list of conditions and the following
//       disclaimer in the documentation and/or other materials provided
//       with the distribution.
// 
//     * Neither the name of Sandia Corporation nor the names of its
//       contributors may be used to endorse or promote products derived
//       from this software without specific prior written permission.
// 
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
// "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
// LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
// A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
// OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
// SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
// LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
// DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
// THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
// OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
// 

#include <ctype.h>                      // for tolower
#include <algorithm>                    // for max
#include <functional>                   // for equal_to
#include <iostream>                     // for ostream, size_t, cout, endl
#include <gtest/gtest.h>
#include <stk_util/util/ci_string.hpp>  // for ci_string, operator<<
#include <string>                       // for basic_string, operator==, etc
#include <utility>                      // for pair
#include "boost/functional/hash/hash.hpp"  // for hash
#include "boost/smart_ptr/shared_array.hpp"  // for shared_array
#include "boost/smart_ptr/shared_ptr.hpp"  // for shared_ptr
#include "boost/unordered/unordered_set.hpp"  // for unordered_set, etc
namespace boost { template <class T> struct hash; }

//boost tr1 headers...
//On the sun, couldn't get '#include <memory>' to work, so we're using the boost
//form instead...



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

TEST(UnitTestBoost, testUnit)
{
  {  
    double* d = new double;
    boost::shared_ptr<double> dptr(d);

    ASSERT_EQ( dptr.get(), d);

    double* d2 = new double[1];
    boost::shared_array<double> dptr2(d2);

    ASSERT_EQ( dptr2.get(), d2);
  }
  
  // Had to comment this out because boost/tr1/array.hpp is incompatible with
  // stk_utest_macros.hpp
  /*
  boost::array<double,5> my_array;

  my_array[0] = 5.0;

  ASSERT_EQ( my_array[0], 5.0 );
  ASSERT_EQ( my_array.size(), (boost::array<double,5>::size_type)5 );
  */

  boost::unordered_set<int> int_set;

  int_set.insert(5);

  ASSERT_EQ( int_set.size(), (boost::unordered_set<int>::size_type)1 );

  boost::unordered_set<ci_string> ci_string_set;

  ci_string_set.insert("Test");
  std::pair<boost::unordered_set<ci_string>::iterator, bool> res = ci_string_set.insert("test");

  ASSERT_EQ( ci_string_set.size(), (boost::unordered_set<ci_string>::size_type)1 );
  ASSERT_EQ( res.second, false );

  ci_string s("This is a test");

  ASSERT_TRUE( s == "this is a test" );
  
  std::cout << s << std::endl;
  
}

