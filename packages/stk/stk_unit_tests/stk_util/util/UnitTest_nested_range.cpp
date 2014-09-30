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


#include <stddef.h>                     // for NULL
#include <gtest/gtest.h>
#include <stk_util/util/nested_range.hpp>  // for nested_range
#include <utility>                      // for pair, make_pair
#include <vector>                       // for vector
#include "boost/optional/optional.hpp"  // for operator!=


namespace vector_vector_int {

typedef std::vector< std::vector<int> > nested_type;
typedef std::vector< std::pair<int*,int*> > nested_pair_type;

typedef stk::util::nested_range< nested_type > range;
typedef stk::util::nested_range< nested_pair_type > pair_range;


}
/// srk 12/20/12 - these tests seem to hang on boost 1.50 / Trilinos build
#if defined(STK_BUILT_IN_SIERRA)

TEST ( nested_range, basic)
{
  using namespace vector_vector_int;

  const int OUTER = 10;
  const int INNER = 10;

  nested_type a(OUTER);
  nested_pair_type ap(OUTER);

  {
    int count = 0;
    for (int i=0; i<OUTER; ++i) {
      a[i].resize(INNER);
      int* int_ptr = NULL;
      ap[i] = std::make_pair(int_ptr,int_ptr);
      for (int j=0; j<INNER; ++j) {
        a[i][j] = ++count;
      }
    }
  }

  range rng(a);
rng.begin();
//  BOOST_FOREACH(int i, rng) {
//   std::cout<<i<<std::endl;
//  }

  pair_range prng(ap);
prng.begin();
//  BOOST_FOREACH(int i, prng) {
//   std::cout<<i<<std::endl;
//  }

}

//namespace vector_vector_int {
//
//struct to_inner_range {
//  typedef std::vector<int> result_type;
//  typedef std::vector<int>* value_type;
//
//  result_type& operator()(value_type& r) const { return *r; }
//  const result_type& operator()(const value_type& r) const { return *r; }
//};
//
//typedef std::vector< std::vector<int>* > nested_ptr_type;
//
//typedef stk::util::nested_range< nested_ptr_type, std::vector<int>, to_inner_range > ptr_range;
//
//}
//
//TEST ( nested_range, nested_ptr)
//{
//  using namespace vector_vector_int;
//
//  const int OUTER = 10;
//  const int INNER = 10;
//
//  nested_ptr_type a(OUTER);
//
//  {
//    int count = 0;
//    for (int i=0; i<OUTER; ++i) {
//      a[i] = new std::vector<int>(INNER);
//      for (int j=0; j<INNER; ++j) {
//        (*a[i])[j] = ++count;
//      }
//    }
//  }
//
//  BOOST_FOREACH(std::vector<int>* vecptr,a) {
//    std::vector<int>& vec = *vecptr;
//    BOOST_FOREACH(int i,vec) {
//      std::cout<<i<<std::endl;
//    }
//  }
//
//  ptr_range rng(a,to_inner_range());
//  BOOST_FOREACH(int i, rng) {
//   std::cout<<i<<std::endl;
//  }
//
//}

#endif
