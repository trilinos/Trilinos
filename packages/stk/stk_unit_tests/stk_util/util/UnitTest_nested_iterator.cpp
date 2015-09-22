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


#include <iterator>                     // for distance
#include <list>                         // for list
#include <gtest/gtest.h>
#include <stk_util/util/nested_iterator.hpp>  // for nested_iterator
#include <stk_util/util/nested_range.hpp>  // for identity
#include <vector>                       // for vector

namespace vector_vector_int {

  typedef std::vector< std::vector<int> > nested_type;

  typedef stk::util::nested_iterator< std::vector<std::vector<int> >,
                                      std::vector<int>,
                                      stk::util::details::identity<std::vector<int> >
                                    > nested_iterator;

  typedef stk::util::nested_iterator< const std::vector<std::vector<int> >,
                                      std::vector<int>,
                                      stk::util::details::identity<std::vector<int> >
                                    > const_nested_iterator;
}

#if defined(STK_BUILT_IN_SIERRA)
TEST ( nested_iterator, vector_vector_int)
{
  using namespace vector_vector_int;

  const int OUTER = 10;
  const int INNER = 10;

  nested_type a(OUTER);

  {
    int count = 0;
    for (int i=0; i<OUTER; ++i) {
      a[i].resize(INNER);
      for (int j=0; j<INNER; ++j) {
        a[i][j] = ++count;
      }
    }
  }

  nested_iterator itr(a);
  const nested_iterator end;

  EXPECT_EQ( OUTER*INNER, std::distance(itr,end));

  {
    int count = 0;
    for(; itr != end; ++itr) {
      *itr *= 3;
      EXPECT_EQ(++count*3,*itr);
    }
  }

}

TEST ( nested_iterator, vector_vector_int_nonconst_to_const)
{
  using namespace vector_vector_int;

  const int OUTER = 10;
  const int INNER = 10;

  nested_type a(OUTER);

  {
    int count = 0;
    for (int i=0; i<OUTER; ++i) {
      a[i].resize(INNER);
      for (int j=0; j<INNER; ++j) {
        a[i][j] = ++count;
      }
    }
  }

  nested_iterator itr(a);

  const_nested_iterator const_itr = itr;

  const_nested_iterator end;

  EXPECT_EQ( OUTER*INNER, std::distance(const_itr,end));

  {
    int count = 0;
    for(; const_itr != end; ++const_itr) {
      EXPECT_EQ(++count,*const_itr);
    }
  }

}

namespace list_vector_int {

  typedef std::list< std::vector<int> > nested_type;

  typedef stk::util::nested_iterator< std::list<std::vector<int> >,
                                      std::vector<int>,
                                      stk::util::details::identity<std::vector<int> >
                                    > nested_iterator;
}

TEST ( nested_iterator, list_vector_int)
{
  using namespace list_vector_int;

  const int OUTER = 10;
  const int INNER = 10;

  nested_type a(OUTER);

  {
    int count = 0;
    for (int i=0; i<OUTER; ++i) {
      std::vector<int> tmp;
      a.push_back(tmp);
      for (int j=0; j<INNER; ++j) {
        a.back().push_back(++count);
      }
    }
  }

  nested_iterator itr(a);
  const nested_iterator end;

  EXPECT_EQ( OUTER*INNER, std::distance(itr,end));

  {
    int count = 0;
    for(; itr != end; ++itr) {
      *itr *= 3;
      EXPECT_EQ(++count*3,*itr);
    }
  }

}

#endif
