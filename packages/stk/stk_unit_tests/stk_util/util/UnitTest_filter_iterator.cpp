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


#include <stddef.h>                     // for size_t
#include <boost/iterator/filter_iterator.hpp>  // for filter_iterator
#include <iostream>                     // for basic_ostream::operator<<, etc
#include <gtest/gtest.h>
#include <vector>                       // for vector, vector<>::iterator
#include "boost/iterator/iterator_facade.hpp"  // for iterator_facade, etc


namespace vector_vector_int {

struct vec_filter {
 vec_filter(size_t sz=0) : m_sz(sz) {}

 bool operator()(const std::vector<int>& vec) const { return vec.size() > m_sz; }

 private:
  size_t m_sz;
};

typedef std::vector< std::vector<int> > nested_type;

}

TEST ( filter_iterator, vector_vector_int)
{
  using namespace vector_vector_int;

  const int OUTER = 10;

  nested_type a(OUTER);

  {
    int count = 0;
    for (int i=0; i<OUTER; ++i) {
      a[i].resize(i+1);
      for (int j=0; j<i+1; ++j) {
        a[i][j] = ++count;
      }
    }
  }

  boost::filter_iterator<vec_filter,std::vector<std::vector<int> >::iterator> itr(vec_filter(4),a.begin(),a.end());
  boost::filter_iterator<vec_filter,std::vector<std::vector<int> >::iterator> itr_end(a.end(), a.end());

  for(; itr!=itr_end; ++itr)
  {
    std::vector<int>& v = *itr;
    for(std::vector<int>::iterator vit=v.begin(),vend=v.end(); vit!=vend; ++vit) {
      std::cout<<*vit<<" ";
    }
    std::cout<<std::endl;
  }

}
