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

#ifndef stk_util_util_SortAndUnique_hpp
#define stk_util_util_SortAndUnique_hpp

#include<algorithm>

namespace stk
{
namespace util
{

template<typename VECTOR, typename COMPARE>
void sort_and_unique(VECTOR &vector, COMPARE compare )
{
    std::sort(vector.begin(), vector.end(), compare);
    auto endIter = std::unique(vector.begin(), vector.end());
    vector.resize(endIter - vector.begin());
}

template<typename VECTOR, typename COMPARE_LESS, typename COMPARE_EQUALS>
void sort_and_unique(VECTOR &vector, COMPARE_LESS compare_less, COMPARE_EQUALS compare_equals  )
{
    std::sort(vector.begin(), vector.end(), compare_less);
    auto endIter = std::unique(vector.begin(), vector.end(), compare_equals);
    vector.resize(endIter - vector.begin());
}

template<typename VECTOR>
void sort_and_unique(VECTOR &vec)
{
    sort_and_unique(vec,std::less<typename VECTOR::value_type>());
}

template<class VECTOR, typename COMPARE>
bool is_sorted_and_unique(const VECTOR& vec, COMPARE compare)
{
    bool sorted_and_unique = true;
    for(size_t i=1; i<vec.size(); ++i) {
        if (!compare(vec[i-1],vec[i])) {
            sorted_and_unique = false;
        }
    }
    return sorted_and_unique;
}

template<class VECTOR>
bool insert_keep_sorted_and_unique(typename VECTOR::value_type p, VECTOR& procs)
{
  typename VECTOR::iterator iter = std::lower_bound(procs.begin(), procs.end(), p);
  if (iter == procs.end() || *iter != p) {
    procs.insert(iter, p);
    return true;
  }
  return false;
}

} //namespace util
} //namespace stk

#endif
