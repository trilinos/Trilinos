 // Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
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
//     * Neither the name of NTESS nor the names of its contributors
//       may be used to endorse or promote products derived from this
//       software without specific prior written permission.
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

#include <algorithm>
#include <functional>
#include <stk_util/util/ReportHandler.hpp>

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

template<class VECTOR, class COMPARE_LESS>
bool is_sorted_and_unique(const VECTOR& vec, COMPARE_LESS cmp)
{
  for(size_t i=1; i<vec.size(); ++i) {
    if (!cmp(vec[i-1], vec[i])) {
      return false;
    }
  }
  return true;
}


template<class VECTOR>
bool is_sorted_and_unique(const VECTOR& vec)
{
  return is_sorted_and_unique(vec, std::less<typename VECTOR::value_type>());
}

template<class VECTOR>
bool insert_keep_sorted_and_unique(typename VECTOR::value_type item, VECTOR& vec)
{
  //This ThrowAssert causes a DetectHinge unit-test to fail! That needs to be debugged so that this can be turned on!
  //ThrowAssertMsg(is_sorted_and_unique(vec), "input vector must be sorted and unique");

  typename VECTOR::iterator iter = std::lower_bound(vec.begin(), vec.end(), item);
  if (iter == vec.end() || *iter != item) {
    vec.insert(iter, item);
    return true;
  }
  return false;
}

template<class VECTOR, typename COMPARE_LESS>
bool insert_keep_sorted_and_unique(typename VECTOR::value_type item, VECTOR& vec, COMPARE_LESS compare_less)
{
  //This ThrowAssert causes a DetectHinge unit-test to fail! That needs to be debugged so that this can be turned on!
  //ThrowAssertMsg(is_sorted_and_unique(vec,compare_less), "input vector must be sorted and unique");

  typename VECTOR::iterator iter = std::lower_bound(vec.begin(), vec.end(), item, compare_less);
  if (iter == vec.end() || *iter != item) {
    vec.insert(iter, item);
    return true;
  }
  return false;
}

template<class VECTOR, typename COMPARE_LESS, typename COMPARE_EQUALS>
bool insert_keep_sorted_and_unique(typename VECTOR::value_type item, VECTOR& vec, COMPARE_LESS compare_less, COMPARE_EQUALS compare_equals)
{
  ThrowAssertMsg(is_sorted_and_unique(vec,compare_less), "input vector must be sorted and unique");

  typename VECTOR::iterator iter = std::lower_bound(vec.begin(), vec.end(), item, compare_less);
  if (iter == vec.end() || !compare_equals(*iter, item)) {
    vec.insert(iter, item);
    return true;
  }
  return false;
}

namespace impl {

template<class ITEM, class COMPARE_LESS>
struct CompareEqual {
  CompareEqual(COMPARE_LESS compareLess) : cmp(compareLess)
  {}

  bool operator()(const ITEM& lhs, const ITEM& rhs) const
  { return (!cmp(lhs,rhs) && !cmp(rhs,lhs)); }

  COMPARE_LESS cmp;
};

} //namespace impl

template<class VECTOR, class COMPARE_LESS>
bool insert_keep_sorted_and_unique(const VECTOR& sortedItemsToInsert, VECTOR& sortedVec, COMPARE_LESS cmp)
{
  if (sortedItemsToInsert.empty()) {
    return false;
  }

  ThrowAssertMsg(is_sorted_and_unique(sortedItemsToInsert, cmp), "input items must be sorted and unique");
  ThrowAssertMsg(is_sorted_and_unique(sortedVec, cmp), "input vector must be sorted and unique");

  const size_t oldLength = sortedVec.size();

  if (oldLength > 0) {
    typename VECTOR::value_type firstNewItem = sortedItemsToInsert[0];
    typename VECTOR::value_type lastNewItem = sortedItemsToInsert.back();

    typename VECTOR::iterator startOfNew = std::lower_bound(sortedVec.begin(), sortedVec.end(),
                                                            firstNewItem, cmp);
    const unsigned num = startOfNew-sortedVec.begin();

    sortedVec.insert(startOfNew, sortedItemsToInsert.begin(), sortedItemsToInsert.end());
    startOfNew = sortedVec.begin()+num;

    typename VECTOR::iterator endOfNew = startOfNew+sortedItemsToInsert.size();
    typename VECTOR::iterator endOfMerge = std::lower_bound(endOfNew, sortedVec.end(), lastNewItem, cmp);

    std::inplace_merge(startOfNew, endOfNew, endOfMerge);

    typename VECTOR::iterator iter = std::unique(startOfNew, sortedVec.end(),
                           impl::CompareEqual<typename VECTOR::value_type,COMPARE_LESS>(cmp));
    sortedVec.resize(iter - sortedVec.begin());
  }
  else {
    sortedVec = sortedItemsToInsert;
  }

  return sortedVec.size() > oldLength;
}

template<class VECTOR>
bool insert_keep_sorted_and_unique(const VECTOR& sortedItemsToInsert, VECTOR& sortedVec)
{
  return insert_keep_sorted_and_unique(sortedItemsToInsert, sortedVec, std::less<typename VECTOR::value_type>());
}

} //namespace util
} //namespace stk

#endif
