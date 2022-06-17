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
//

#ifndef _CONVEXGROUP_H_
#define _CONVEXGROUP_H_

#include "stk_util/util/SortAndUnique.hpp"

#include <vector>
#include <algorithm>

namespace stk {
namespace tools {
namespace impl {


template <typename PAIR, typename ID, typename COMPARE>
class ConvexGroup {
public:
  ConvexGroup(COMPARE comparator, ID getId) : m_comparator(comparator), m_idGetter(getId) {}
  ConvexGroup(COMPARE comparator, ID getId, typename PAIR::UNIT block) : m_comparator(comparator), m_idGetter(getId)
  {
    m_members.push_back(block);
  }
  ConvexGroup(COMPARE comparator, ID getId, typename PAIR::UNIT block1, typename PAIR::UNIT block2) : m_comparator(comparator), m_idGetter(getId)
  {
    m_members.push_back(block1);
    stk::util::insert_keep_sorted_and_unique(block2, m_members, m_comparator);
  }

  const std::vector<typename PAIR::UNIT>& get_members() const { return m_members; }
  COMPARE& get_comparator() { return m_comparator; }

  unsigned get_id() const {
    unsigned groupId = std::numeric_limits<unsigned>::max();

    for(const typename PAIR::UNIT& member : m_members) {
      groupId = std::min(groupId, m_idGetter(member));
    }
    return groupId;
  }

private:
  std::vector<typename PAIR::UNIT> m_members;
  COMPARE m_comparator;
  ID m_idGetter;
};

template <typename PAIR, typename ID, typename COMPARE>
using ConvexGroupVector = std::vector<ConvexGroup<PAIR,ID,COMPARE>>;

template <typename PAIR, typename ID, typename COMPARE>
int find_in_groups(const ConvexGroupVector<PAIR,ID,COMPARE>& groupVec, typename PAIR::UNIT& part, COMPARE& comparator)
{
  for(unsigned i = 0; i < groupVec.size(); i++) {
    const std::vector<typename PAIR::UNIT> members = groupVec[i].get_members();
    if(std::binary_search(members.begin(), members.end(), part, comparator)) {
      return i;
    }
  }
  return -1;
}

template <typename PAIR, typename ID, typename COMPARE>
void merge_groups(ConvexGroupVector<PAIR,ID,COMPARE>& groupVec, int groupIndex1, int groupIndex2, COMPARE& comparator)
{
  if(groupIndex1 == groupIndex2 || groupIndex1 >= (int)groupVec.size() || groupIndex2 >= (int)groupVec.size() ||
     groupIndex1 < 0 || groupIndex2 < 0) { return; }

  auto it1 = groupVec.begin()+groupIndex1;
  auto it2 = groupVec.begin()+groupIndex2;

  for(typename PAIR::UNIT& unit : it2->get_members()) {
    stk::util::insert_keep_sorted_and_unique(unit, it1->get_members(), comparator);
  }
  groupVec.erase(it2);
}

template <typename PAIR, typename ID, typename COMPARE>
void populate_convex_group(const PAIR& pair, ConvexGroupVector<PAIR,ID,COMPARE>& groupVec, ID& getId, COMPARE& comparator)
{
  const typename PAIR::UNIT& first = pair.get_first();
  const typename PAIR::UNIT& second = pair.get_second();

  int groupIndex1 = find_in_groups(groupVec, first, comparator);
  int groupIndex2 = find_in_groups(groupVec, second, comparator);

  if(pair.is_adjacent()) {
    if(groupIndex1 < 0 && groupIndex2 < 0) {
      groupVec.push_back(ConvexGroup<PAIR,ID,COMPARE>(comparator, getId, first, second));
    } else if(groupIndex1 >= 0 && groupIndex2 >= 0) {
      merge_groups(groupVec, groupIndex1, groupIndex2, comparator);
    } else if(groupIndex1 < 0 && groupIndex2 >= 0) {
      stk::util::insert_keep_sorted_and_unique(first, groupVec[groupIndex2].get_members(), comparator);
    } else {
      stk::util::insert_keep_sorted_and_unique(second, groupVec[groupIndex1].get_members(), comparator);
    }
  } else {
    if(groupIndex1 < 0) {
      groupVec.push_back(ConvexGroup<PAIR,ID,COMPARE>(comparator, getId, first));
    }
    if(groupIndex2 < 0) {
      groupVec.push_back(ConvexGroup<PAIR,ID,COMPARE>(comparator, getId, second));
    }
  }
}

template <typename PAIR, typename ID, typename COMPARE>
ConvexGroupVector<PAIR,ID,COMPARE> get_convex_groups(const std::vector<PAIR>& pairVector, ID& getId, COMPARE& comparator)
{
  ConvexGroupVector<PAIR,ID,COMPARE> groupVec;

  for(const PAIR& pair : pairVector) {
    populate_convex_group(pair, groupVec, getId, comparator);
  }

  return groupVec;
}

}
}
}

#endif
