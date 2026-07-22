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

#include "gtest/gtest.h"
#include "stk_util/util/SortedIntersect.hpp"
#include "stk_util/util/SortedContains.hpp"
#include <vector>

TEST( SortedContains, bothEmpty)
{
  std::vector<int> sortedSuperset;
  std::vector<int> sortedSubset;

  EXPECT_TRUE(stk::util::sorted_contains_all(sortedSuperset.begin(), sortedSuperset.end(),
                                              sortedSubset.begin(), sortedSubset.end()));
}

TEST( SortedContains, firstEmpty)
{
  std::vector<int> sortedSuperset;
  std::vector<int> sortedSubset = {1, 2, 3};

  EXPECT_FALSE(stk::util::sorted_contains_all(sortedSuperset.begin(), sortedSuperset.end(),
                                              sortedSubset.begin(), sortedSubset.end()));
}

TEST( SortedContains, secondEmpty)
{
  std::vector<int> sortedSuperset = {1, 2, 3};
  std::vector<int> sortedSubset;

  EXPECT_TRUE(stk::util::sorted_contains_all(sortedSuperset.begin(), sortedSuperset.end(),
                                              sortedSubset.begin(), sortedSubset.end()));
}

TEST( SortedContains, firstFoundInSecond)
{
  std::vector<int> sortedSuperset = {1, 2};
  std::vector<int> sortedSubset = {0, 1, 3, 5};

  EXPECT_FALSE(stk::util::sorted_contains_all(sortedSuperset.begin(), sortedSuperset.end(),
                                             sortedSubset.begin(), sortedSubset.end()));
}

TEST( SortedContains, firstEndSmallerThanSecondBegin)
{
  std::vector<int> sortedSuperset = {1, 2};
  std::vector<int> sortedSubset = {3, 4, 5, 6};

  EXPECT_FALSE(stk::util::sorted_contains_all(sortedSuperset.begin(), sortedSuperset.end(),
                                              sortedSubset.begin(), sortedSubset.end()));
}

TEST( SortedContains, firstBeginBiggerThanSecondEnd)
{
  std::vector<int> sortedSuperset = {3, 4, 5, 6};
  std::vector<int> sortedSubset = {1, 2};

  EXPECT_FALSE(stk::util::sorted_contains_all(sortedSuperset.begin(), sortedSuperset.end(),
                                              sortedSubset.begin(), sortedSubset.end()));
}

TEST( SortedContains, firstContainsSecond)
{
  std::vector<int> sortedSuperset = {1, 2, 3, 4, 5, 6};
  std::vector<int> sortedSubset = {2, 3};

  EXPECT_TRUE(stk::util::sorted_contains_all(sortedSuperset.begin(), sortedSuperset.end(),
                                              sortedSubset.begin(), sortedSubset.end()));
}

TEST( SortedContains, firstEqualsSecond)
{
  std::vector<int> sortedSuperset = {1, 2, 3, 4, 5, 6};
  std::vector<int> sortedSubset = {1, 2, 3, 4, 5, 6};

  EXPECT_TRUE(stk::util::sorted_contains_all(sortedSuperset.begin(), sortedSuperset.end(),
                                              sortedSubset.begin(), sortedSubset.end()));
}

