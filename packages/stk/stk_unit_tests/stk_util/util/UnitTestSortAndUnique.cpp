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
#include "stk_util/util/SortAndUnique.hpp"  // for insert_keep_sorted_and_unique
#include <vector>                           // for vector


TEST( SortAndUnique, insert_keep_sorted_and_unique_overlapping_vecs)
{
  std::vector<int> sortedVec = {0, 2, 4, 6};
  std::vector<int> itemsToInsert = {5, 7, 9};

  std::vector<int> expectedResult = {0, 2, 4, 5, 6, 7, 9};

  stk::util::insert_keep_sorted_and_unique(itemsToInsert, sortedVec);

  EXPECT_EQ(expectedResult, sortedVec);
}

TEST( SortAndUnique, insert_keep_sorted_and_unique_overlapping_vecs_high_capacity)
{
  std::vector<int> sortedVec = {0, 2, 4, 6};
  sortedVec.reserve(16);
  std::vector<int> itemsToInsert = {5, 7, 9};

  std::vector<int> expectedResult = {0, 2, 4, 5, 6, 7, 9};

  stk::util::insert_keep_sorted_and_unique(itemsToInsert, sortedVec);

  EXPECT_EQ(expectedResult, sortedVec);
}

TEST( SortAndUnique, insert_keep_sorted_and_unique_vec_items_at_beginning)
{
  std::vector<int> sortedVec = {2, 4, 6};
  std::vector<int> itemsToInsert = {0, 1};

  std::vector<int> expectedResult = {0, 1, 2, 4, 6};

  stk::util::insert_keep_sorted_and_unique(itemsToInsert, sortedVec);

  EXPECT_EQ(expectedResult, sortedVec);
}

TEST( SortAndUnique, insert_keep_sorted_and_unique_vec_items_at_beginning_high_capacity)
{
  std::vector<int> sortedVec = {2, 4, 6};
  sortedVec.reserve(16);
  std::vector<int> itemsToInsert = {0, 1};

  std::vector<int> expectedResult = {0, 1, 2, 4, 6};

  stk::util::insert_keep_sorted_and_unique(itemsToInsert, sortedVec);

  EXPECT_EQ(expectedResult, sortedVec);
}

TEST( SortAndUnique, insert_keep_sorted_and_unique_vec_items_already_present)
{
  std::vector<int> sortedVec = {0, 2, 4, 6};
  std::vector<int> itemsToInsert = {4, 6};

  std::vector<int> expectedResult = {0, 2, 4, 6};

  stk::util::insert_keep_sorted_and_unique(itemsToInsert, sortedVec);

  EXPECT_EQ(expectedResult, sortedVec);
}

TEST( SortAndUnique, insert_keep_sorted_and_unique_vec_items_already_present_high_capacity)
{
  std::vector<int> sortedVec = {0, 2, 4, 6};
  sortedVec.reserve(16);
  std::vector<int> itemsToInsert = {4, 6};

  std::vector<int> expectedResult = {0, 2, 4, 6};

  stk::util::insert_keep_sorted_and_unique(itemsToInsert, sortedVec);

  EXPECT_EQ(expectedResult, sortedVec);
}

TEST( SortAndUnique, insert_keep_sorted_and_unique_vec_items_empty)
{
  std::vector<int> sortedVec = {0, 2, 4, 6};
  std::vector<int> itemsToInsert;

  std::vector<int> expectedResult = {0, 2, 4, 6};

  stk::util::insert_keep_sorted_and_unique(itemsToInsert, sortedVec);

  EXPECT_EQ(expectedResult, sortedVec);
}

TEST( SortAndUnique, insert_keep_sorted_and_unique_vec_empty)
{
  std::vector<int> sortedVec;
  std::vector<int> itemsToInsert = {0, 2, 4, 6};

  std::vector<int> expectedResult = {0, 2, 4, 6};

  stk::util::insert_keep_sorted_and_unique(itemsToInsert, sortedVec);

  EXPECT_EQ(expectedResult, sortedVec);
}

struct MyTestLess {
  bool operator()(int lhs, int rhs) const { return lhs < rhs; }
};

TEST( SortAndUnique, insert_keep_sorted_and_unique_overlapping_vecs_myless)
{
  std::vector<int> sortedVec = {0, 2, 4, 6};
  std::vector<int> itemsToInsert = {5, 7, 9};

  std::vector<int> expectedResult = {0, 2, 4, 5, 6, 7, 9};

  stk::util::insert_keep_sorted_and_unique(itemsToInsert, sortedVec, MyTestLess());

  EXPECT_EQ(expectedResult, sortedVec);
}

TEST( SortAndUnique, insert_keep_sorted_and_unique_overlapping_vecs_myless_high_capacity)
{
  std::vector<int> sortedVec = {0, 2, 4, 6};
  sortedVec.reserve(16);
  std::vector<int> itemsToInsert = {5, 7, 9};

  std::vector<int> expectedResult = {0, 2, 4, 5, 6, 7, 9};

  stk::util::insert_keep_sorted_and_unique(itemsToInsert, sortedVec, MyTestLess());

  EXPECT_EQ(expectedResult, sortedVec);
}

TEST( SortAndUnique, insert_keep_sorted_and_unique_overlapping_vecs_items_at_end_myless_high_capacity)
{
  std::vector<int> sortedVec = {0, 2, 4, 6};
  sortedVec.reserve(16);
  std::vector<int> itemsToInsert = {8, 9, 10};

  std::vector<int> expectedResult = {0, 2, 4, 6, 8, 9, 10};

  stk::util::insert_keep_sorted_and_unique(itemsToInsert, sortedVec, MyTestLess());

  EXPECT_EQ(expectedResult, sortedVec);
}

#ifndef NDEBUG
TEST( SortAndUnique, insert_keep_sorted_and_unique_items_not_sorted)
{
  std::vector<int> sortedVec = {0, 4, 2, 6};
  std::vector<int> itemsToInsert = {4, 6};

  std::vector<int> expectedResult = {0, 2, 4, 6};

  EXPECT_THROW(stk::util::insert_keep_sorted_and_unique(itemsToInsert, sortedVec), std::logic_error);
}
#endif

