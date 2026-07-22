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
#include "stk_util/util/MCSR.hpp"
#include <vector>


TEST( MCSR, basic)
{
  constexpr unsigned numRows = 4;
  constexpr int invalidItem = -1;
  stk::util::MCSR<int> mcsr(numRows, invalidItem);

  for(unsigned i=0; i<numRows; ++i) {
    EXPECT_EQ(0u,  mcsr.size(i));
    EXPECT_EQ(mcsr.begin(i), mcsr.end(i));
    EXPECT_TRUE(mcsr.items(i).empty());
  }

#ifndef NDEBUG
  EXPECT_ANY_THROW(mcsr.size(numRows));
  EXPECT_ANY_THROW(mcsr.begin(numRows));
#endif
}

TEST(MCSR, find_sorted_insertion_index_empty)
{
  std::vector<int> items;
  stk::util::IndexRange indices(0,0);
  const int item = 42;
  EXPECT_EQ(0, stk::util::find_sorted_insertion_index(items, indices, item));
}

TEST(MCSR, addItem)
{
  constexpr unsigned numRows = 3;
  constexpr int invalidItem = -1;
  stk::util::MCSR<int> mcsr(numRows, invalidItem);

  unsigned row = 0;
  int item = 1;

  EXPECT_TRUE(mcsr.add_item(row, item));
  EXPECT_EQ(1u, mcsr.size(row));
  ASSERT_EQ(1u, std::distance(mcsr.begin(row), mcsr.end(row)));
  EXPECT_EQ(1u, mcsr.items(row).size());
  EXPECT_EQ(item, mcsr.begin(row)[0]);
}

TEST(MCSR, addTwoSeparateItems)
{
  constexpr unsigned numRows = 3;
  constexpr int invalidItem = -1;
  stk::util::MCSR<int> mcsr(numRows, invalidItem);

  unsigned row0 = 0, row1 = 1;
  int item1 = 1, item11 = 11;
  EXPECT_TRUE(mcsr.add_item(row0, item1));
  EXPECT_TRUE(mcsr.add_item(row1, item11));

  EXPECT_EQ(1u, mcsr.size(row0));
  ASSERT_EQ(1u, std::distance(mcsr.begin(row0), mcsr.end(row0)));
  EXPECT_EQ(item1, mcsr.begin(row0)[0]);

  EXPECT_EQ(1u, mcsr.size(row1));
  ASSERT_EQ(1u, std::distance(mcsr.begin(row1), mcsr.end(row1)));
  EXPECT_EQ(item11, mcsr.begin(row1)[0]);
}

TEST(MCSR, addTwoSeparateItemsInReverseOrder)
{
  constexpr unsigned numRows = 3;
  constexpr int invalidItem = -1;
  stk::util::MCSR<int> mcsr(numRows, invalidItem);

  unsigned row0 = 0, row1 = 1;
  int item1 = 1, item11 = 11;
  EXPECT_TRUE(mcsr.add_item(row1, item11));
  EXPECT_TRUE(mcsr.add_item(row0, item1));

  EXPECT_EQ(1u, mcsr.size(row0));
  ASSERT_EQ(1u, std::distance(mcsr.begin(row0), mcsr.end(row0)));
  EXPECT_EQ(item1, mcsr.begin(row0)[0]);

  EXPECT_EQ(1u, mcsr.size(row1));
  ASSERT_EQ(1u, std::distance(mcsr.begin(row1), mcsr.end(row1)));
  EXPECT_EQ(item11, mcsr.begin(row1)[0]);
}

TEST(MCSR, addItemsOutOfOrderForSingleRow)
{
  constexpr unsigned numRows = 3;
  constexpr int invalidItem = -1;
  stk::util::MCSR<int> mcsr(numRows, invalidItem);

  unsigned row0 = 0;

  std::vector<int> items = {5, 9, 7};

  for(int item : items) {
    EXPECT_TRUE(mcsr.add_item(row0, item));
  }

  EXPECT_EQ(3u, mcsr.size(row0));
  ASSERT_EQ(3u, std::distance(mcsr.begin(row0), mcsr.end(row0)));
  const int* rowItems = mcsr.begin(row0);
  std::sort(items.begin(), items.end());
  for(unsigned i=0; i<3u; ++i) {
    EXPECT_EQ(items[i], rowItems[i]);
  }
}

TEST(MCSR, addDuplicates_noOp)
{
  constexpr unsigned numRows = 3;
  constexpr int invalidItem = -1;
  stk::util::MCSR<int> mcsr(numRows, invalidItem);

  unsigned row0 = 0;

  std::vector<int> items = {5, 9, 7};

  for(int item : items) {
    EXPECT_TRUE(mcsr.add_item(row0, item));
  }

  EXPECT_FALSE(mcsr.add_item(row0, items[1]));

  EXPECT_EQ(3u, mcsr.size(row0));
  ASSERT_EQ(3u, std::distance(mcsr.begin(row0), mcsr.end(row0)));
  EXPECT_EQ(3u, mcsr.items(row0).size());
  const int* rowItems = mcsr.begin(row0);
  std::sort(items.begin(), items.end());
  for(unsigned i=0; i<3u; ++i) {
    EXPECT_EQ(items[i], rowItems[i]);
  }
}

TEST(MCSR, addTwoSeparateItemsThenAnotherForFirstRow)
{
  constexpr unsigned numRows = 3;
  constexpr int invalidItem = -1;
  stk::util::MCSR<int> mcsr(numRows, invalidItem);

  unsigned row0 = 0, row1 = 1;

  int item1 = 1, item11 = 11;
  EXPECT_TRUE(mcsr.add_item(row0, item1));
  EXPECT_TRUE(mcsr.add_item(row1, item11));

  int item2 = 2;
  EXPECT_TRUE(mcsr.add_item(row0, item2));

  EXPECT_EQ(2u, mcsr.size(row0));
  ASSERT_EQ(2u, std::distance(mcsr.begin(row0), mcsr.end(row0)));
  EXPECT_EQ(item1, mcsr.begin(row0)[0]);
  EXPECT_EQ(item2, mcsr.begin(row0)[1]);

  EXPECT_EQ(1u, mcsr.size(row1));
  ASSERT_EQ(1u, std::distance(mcsr.begin(row1), mcsr.end(row1)));
  EXPECT_EQ(item11, mcsr.begin(row1)[0]);

  EXPECT_EQ(3u, mcsr.total_num_items());
  EXPECT_EQ(1u, mcsr.num_unused_entries());

  mcsr.compress();

  EXPECT_EQ(3u, mcsr.total_num_items());
  EXPECT_EQ(0u, mcsr.num_unused_entries());
}

TEST(MCSR, removeTheOnlyConnectivity)
{
  constexpr unsigned numRows = 3;
  constexpr int invalidItem = -1;
  stk::util::MCSR<int> mcsr(numRows, invalidItem);

  unsigned row0 = 0;

  int item1 = 1;

  EXPECT_TRUE(mcsr.add_item(row0, item1));
  EXPECT_TRUE(mcsr.remove_item(row0, item1));

  EXPECT_EQ(0u, mcsr.size(row0));
  ASSERT_EQ(0u, std::distance(mcsr.begin(row0), mcsr.end(row0)));
  EXPECT_EQ(1u, mcsr.num_unused_entries());
  std::vector<int> testItems;
  mcsr.clear(0, testItems);
  EXPECT_EQ(1u, std::count_if(testItems.begin(), testItems.end(),
          [&](const int& item){return item == invalidItem;}));
}

TEST(MCSR, addItemsInWeirdOrder)
{
  constexpr unsigned numRows = 3;
  constexpr int invalidItem = -1;
  stk::util::MCSR<int> mcsr(numRows, invalidItem);

  unsigned row0 = 0, row1 = 1, row2 = 2;

  int item1 = 1, item2 = 2, item11 = 11, item12 = 12;

  EXPECT_TRUE(mcsr.add_item(row0, item1));
  EXPECT_TRUE(mcsr.add_item(row1, item11));
  EXPECT_TRUE(mcsr.add_item(row0, item2));
  EXPECT_TRUE(mcsr.add_item(row2, item12));

  EXPECT_EQ(2u, mcsr.size(row0));
  EXPECT_EQ(item1, mcsr.begin(row0)[0]);
  EXPECT_EQ(item2, mcsr.begin(row0)[1]);

  EXPECT_EQ(1u, mcsr.size(row1));
  EXPECT_EQ(item11, mcsr.begin(row1)[0]);

  EXPECT_EQ(1u, mcsr.size(row2));
  EXPECT_EQ(item12, mcsr.begin(row2)[0]);
}

TEST(MCSR, addAndRemoveMultipleItems)
{
  constexpr unsigned numRows = 3;
  constexpr int invalidItem = -1;
  stk::util::MCSR<int> mcsr(numRows, invalidItem);

  unsigned row0 = 0, row1 = 1;

  int entity1(1), entity2(2), entity11(11), entity12(12);

  EXPECT_TRUE(mcsr.add_item(row0, entity1));
  EXPECT_TRUE(mcsr.add_item(row1, entity11));
  EXPECT_TRUE(mcsr.add_item(row0, entity2));
  EXPECT_TRUE(mcsr.add_item(row1, entity12));
  EXPECT_EQ(2u, mcsr.num_unused_entries());

  EXPECT_TRUE(mcsr.remove_item(row0, entity2));
  EXPECT_EQ(1u, mcsr.size(row0));
  EXPECT_EQ(entity1, mcsr.begin(row0)[0]);

  EXPECT_TRUE(mcsr.remove_item(row1, entity11));
  EXPECT_EQ(1u, mcsr.size(row1));
  EXPECT_EQ(entity12, mcsr.begin(1)[0]);

  EXPECT_TRUE(mcsr.remove_item(row1, entity12));
  EXPECT_EQ(0u, mcsr.size(row1));
  EXPECT_EQ(mcsr.begin(row1), mcsr.end(row1));

  EXPECT_FALSE(mcsr.remove_item(row0, entity2));
  EXPECT_EQ(1u, mcsr.size(row0));
  EXPECT_EQ(entity1, mcsr.begin(row0)[0]);

  EXPECT_EQ(5u, mcsr.num_unused_entries());
  std::vector<int> testItems;
  mcsr.clear(0, testItems);
  EXPECT_EQ(5u, std::count_if(testItems.begin(), testItems.end(),
          [&](const int& item){return item == invalidItem;}));
}

TEST(MCSR, removeAllItemsForRow)
{
  constexpr unsigned numRows = 3;
  constexpr int invalidItem = -1;
  stk::util::MCSR<int> mcsr(numRows, invalidItem);

  unsigned row0 = 0, row1 = 1;

  int entity1(1), entity2(2), entity11(11), entity12(12);

  EXPECT_TRUE(mcsr.add_item(row0, entity1));
  EXPECT_TRUE(mcsr.add_item(row1, entity11));
  EXPECT_TRUE(mcsr.add_item(row0, entity2));
  EXPECT_TRUE(mcsr.add_item(row1, entity12));

  EXPECT_TRUE(mcsr.remove_items(row0));
  EXPECT_EQ(0u, mcsr.size(row0));
  EXPECT_EQ(mcsr.begin(0), mcsr.end(0));

  EXPECT_EQ(2u, mcsr.size(row1));
  EXPECT_EQ(entity11, mcsr.begin(row1)[0]);
  EXPECT_EQ(entity12, mcsr.begin(row1)[1]);
}

TEST(MCSR, makeHoleThenFillHoleWithoutRaisingCapacity)
{
  constexpr unsigned numRows = 3;
  constexpr int invalidItem = -1;
  stk::util::MCSR<int> mcsr(numRows, invalidItem);

  unsigned row0 = 0, row1 = 1, row2 = 2;

  int entity1(1), entity2(2), entity11(11), entity12(12);
  int entity21(21), entity22(22);

  EXPECT_TRUE(mcsr.add_item(row0, entity1));
  EXPECT_TRUE(mcsr.add_item(row0, entity2));
  EXPECT_TRUE(mcsr.add_item(row1, entity11));
  EXPECT_TRUE(mcsr.add_item(row1, entity12));
  EXPECT_TRUE(mcsr.add_item(row2, entity21));
  EXPECT_TRUE(mcsr.add_item(row2, entity22));

  size_t totalCapacity = mcsr.total_capacity();
  size_t totalNumConnectivity = mcsr.total_num_items();
  size_t numUnused = mcsr.num_unused_entries();

  EXPECT_TRUE(mcsr.remove_items(row1));
  EXPECT_EQ(0u, mcsr.size(row1));
  EXPECT_EQ(mcsr.begin(row1), mcsr.end(row1));

  EXPECT_TRUE(mcsr.add_item(row1, entity11));
  EXPECT_TRUE(mcsr.add_item(row1, entity12));

  EXPECT_TRUE(totalCapacity == mcsr.total_capacity());
  EXPECT_TRUE(totalNumConnectivity == mcsr.total_num_items());
  EXPECT_TRUE(numUnused == mcsr.num_unused_entries());
}

TEST(MCSR, removeItemsIf)
{
  constexpr unsigned numRows = 3;
  constexpr int invalidItem = -1;
  stk::util::MCSR<int> mcsr(numRows, invalidItem);

  unsigned row0 = 0, row1 = 1, row2 = 2;

  int entity1(1), entity2(2), entity11(11), entity12(12), entity13(13), entity21(21), entity22(22);

  EXPECT_TRUE(mcsr.add_item(row0, entity1));
  EXPECT_TRUE(mcsr.add_item(row0, entity2));
  EXPECT_TRUE(mcsr.add_item(row1, entity11));
  EXPECT_TRUE(mcsr.add_item(row1, entity12));
  EXPECT_TRUE(mcsr.add_item(row1, entity13));
  EXPECT_TRUE(mcsr.add_item(row2, entity21));
  EXPECT_TRUE(mcsr.add_item(row2, entity22));

  EXPECT_TRUE(mcsr.remove_items_if(row1, [&](const int& i) { return i==entity12;}));
  EXPECT_EQ(2u, mcsr.size(row1));
  EXPECT_EQ(entity11, mcsr.begin(row1)[0]);
  EXPECT_EQ(entity13, mcsr.begin(row1)[1]);
}

TEST(MCSR, removeItemsGreaterEqual)
{
  constexpr unsigned numRows = 3;
  constexpr int invalidItem = -1;
  stk::util::MCSR<int> mcsr(numRows, invalidItem);

  unsigned row0 = 0, row1 = 1, row2 = 2;

  int entity1(1), entity2(2), entity11(11), entity12(12), entity13(13), entity21(21), entity22(22);

  EXPECT_TRUE(mcsr.add_item(row0, entity1));
  EXPECT_TRUE(mcsr.add_item(row0, entity2));
  EXPECT_TRUE(mcsr.add_item(row1, entity11));
  EXPECT_TRUE(mcsr.add_item(row1, entity12));
  EXPECT_TRUE(mcsr.add_item(row1, entity13));
  EXPECT_TRUE(mcsr.add_item(row2, entity21));
  EXPECT_TRUE(mcsr.add_item(row2, entity22));

  EXPECT_TRUE(mcsr.remove_items_if(row1, [&](const int& i) { return i >= entity12; }));
  EXPECT_EQ(1u, mcsr.size(row1));
  EXPECT_EQ(entity11, mcsr.begin(row1)[0]);
}

