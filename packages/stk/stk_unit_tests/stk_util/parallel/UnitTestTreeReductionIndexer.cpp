#include "stk_util/parallel/TreeReductionIndexer.hpp"
#include "gtest/gtest.h"

TEST(TreeReductionIndexer, Getters)
{
  stk::impl::TreeReductionIndexer indexer(36, 3);

  EXPECT_EQ(indexer.get_tree_width(), 36);
  EXPECT_EQ(indexer.get_tree_depth(), 5);
  EXPECT_EQ(indexer.get_num_children(), 3);
}

TEST(TreeReductionIndexer, TreeDepth2Children)
{

  std::vector<int> tree_widths = {1, 2, 3, 4, 5, 6, 7, 8};
  std::vector<int> tree_depths = {1, 2, 3, 3, 4, 4, 4, 4};

  for (size_t i=0; i < tree_widths.size(); ++i)
  {
    stk::impl::TreeReductionIndexer indexer(tree_widths[i], 2);
    EXPECT_EQ(indexer.get_tree_depth(), tree_depths[i]);
  }
}

TEST(TreeReductionIndexer, TreeDepth3Children)
{

  std::vector<int> tree_widths = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14};
  std::vector<int> tree_depths = {1, 2, 2, 3, 3, 3, 3, 3, 3, 4,  4, 4, 4, 4};

  for (size_t i=0; i < tree_widths.size(); ++i)
  {
    stk::impl::TreeReductionIndexer indexer(tree_widths[i], 3);
    EXPECT_EQ(indexer.get_tree_depth(), tree_depths[i]);
  }
}

TEST(TreeReductionIndexer, NumChildrenErrors)
{
  EXPECT_ANY_THROW(stk::impl::TreeReductionIndexer indexer(2, -1));
  EXPECT_ANY_THROW(stk::impl::TreeReductionIndexer indexer(2, 0));
  EXPECT_ANY_THROW(stk::impl::TreeReductionIndexer indexer(2, 1));
}


TEST(TreeReductionIndexer, ParentRank3Children)
{
  int treeWidth = 13;
  stk::impl::TreeReductionIndexer indexer(treeWidth, 3);
  std::vector<int> levelZeroParents = {0, 0, 0, 3, 3, 3, 6, 6, 6, 9, 9, 9, 12};
  for (int i=0; i < treeWidth; ++i)
  {
    EXPECT_EQ(indexer.get_parent_rank(1, i), levelZeroParents[i]);  
  }

  std::vector<int> levelOneParents = {0, 0, 0, 0, 0, 0, 0, 0, 0, 9, 9, 9, 9};
  for (int i=0; i < treeWidth; ++i)
  {
    EXPECT_EQ(indexer.get_parent_rank(2, i), levelOneParents[i]);  
  }

  std::vector<int> levelTwoParents(treeWidth, 0);
  for (int i=0; i < treeWidth; ++i)
  {
    EXPECT_EQ(indexer.get_parent_rank(3, i), levelTwoParents[i]);  
  }
}

TEST(TreeReductionIndexer, ParentRankErrors)
{

  int treeWidth = 13;
  stk::impl::TreeReductionIndexer indexer(treeWidth, 3);
  EXPECT_ANY_THROW(indexer.get_parent_rank(indexer.get_tree_depth() + 1, 0));
  EXPECT_ANY_THROW(indexer.get_parent_rank(-1, 0));
}

namespace {
void testChildren(stk::impl::TreeReductionIndexer& indexer, int level, std::vector< std::vector<int> >& levelChildren)
{ 
  int idx = 0;
  std::vector<int> children;
  for (int i=0; i < indexer.get_tree_width(); ++i)
  {
    if (indexer.is_rank_on_level(level+1, i))
    {
      indexer.get_child_ranks(level, i, children);
      EXPECT_EQ(children.size(), levelChildren[idx].size());
      for (size_t j=0; j < levelChildren[idx].size(); ++j)
      {
        EXPECT_EQ(children[j], levelChildren[idx][j]);
      }
      idx++;
    } else
    {
      EXPECT_ANY_THROW(indexer.get_child_ranks(level, i, children));
    }
  }  
}
}

TEST(TreeReductionIndexer, ChildRanks3Children)
{
  int treeWidth = 13;
  stk::impl::TreeReductionIndexer indexer(treeWidth, 3);

  std::vector< std::vector<int> > levelZeroChildren = { {1, 2}, {4, 5}, {7, 8}, {10, 11}, {} };
  std::vector< std::vector<int> > levelOneChildren  = { {3, 6}, {12}};

  testChildren(indexer, 0, levelZeroChildren);

  testChildren(indexer, 1, levelOneChildren);
}

TEST(TreeReductionIndexer, ChildRankErrors)
{
  int treeWidth = 13;
  stk::impl::TreeReductionIndexer indexer(treeWidth, 3);

  std::vector<int> invalidRanks = {1, 2, 4, 5, 7, 8, 10, 11};
  std::vector<int> children;
  for (auto& rank : invalidRanks)
  {
    EXPECT_ANY_THROW(indexer.get_child_ranks(rank, 1, children));
  }
}


TEST(TreeReductionIndexer, IsRankOnLevel)
{
  int treeWidth = 13;
  stk::impl::TreeReductionIndexer indexer(treeWidth, 3);

  std::vector<bool> isOnLevelOne = {true, false, false, true, false, false, true, false, false, true, false, false, true};
  std::vector<bool> isOnLevelTwo = {true, false, false, false, false, false, false, false, false, true, false, false, false};

  for (int i=0; i < treeWidth; ++i)
  {
    EXPECT_TRUE(indexer.is_rank_on_level(0, i));
    EXPECT_EQ(indexer.is_rank_on_level(1, i), isOnLevelOne[i]);
    EXPECT_EQ(indexer.is_rank_on_level(2, i), isOnLevelTwo[i]);
  }

}