#ifndef stk_util_parallel_TreeReductionIndexer
#define stk_util_parallel_TreeReductionIndexer

#include "stk_util/util/ReportHandler.hpp"

#include <cmath>
#include <vector>

namespace stk {
namespace impl {

// This class computes the indexing needed to do tree
// reductions.  See the diagram below, for a tree of width
// 12 where the width of the tree reduces by a factor of 3
// on each level (number of children)
//               ^
//               | x
//               | x                 x
// Depth (level) | x     x     x     x
//               | x x x x x x x x x x x x
//               ------------------>
//                 Width (rank)
class TreeReductionIndexer
{
  public:
    TreeReductionIndexer(int treeWidth, int numChildren) :
      m_treeWidth(treeWidth),
      m_numChildren(numChildren),
      m_treeDepth(compute_tree_depth())
    {}

    int get_tree_depth() const { return m_treeDepth; }

    int get_tree_width() const { return m_treeWidth; }

    int get_num_children() const { return m_numChildren; }

    int get_parent_rank(int level, int myrank) const
    {
      check_level(level);
      check_rank(myrank);
      STK_ThrowRequireMsg(level < m_treeDepth, "rank on top level has no parent");

      int rankOnLevel = myrank / int_pow(m_numChildren, level);
      return get_rank_on_level_zero(level, rankOnLevel);
    }

    int get_child_ranks(int level, int myrank, std::vector<int>& childRanks) const
    {
      check_level(level);
      check_rank(myrank);
      STK_ThrowRequireMsg(level < m_treeDepth - 1, "root rank is not anyones child");
      STK_ThrowRequireMsg(is_rank_on_level(level+1, myrank), "rank is not present above given level, cannot compute children");

      int treeWidthOnLevel = get_rank(level, get_parent_rank(level, m_treeWidth-1)) + 1;
      int rankOnLevel = get_rank(level, myrank);
      int numChildren = std::min(m_numChildren - 1, treeWidthOnLevel - rankOnLevel - 1);
      childRanks.resize(numChildren);
      for (int i=0; i < numChildren; ++i)
      {
        childRanks[i] = get_rank_on_level_zero(level, rankOnLevel + i + 1);
      }

      return numChildren;
    }

    bool is_rank_on_level(int level, int rankOnLevelZero) const
    {
      return rankOnLevelZero % int_pow(m_numChildren, level) == 0;
    }
  private:

    int get_rank_on_level_zero(int level, int rankOnLevel) const
    {
      return rankOnLevel * int_pow(m_numChildren, level);
    }

    int get_rank(int level, int rankOnLevelZero) const
    {
      STK_ThrowRequireMsg(level >= 0, "Level must be positive");
      STK_ThrowRequireMsg(rankOnLevelZero >= 0, "rank on level zero must be positive");

      return rankOnLevelZero / int_pow(m_numChildren, level);
    }

    int int_pow(int base, int exponent) const
    {
      if (exponent == 0)
        return 1;
      else if (exponent == 1)
        return base;
      else
      {
        int tmp = int_pow(base, exponent/2);
        return exponent % 2 == 0 ? tmp*tmp : base * tmp*tmp;
      }
    }

    int compute_tree_depth() const
    {
      STK_ThrowRequireMsg(m_numChildren > 1, "Tree reduction requires number of children greater than one");
      
      int width = m_treeWidth;
      int depth = 1;
      while (width > 1)
      {
        int extra = (width % m_numChildren == 0) ? 0 : 1;
        width = width / m_numChildren + extra;
        depth++;
      }

      return depth;
    }

    void check_level(int level) const
    {
      STK_ThrowRequireMsg(level >= 0 && level < m_treeDepth, "must have 0 <= level < tree depth");
    }

    void check_rank(int rank) const
    {
      STK_ThrowRequireMsg(rank >= 0 && rank < m_treeWidth, "must have 0 <= rank < tree width");

    }

    int m_treeWidth;
    int m_numChildren;
    int m_treeDepth;

};

}
}

#endif