#include "stk_middle_mesh/for_reverser.hpp"
#include "gtest/gtest.h"

namespace stk {
namespace middle_mesh {
namespace impl {

TEST(ForReverser, Forward)
{
  std::vector<int> vals(3, 0);
  utils::impl::ForReverser r(vals.size());

  for (std::size_t i = r.init(); r.condition(i); r.increment(i))
  {
    vals[i] += 1;
  }

  for (int i = 0; i < 3; ++i)
    EXPECT_EQ(vals[i], 1);
}

TEST(ForReverser, Reversed)
{
  std::vector<int> vals(3, 0);
  utils::impl::ForReverser r(vals.size(), true);

  for (std::size_t i = r.init(); r.condition(i); r.increment(i))
  {
    vals[i] += 1;
  }

  for (int i = 0; i < 3; ++i)
    EXPECT_EQ(vals[i], 1);
}

} // namespace impl
} // namespace middle_mesh
} // namespace stk
