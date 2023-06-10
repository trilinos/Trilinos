#include "stk_middle_mesh/backtracking_line_search.hpp"
#include "gtest/gtest.h"

namespace stk {
namespace middle_mesh {
namespace impl {

TEST(LineSearch, Quadratic)
{
  auto f = [](const double x) { return (x - 1) * (x - 1); };

  double x0    = 1;
  double pk    = 0.5;
  double grad0 = 2 * (x0 - 1) * pk;

  opt::impl::BacktrackingLineSearch linesearch(0.5, 0.9);

  double alpha = linesearch.search(f, grad0);
  EXPECT_FLOAT_EQ(alpha, 1);
}

} // namespace impl
} // namespace middle_mesh
} // namespace stk
