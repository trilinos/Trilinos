#include <cmath>

#include "stk_middle_mesh/newton2.hpp"
#include "gtest/gtest.h"

namespace stk {
namespace middle_mesh {
namespace impl {

TEST(Newton2, Solve)
{
  auto func = [](const double x[2], double rhs[2]) {
    rhs[0] = std::pow(x[0], 4);
    rhs[1] = std::pow(x[1] - 1, 4);
  };
  auto jacTest = [](const double x[2], utils::impl::Mat2x2<double>& jac) {
    jac(0, 0) = 4 * std::pow(x[0], 3);
    jac(0, 1) = 0;
    jac(1, 0) = 0;
    jac(1, 1) = 4 * std::pow(x[1] - 1, 3);
  };

  double x0[2] = {2, 2};

  utils::impl::Newton2 newton(1e-60, 1000);
  newton.solve(func, jacTest, x0);

  EXPECT_NEAR(x0[0], 0.0, 1e-13);
  EXPECT_NEAR(x0[1], 1.0, 1e-13);
}

} // namespace impl
} // namespace middle_mesh
} // namespace stk
