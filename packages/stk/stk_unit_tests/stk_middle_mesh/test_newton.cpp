#include "stk_middle_mesh/newton.hpp"
#include "gtest/gtest.h"

namespace stk {
namespace middle_mesh {
namespace impl {

TEST(Newton, Solve)
{
  auto func = [](const std::vector<double>& x, std::vector<double>& rhs) {
    rhs[0] = std::pow(x[0], 4);
    rhs[1] = std::pow(x[1] - 1, 4);
  };
  auto jacTest = [](const std::vector<double>& x, utils::impl::Matrix<double>& jac) {
    jac(0, 0) = 4 * std::pow(x[0], 3);
    jac(0, 1) = 0;
    jac(1, 0) = 0;
    jac(1, 1) = 4 * std::pow(x[1] - 1, 3);
  };

  std::vector<double> x0 = {2, 2};

  utils::impl::Newton newton(2, 1e-60, 1000);
  int retVal = newton.solve(func, jacTest, x0);

  EXPECT_EQ(retVal, 0);
  EXPECT_NEAR(x0[0], 0.0, 1e-13);
  EXPECT_NEAR(x0[1], 1.0, 1e-13);
}
} // namespace impl
} // namespace middle_mesh
} // namespace stk
