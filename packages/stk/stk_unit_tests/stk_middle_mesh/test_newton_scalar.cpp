#include <cmath>

#include "stk_middle_mesh/newton_scalar.hpp"
#include "gtest/gtest.h"

namespace stk {
namespace middle_mesh {
namespace impl {

TEST(NewtonScalar, Solve)
{
  auto func = [](double x) { return std::pow(x, 4); };
  auto jac  = [](double x) { return 4 * std::pow(x, 3); };

  double x0 = 2;

  utils::impl::NewtonScalar newton(1e-60, 1000);
  newton.solve(func, jac, x0);

  EXPECT_NEAR(x0, 0.0, 1e-13);
}

} // namespace impl
} // namespace middle_mesh
} // namespace stk
