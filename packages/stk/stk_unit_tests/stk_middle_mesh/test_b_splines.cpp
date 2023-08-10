#include "gtest/gtest.h"

#include "stk_middle_mesh/b_splines.hpp"

namespace stk {
namespace middle_mesh {
namespace impl {

TEST(CubicBSplines, Values)
{
  utils::impl::CubicBSplines bsplines;
  std::array<double, 4> vals;

  bsplines.eval(0, vals);
  EXPECT_NEAR(vals[0], 1.0 / 6, 1e-13);
  EXPECT_NEAR(vals[1], 4.0 / 6, 1e-13);
  EXPECT_NEAR(vals[2], 1.0 / 6, 1e-13);
  EXPECT_NEAR(vals[3], 0, 1e-13);

  bsplines.eval(1, vals);
  EXPECT_NEAR(vals[0], 0.0 / 6, 1e-13);
  EXPECT_NEAR(vals[1], 1.0 / 6, 1e-13);
  EXPECT_NEAR(vals[2], 4.0 / 6, 1e-13);
  EXPECT_NEAR(vals[3], 1.0 / 6, 1e-13);

  bsplines.eval(0.5, vals);
  EXPECT_NEAR(vals[0], 1.0 / 48, 1e-13);
  EXPECT_NEAR(vals[1], (-9.0 / 8 + 4) / 6, 1e-13);
  EXPECT_NEAR(vals[2], (15.0 / 8 + 1) / 6, 1e-13);
  EXPECT_NEAR(vals[3], 1.0 / 48, 1e-13);
}
} // namespace impl
} // namespace middle_mesh
} // namespace stk
