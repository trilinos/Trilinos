#include "stk_middle_mesh/taylor_patch.hpp"
#include "gtest/gtest.h"

namespace stk {
namespace middle_mesh {
namespace impl {

namespace {
void expect_float_eq(const utils::Point& pt1, const utils::Point& pt2, double tol = -1)
{
  if (tol < 0)
  {
    EXPECT_FLOAT_EQ(pt1.x, pt2.x);
    EXPECT_FLOAT_EQ(pt1.y, pt2.y);
    EXPECT_FLOAT_EQ(pt1.z, pt2.z);
  } else
  {
    EXPECT_NEAR(pt1.x, pt2.x, tol);
    EXPECT_NEAR(pt1.y, pt2.y, tol);
    EXPECT_NEAR(pt1.z, pt2.z, tol);
  }
}

template <typename T>
void expect_float_eq(const utils::impl::Matrix<T>& a, const utils::impl::Matrix<T>& b, double tol = -1)
{
  ASSERT_EQ(a.extent0(), b.extent0());
  ASSERT_EQ(a.extent1(), b.extent1());
  if (tol < 0)
    for (int i = 0; i < a.extent0(); ++i)
      for (int j = 0; j < a.extent1(); ++j)
        EXPECT_FLOAT_EQ(a(i, j), b(i, j));
  else
    for (int i = 0; i < a.extent0(); ++i)
      for (int j = 0; j < a.extent1(); ++j)
        EXPECT_NEAR(a(i, j), b(i, j), tol);
}

} // namespace

TEST(TaylorPatch, Factorial)
{
  EXPECT_EQ(utils::impl::factorial(0), 1);
  EXPECT_EQ(utils::impl::factorial(1), 1);
  EXPECT_EQ(utils::impl::factorial(2), 2);
  EXPECT_EQ(utils::impl::factorial(3), 6);
  EXPECT_EQ(utils::impl::factorial(4), 24);
}

TEST(TaylorPatch, Plane)
{
  // test plane through z=0
  std::vector<utils::Point> pts = {utils::Point(0, 0, 0), utils::Point(1, 0, 0), utils::Point(0, 1, 0)};
  utils::impl::TaylorPatch patch;
  patch.construct_patch(pts);

  // for this case, the (u, v, w) coordinate system is aligned with (x,y,z)
  expect_float_eq(patch.eval_point(0, 0), utils::Point(0, 0));
  expect_float_eq(patch.eval_point(1, 0), utils::Point(1, 0));
  expect_float_eq(patch.eval_point(0, 1), utils::Point(0, 1));

  double derivs[2];
  double derivsEx[2] = {0, 0};
  patch.eval_deriv(0, 0, derivs);
  EXPECT_FLOAT_EQ(derivs[0], derivsEx[0]);
  EXPECT_FLOAT_EQ(derivs[1], derivsEx[1]);

  patch.eval_deriv(1, 0, derivs);
  EXPECT_FLOAT_EQ(derivs[0], derivsEx[0]);
  EXPECT_FLOAT_EQ(derivs[1], derivsEx[1]);

  patch.eval_deriv(0, 1, derivs);
  EXPECT_FLOAT_EQ(derivs[0], derivsEx[0]);
  EXPECT_FLOAT_EQ(derivs[1], derivsEx[1]);
}

TEST(TaylorPatch, OffsetPlane)
{
  // test plane through z=1
  std::vector<utils::Point> pts = {utils::Point(0, 0, 1), utils::Point(1, 0, 1), utils::Point(0, 1, 1)};
  utils::impl::TaylorPatch patch;
  patch.construct_patch(pts);

  // for this case, the (u, v, w) coordinate system is aligned with (x,y,z)
  expect_float_eq(patch.eval_point(0, 0), utils::Point(0, 0, 1));
  expect_float_eq(patch.eval_point(1, 0), utils::Point(1, 0, 1));
  expect_float_eq(patch.eval_point(0, 1), utils::Point(0, 1, 1));

  double derivs[2];
  double derivsEx[2] = {0, 0};
  patch.eval_deriv(0, 0, derivs);
  EXPECT_FLOAT_EQ(derivs[0], derivsEx[0]);
  EXPECT_FLOAT_EQ(derivs[1], derivsEx[1]);

  patch.eval_deriv(1, 0, derivs);
  EXPECT_FLOAT_EQ(derivs[0], derivsEx[0]);
  EXPECT_FLOAT_EQ(derivs[1], derivsEx[1]);

  patch.eval_deriv(0, 1, derivs);
  EXPECT_FLOAT_EQ(derivs[0], derivsEx[0]);
  EXPECT_FLOAT_EQ(derivs[1], derivsEx[1]);
}

TEST(TaylorPatch, GeneralPlane)
{
  // test plane through x + y + z = 1
  std::vector<utils::Point> pts = {utils::Point(0, 0, 1), utils::Point(1, 0, 0), utils::Point(0, 1, 0)};
  utils::impl::TaylorPatch patch;
  patch.construct_patch(pts);

  // for this case, the (u, v, w) coordinate system is aligned with (x,y,z)
  double tol = 1e-15;
  expect_float_eq(patch.eval_point(0, 0), utils::Point(0, 0, 1), tol);
  expect_float_eq(patch.eval_point(1, 0), utils::Point(1, 0, 0), tol);
  expect_float_eq(patch.eval_point(0, 1), utils::Point(0, 1, 0), tol);
  expect_float_eq(patch.eval_point(1, 1), utils::Point(1, 1, -1), tol);

  double derivs[2];
  double derivsEx[2] = {-1, -1};
  patch.eval_deriv(0, 0, derivs);
  EXPECT_FLOAT_EQ(derivs[0], derivsEx[0]);
  EXPECT_FLOAT_EQ(derivs[1], derivsEx[1]);

  patch.eval_deriv(1, 0, derivs);
  EXPECT_FLOAT_EQ(derivs[0], derivsEx[0]);
  EXPECT_FLOAT_EQ(derivs[1], derivsEx[1]);

  patch.eval_deriv(0, 1, derivs);
  EXPECT_FLOAT_EQ(derivs[0], derivsEx[0]);
  EXPECT_FLOAT_EQ(derivs[1], derivsEx[1]);

  patch.eval_deriv(1, 1, derivs);
  EXPECT_FLOAT_EQ(derivs[0], derivsEx[0]);
  EXPECT_FLOAT_EQ(derivs[1], derivsEx[1]);
}

TEST(TaylorPatch, GeneralPlaneCenteredElsewhere)
{
  // test plane through x + y + z = 1
  std::vector<utils::Point> pts = {utils::Point(0, 0, 1), utils::Point(1, 0, 0), utils::Point(0, 1, 0)};
  utils::impl::TaylorPatch patch;
  patch.construct_patch(pts, utils::Point(2, 2));

  // for this case, the (u, v, w) coordinate system is aligned with (x,y,z)
  double tol = 1e-14;
  expect_float_eq(patch.eval_point(0, 0), utils::Point(0, 0, 1), tol);
  expect_float_eq(patch.eval_point(1, 0), utils::Point(1, 0, 0), tol);
  expect_float_eq(patch.eval_point(0, 1), utils::Point(0, 1, 0), tol);
  expect_float_eq(patch.eval_point(1, 1), utils::Point(1, 1, -1), tol);

  double derivs[2];
  double derivsEx[2] = {-1, -1};
  patch.eval_deriv(0, 0, derivs);
  EXPECT_FLOAT_EQ(derivs[0], derivsEx[0]);
  EXPECT_FLOAT_EQ(derivs[1], derivsEx[1]);

  patch.eval_deriv(1, 0, derivs);
  EXPECT_FLOAT_EQ(derivs[0], derivsEx[0]);
  EXPECT_FLOAT_EQ(derivs[1], derivsEx[1]);

  patch.eval_deriv(0, 1, derivs);
  EXPECT_FLOAT_EQ(derivs[0], derivsEx[0]);
  EXPECT_FLOAT_EQ(derivs[1], derivsEx[1]);

  patch.eval_deriv(1, 1, derivs);
  EXPECT_FLOAT_EQ(derivs[0], derivsEx[0]);
  EXPECT_FLOAT_EQ(derivs[1], derivsEx[1]);
}

TEST(TaylorPatch, GeneralPlaneShifted)
{
  // test plane through (x - 1) + (y - 1) + (z - 1) = 1
  std::vector<utils::Point> pts = {utils::Point(1, 1, 2), utils::Point(2, 1, 1), utils::Point(1, 2, 1)};
  utils::impl::TaylorPatch patch;
  patch.construct_patch(pts);

  double tol = 1e-15;
  expect_float_eq(patch.eval_point(1, 1), utils::Point(1, 1, 2), tol);
  expect_float_eq(patch.eval_point(2, 1), utils::Point(2, 1, 1), tol);
  expect_float_eq(patch.eval_point(1, 2), utils::Point(1, 2, 1), tol);
  expect_float_eq(patch.eval_point(2, 2), utils::Point(2, 2, 0), tol);

  double derivs[2];
  double derivsEx[2] = {-1, -1};
  patch.eval_deriv(1, 1, derivs);
  EXPECT_FLOAT_EQ(derivs[0], derivsEx[0]);
  EXPECT_FLOAT_EQ(derivs[1], derivsEx[1]);

  patch.eval_deriv(2, 1, derivs);
  EXPECT_FLOAT_EQ(derivs[0], derivsEx[0]);
  EXPECT_FLOAT_EQ(derivs[1], derivsEx[1]);

  patch.eval_deriv(1, 2, derivs);
  EXPECT_FLOAT_EQ(derivs[0], derivsEx[0]);
  EXPECT_FLOAT_EQ(derivs[1], derivsEx[1]);

  patch.eval_deriv(2, 2, derivs);
  EXPECT_FLOAT_EQ(derivs[0], derivsEx[0]);
  EXPECT_FLOAT_EQ(derivs[1], derivsEx[1]);
}

TEST(TaylorPatch, NonOrthogonalPlane)
{
  // test plane through (x - 1) + (y - 1) + (z - 1) = 1
  std::vector<utils::Point> pts = {utils::Point(1, 1, 2), utils::Point(2, 1, 1), utils::Point(2, 2, 0)};
  utils::impl::TaylorPatch patch;
  patch.construct_patch(pts);

  // for this case, the (u, v, w) coordinate system is aligned with (x,y,z)
  double tol = 1e-14;
  expect_float_eq(patch.eval_point(1, 1), utils::Point(1, 1, 2), tol);
  expect_float_eq(patch.eval_point(2, 1), utils::Point(2, 1, 1), tol);
  expect_float_eq(patch.eval_point(1, 2), utils::Point(1, 2, 1), tol);
  expect_float_eq(patch.eval_point(2, 2), utils::Point(2, 2, 0), tol);

  double derivs[2];
  double derivsEx[2] = {-1, -1};
  patch.eval_deriv(1, 1, derivs);
  EXPECT_FLOAT_EQ(derivs[0], derivsEx[0]);
  EXPECT_FLOAT_EQ(derivs[1], derivsEx[1]);

  patch.eval_deriv(2, 1, derivs);
  EXPECT_FLOAT_EQ(derivs[0], derivsEx[0]);
  EXPECT_FLOAT_EQ(derivs[1], derivsEx[1]);

  patch.eval_deriv(1, 2, derivs);
  EXPECT_FLOAT_EQ(derivs[0], derivsEx[0]);
  EXPECT_FLOAT_EQ(derivs[1], derivsEx[1]);

  patch.eval_deriv(2, 2, derivs);
  EXPECT_FLOAT_EQ(derivs[0], derivsEx[0]);
  EXPECT_FLOAT_EQ(derivs[1], derivsEx[1]);
}

template <typename Tx, typename Ty>
void test_deriv(const double x, const double y, utils::impl::TaylorPatch& patch, Tx fDx, Ty fDy)
{
  double derivs[2];
  double derivsEx[2] = {fDx(x, y), fDy(x, y)};
  patch.eval_deriv(x, y, derivs);
  EXPECT_FLOAT_EQ(derivs[0], derivsEx[0]);
  EXPECT_FLOAT_EQ(derivs[1], derivsEx[1]);
}

TEST(TaylorPatch, Quadratic)
{
  // use quadratic ax + by + cz + dx*x + e*x*y + f*y*y = g
  const double a = 1, b = 2, c = 3, d = 4, e = 5, f = 6, g = 7;

  auto fplus = [&](const double x, const double y) {
    return (g - a * x - b * y - d * x * x - e * x * y - f * y * y) / c;
  };

  auto fplusDx = [&](const double x, const double y) { return (-a - 2 * d * x - e * y) / c; };

  auto fplusDy = [&](const double x, const double y) { return (-b - e * x - 2 * f * y) / c; };

  std::vector<utils::Point> pts = {utils::Point(1, 1, fplus(1, 1)), utils::Point(2, 1, fplus(2, 1)),
                                   utils::Point(2, 2, fplus(2, 2)), utils::Point(0, 1, fplus(0, 1)),
                                   utils::Point(1, 0, fplus(1, 0)), utils::Point(0, 0, fplus(0, 0))};

  utils::impl::TaylorPatch patch;
  std::cout << "constructing patch" << std::endl;
  patch.construct_patch(pts);

  std::cout << "\nevaluating patch" << std::endl;
  double tol = 1e-13;
  expect_float_eq(patch.eval_point(1, 1), utils::Point(1, 1, fplus(1, 1)), tol);
  expect_float_eq(patch.eval_point(2, 1), utils::Point(2, 1, fplus(2, 1)), tol);
  expect_float_eq(patch.eval_point(2, 2), utils::Point(2, 2, fplus(2, 2)), tol);
  expect_float_eq(patch.eval_point(0, 1), utils::Point(0, 1, fplus(0, 1)), tol);
  expect_float_eq(patch.eval_point(1, 0), utils::Point(1, 0, fplus(1, 0)), tol);
  expect_float_eq(patch.eval_point(0, 0), utils::Point(0, 0, fplus(0, 0)), tol);
  expect_float_eq(patch.eval_point(1.5, 1.5), utils::Point(1.5, 1.5, fplus(1.5, 1.5)), tol);
  expect_float_eq(patch.eval_point(1.5, 1.3), utils::Point(1.5, 1.3, fplus(1.5, 1.3)), tol);
  expect_float_eq(patch.eval_point(-1, -1), utils::Point(-1, -1, fplus(-1, -1)), tol);

  test_deriv(1, 1, patch, fplusDx, fplusDy);
  test_deriv(2, 1, patch, fplusDx, fplusDy);
  test_deriv(2, 2, patch, fplusDx, fplusDy);
  test_deriv(0, 1, patch, fplusDx, fplusDy);
  test_deriv(1, 0, patch, fplusDx, fplusDy);
  test_deriv(0, 0, patch, fplusDx, fplusDy);
  test_deriv(-1, -1, patch, fplusDx, fplusDy);
}

TEST(TaylorPatch, SinglePoint)
{
  std::vector<utils::Point> pts{utils::Point(1, 2, 3)};

  utils::impl::TaylorPatch patch;
  patch.construct_patch(pts);

  double tol = 1e-14;
  expect_float_eq(patch.eval_point(1, 2), utils::Point(1, 2, 3), tol);
  expect_float_eq(patch.eval_point(2, 3), utils::Point(2, 3, 3), tol);
  {
    double derivs[2];
    double derivsEx[2] = {0, 0};

    patch.eval_deriv(2, 1, derivs);
    EXPECT_FLOAT_EQ(derivs[0], derivsEx[0]);
    EXPECT_FLOAT_EQ(derivs[1], derivsEx[1]);
  }

  {
    double derivs[2];
    double derivsEx[2] = {0, 0};
    patch.eval_deriv(1, 2, derivs);
    EXPECT_FLOAT_EQ(derivs[0], derivsEx[0]);
    EXPECT_FLOAT_EQ(derivs[1], derivsEx[1]);
  }
}

TEST(TaylorPatch, Perturbation)
{
  // the z coordinates of some points are perturbed by machine epsilon
  // If the method used to compute the taylor patch is not stable, this small
  // perturbation can cause large variations in the computed surface
  utils::Point pt0(0.8000000000000002, 0.7999999999999998, 0.9999999999999996);
  utils::Point pt1(1, 0.7999999999999999, 0.9999999999999997);
  utils::Point pt2(0.8000000000000002, 0.6, 0.9999999999999996);
  utils::Point pt3(1, 0.6, 0.9999999999999994);
  utils::Point pt4(0.6000000000000001, 0.6, 0.9999999999999996);
  utils::Point pt5(0.6000000000000001, 0.7999999999999998, 0.9999999999999996);
  utils::Point pt6(1, 1, 0.9999999999999997);
  utils::Point pt7(0.8, 1, 0.9999999999999996);
  utils::Point pt8(0.6000000000000001, 1, 0.9999999999999996);

  std::vector<utils::Point> pts = {pt0, pt1, pt2, pt3, pt4, pt5, pt6, pt7, pt8};

  utils::impl::TaylorPatch patch;
  patch.construct_patch(pts);

  double tol = 1e-14;
  EXPECT_NEAR(patch.eval_point(0.8,  0.8).z,  1.0, tol);
  EXPECT_NEAR(patch.eval_point(0.85, 0.8).z,  1.0, tol);
  EXPECT_NEAR(patch.eval_point(0.8,  0.85).z, 1.0, tol);
  EXPECT_NEAR(patch.eval_point(0.85, 0.85).z, 1.0, tol);
  EXPECT_NEAR(patch.eval_point(0.75, 0.8).z,  1.0, tol);
  EXPECT_NEAR(patch.eval_point(0.8,  0.75).z, 1.0, tol);
  EXPECT_NEAR(patch.eval_point(0.75, 0.75).z, 1.0, tol);
}

// TODO: test single point case

} // namespace impl
} // namespace middle_mesh
} // namespace stk
