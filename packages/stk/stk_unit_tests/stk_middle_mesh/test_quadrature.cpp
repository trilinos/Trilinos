#include "stk_middle_mesh/quadrature.hpp"
#include "gtest/gtest.h"

namespace stk {
namespace middle_mesh {
namespace impl {

namespace {

double integrate(disc1d::impl::Quadrature& quad, disc1d::impl::LegendrePolynomials poly, int degree)
{
  double val = 0;
  for (int i = 0; i < quad.get_num_points(); ++i)
    val += quad.get_weights()[i] * poly.compute(degree, quad.get_points()[i]);

  return val;
}

double get_exact_integral(disc1d::impl::Quadrature& quad, disc1d::impl::LegendrePolynomials poly, int degree)
{
  double valR = poly.compute_antideriv(degree, quad.get_xi_max());
  double valL = poly.compute_antideriv(degree, quad.get_xi_min());

  return valR - valL;
}
} // namespace

TEST(Quadrature, PolynomialExactness)
{
  disc1d::impl::LegendrePolynomials poly;
  std::vector<int> quadDegrees = {1, 3, 5};
  for (int quadDegree : quadDegrees)
  {
    disc1d::impl::Quadrature quad(quadDegree);
    for (int polyDegree = 0; polyDegree <= quadDegree; ++polyDegree)
    {
      double val   = integrate(quad, poly, polyDegree);
      double valEx = get_exact_integral(quad, poly, polyDegree);
      EXPECT_NEAR(val, valEx, 1e-13);
    }
  }
}
} // namespace impl
} // namespace middle_mesh
} // namespace stk
