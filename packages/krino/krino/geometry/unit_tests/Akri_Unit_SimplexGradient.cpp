// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include <gtest/gtest.h>
#include <Akri_SimplexGradient.hpp>
#include <Akri_UnitTestUtils.hpp>

#include <stk_math/StkVector.hpp>

#include <array>
#include <cmath>
#include <vector>

namespace krino {

// SimplexGradient computes the (constant) gradient of a linear field
// interpolated over a simplex, together with the element area/volume. The
// defining property: for a linear field d(x) = g.x + c the recovered gradient
// is exactly g, independent of the element's shape or placement. The tests use
// linear fields with known g and check both the reference simplex (identity
// Jacobian) and a scaled/translated one (non-trivial Jacobian).

//--------------------------------------------------------------------------
// Triangle (2D)
//--------------------------------------------------------------------------

// Reference triangle (0,0),(1,0),(0,1) with field d = 2x + 3y. Gradient (2,3,0),
// area 1/2.
TEST(SimplexGradient, triangle2d_referenceElement_gradientAndArea)
{
  const std::array<stk::math::Vector3d,3> x = {{ {0,0,0}, {1,0,0}, {0,1,0} }};
  const std::array<double,3> d = {{ 0.0, 2.0, 3.0 }};

  stk::math::Vector3d grad;
  double area;
  calculate_triangle2d_gradient_and_area(x, d, grad, area);

  expect_eq(stk::math::Vector3d(2.0, 3.0, 0.0), grad);
  EXPECT_DOUBLE_EQ(0.5, area);
}

// Same field on a scaled + translated triangle (1,1),(3,1),(1,4): the gradient
// is still exactly (2,3,0), and the area is the true triangle area (0.5*2*3=3).
TEST(SimplexGradient, triangle2d_transformedElement_gradientIsInvariant)
{
  const std::array<stk::math::Vector3d,3> x = {{ {1,1,0}, {3,1,0}, {1,4,0} }};
  const auto field = [](const stk::math::Vector3d & p) { return 2.0*p[0] + 3.0*p[1]; };
  const std::array<double,3> d = {{ field(x[0]), field(x[1]), field(x[2]) }};

  stk::math::Vector3d grad;
  double area;
  calculate_triangle2d_gradient_and_area(x, d, grad, area);

  expect_eq(stk::math::Vector3d(2.0, 3.0, 0.0), grad);
  EXPECT_DOUBLE_EQ(3.0, area);
}

// A constant field has zero gradient.
TEST(SimplexGradient, triangle2d_constantField_zeroGradient)
{
  const std::array<stk::math::Vector3d,3> x = {{ {0,0,0}, {1,0,0}, {0,1,0} }};
  const std::array<double,3> d = {{ 5.0, 5.0, 5.0 }};
  expect_eq(stk::math::Vector3d(0.0, 0.0, 0.0), calculate_triangle2d_gradient(x, d));
}

TEST(SimplexGradient, triangle2d_gradientMagnitude)
{
  const std::array<stk::math::Vector3d,3> x = {{ {0,0,0}, {1,0,0}, {0,1,0} }};
  const std::array<double,3> d = {{ 0.0, 2.0, 3.0 }};
  EXPECT_NEAR(std::sqrt(13.0), calculate_triangle2d_gradient_magnitude(x, d), 1.e-12);
}

//--------------------------------------------------------------------------
// Tetrahedron (3D)
//--------------------------------------------------------------------------

// Reference tet with field d = 2x + 3y + 4z. Gradient (2,3,4), volume 1/6.
TEST(SimplexGradient, tetrahedron_referenceElement_gradientAndVolume)
{
  const std::array<stk::math::Vector3d,4> x = {{ {0,0,0}, {1,0,0}, {0,1,0}, {0,0,1} }};
  const std::array<double,4> d = {{ 0.0, 2.0, 3.0, 4.0 }};

  stk::math::Vector3d grad;
  double vol;
  calculate_tetrahedron_gradient_and_volume(x, d, grad, vol);

  expect_eq(stk::math::Vector3d(2.0, 3.0, 4.0), grad);
  EXPECT_DOUBLE_EQ(1.0/6.0, vol);
}

// Same field on a scaled + translated tet: gradient is still exactly (2,3,4),
// and the volume is 1/6 * (2*3*4) = 4.
TEST(SimplexGradient, tetrahedron_transformedElement_gradientIsInvariant)
{
  const std::array<stk::math::Vector3d,4> x = {{ {1,1,1}, {3,1,1}, {1,4,1}, {1,1,5} }};
  const auto field = [](const stk::math::Vector3d & p) { return 2.0*p[0] + 3.0*p[1] + 4.0*p[2]; };
  const std::array<double,4> d = {{ field(x[0]), field(x[1]), field(x[2]), field(x[3]) }};

  stk::math::Vector3d grad;
  double vol;
  calculate_tetrahedron_gradient_and_volume(x, d, grad, vol);

  expect_eq(stk::math::Vector3d(2.0, 3.0, 4.0), grad);
  EXPECT_DOUBLE_EQ(4.0, vol);
}

TEST(SimplexGradient, tetrahedron_constantField_zeroGradient)
{
  const std::array<stk::math::Vector3d,4> x = {{ {0,0,0}, {1,0,0}, {0,1,0}, {0,0,1} }};
  const std::array<double,4> d = {{ -2.0, -2.0, -2.0, -2.0 }};
  expect_eq(stk::math::Vector3d(0.0, 0.0, 0.0), calculate_tetrahedron_gradient(x, d));
}

TEST(SimplexGradient, tetrahedron_gradientMagnitude)
{
  const std::array<stk::math::Vector3d,4> x = {{ {0,0,0}, {1,0,0}, {0,1,0}, {0,0,1} }};
  const std::array<double,4> d = {{ 0.0, 2.0, 3.0, 4.0 }};
  EXPECT_NEAR(std::sqrt(29.0), calculate_tetrahedron_gradient_magnitude(x, d), 1.e-12);
}

//--------------------------------------------------------------------------
// calculate_simplex_gradient: dispatches on the point count (3 -> triangle,
// 4 -> tet).
//--------------------------------------------------------------------------

TEST(SimplexGradient, simplexGradient_threePoints_dispatchesToTriangle)
{
  const std::vector<stk::math::Vector3d> x = {{ {0,0,0}, {1,0,0}, {0,1,0} }};
  const std::vector<double> d = {{ 0.0, 2.0, 3.0 }};
  expect_eq(stk::math::Vector3d(2.0, 3.0, 0.0), calculate_simplex_gradient(x, d));
}

TEST(SimplexGradient, simplexGradient_fourPoints_dispatchesToTetrahedron)
{
  const std::vector<stk::math::Vector3d> x = {{ {0,0,0}, {1,0,0}, {0,1,0}, {0,0,1} }};
  const std::vector<double> d = {{ 0.0, 2.0, 3.0, 4.0 }};
  expect_eq(stk::math::Vector3d(2.0, 3.0, 4.0), calculate_simplex_gradient(x, d));
}

}
