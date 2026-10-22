// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include <gtest/gtest.h>
#include <Akri_ContourUtils.hpp>
#include <Akri_UnitTestUtils.hpp>

#include <array>

namespace krino {

// compute_linear_edge_crossing and compute_quadratic_edge_crossing are the pure
// interpolation helpers the ContourTri/ContourTet decompositions use to place a
// vertex on a signed edge. They take vertex coordinates plus per-node distances
// and return the point where the distance crosses zero -- no mesh required.

//--------------------------------------------------------------------------
// Linear crossing: loc = d0/(d0-d1); the returned point is
// (1-loc)*coords[i0] + loc*coords[i1].
//--------------------------------------------------------------------------

// Symmetric distances put the crossing at the edge midpoint.
TEST(ContourUtils, linearCrossing_symmetricDistances_isMidpoint)
{
  const std::array<stk::math::Vector3d,3> coords = {{ {0.0, 0.0, 0.0}, {1.0, 0.0, 0.0}, {0.0, 1.0, 0.0} }};
  const std::array<double,3> dist = {{ 2.0, -2.0, 5.0 }};

  const stk::math::Vector3d x = compute_linear_edge_crossing(coords, dist, 0u, 1u);
  expect_eq(stk::math::Vector3d(0.5, 0.0, 0.0), x);
}

// Asymmetric distances weight the crossing toward the smaller-magnitude node.
// d0=1, d1=-3 -> loc = 1/(1-(-3)) = 0.25.
TEST(ContourUtils, linearCrossing_asymmetricDistances_weightsTowardSmallerMagnitude)
{
  const std::array<stk::math::Vector3d,3> coords = {{ {0.0, 0.0, 0.0}, {1.0, 0.0, 0.0}, {0.0, 1.0, 0.0} }};
  const std::array<double,3> dist = {{ 1.0, -3.0, 5.0 }};

  const stk::math::Vector3d x = compute_linear_edge_crossing(coords, dist, 0u, 1u);
  expect_eq(stk::math::Vector3d(0.25, 0.0, 0.0), x);
}

// The crossing is orientation-symmetric: swapping the two endpoints returns the
// same physical point.
TEST(ContourUtils, linearCrossing_isOrientationSymmetric)
{
  const std::array<stk::math::Vector3d,3> coords = {{ {0.0, 0.0, 0.0}, {1.0, 0.0, 0.0}, {0.0, 1.0, 0.0} }};
  const std::array<double,3> dist = {{ 1.0, -3.0, 5.0 }};

  const stk::math::Vector3d forward = compute_linear_edge_crossing(coords, dist, 0u, 1u);
  const stk::math::Vector3d reverse = compute_linear_edge_crossing(coords, dist, 1u, 0u);
  expect_eq(forward, reverse);
}

// A non-axis edge (nodes 1 and 2) interpolates in both components.
TEST(ContourUtils, linearCrossing_diagonalEdge_interpolatesBothComponents)
{
  const std::array<stk::math::Vector3d,3> coords = {{ {0.0, 0.0, 0.0}, {1.0, 0.0, 0.0}, {0.0, 1.0, 0.0} }};
  const std::array<double,3> dist = {{ 5.0, 1.0, -1.0 }};

  // loc = 1/(1-(-1)) = 0.5 -> midpoint of the (1,0)-(0,1) edge.
  const stk::math::Vector3d x = compute_linear_edge_crossing(coords, dist, 1u, 2u);
  expect_eq(stk::math::Vector3d(0.5, 0.5, 0.0), x);
}

// The 4-coordinate (tet) instantiation behaves the same.
TEST(ContourUtils, linearCrossing_tetInstantiation_interpolatesEdge)
{
  const std::array<stk::math::Vector3d,4> coords = {{ {0.0, 0.0, 0.0}, {0.0, 0.0, 1.0}, {1.0, 0.0, 0.0}, {0.0, 1.0, 0.0} }};
  const std::array<double,4> dist = {{ 1.0, -1.0, 5.0, 5.0 }};

  const stk::math::Vector3d x = compute_linear_edge_crossing(coords, dist, 0u, 1u);
  expect_eq(stk::math::Vector3d(0.0, 0.0, 0.5), x);
}

//--------------------------------------------------------------------------
// Quadratic crossing: uses the vertex distances at i0, i1 and the mid-edge
// distance at i2, then interpolates between coords[i0] and coords[i1].
//--------------------------------------------------------------------------

// A midside distance that lies on the straight line between the endpoints
// reduces to the linear result: d0=2, d1=-2, mid=0 -> midpoint.
TEST(ContourUtils, quadraticCrossing_linearProfile_isMidpoint)
{
  const std::array<stk::math::Vector3d,3> coords = {{ {0.0, 0.0, 0.0}, {1.0, 0.0, 0.0}, {0.0, 1.0, 0.0} }};
  // tri6 distance array: [v0, v1, v2, m01, m12, m20]
  const std::array<double,6> dist = {{ 2.0, -2.0, 5.0, 0.0, 5.0, 5.0 }};

  // Edge 0-1 has mid-edge node ordinal 3.
  const stk::math::Vector3d x = compute_quadratic_edge_crossing(coords, dist, 0u, 1u, 3u);
  expect_eq(stk::math::Vector3d(0.5, 0.0, 0.0), x);
}

// A genuinely curved profile: vertex distances 1, -3 with midside -1 solve to
// loc = 0.25 (the quadratic root along the edge).
TEST(ContourUtils, quadraticCrossing_curvedProfile_solvesQuadraticRoot)
{
  const std::array<stk::math::Vector3d,3> coords = {{ {0.0, 0.0, 0.0}, {1.0, 0.0, 0.0}, {0.0, 1.0, 0.0} }};
  const std::array<double,6> dist = {{ 1.0, -3.0, 5.0, -1.0, 5.0, 5.0 }};

  const stk::math::Vector3d x = compute_quadratic_edge_crossing(coords, dist, 0u, 1u, 3u);
  expect_eq(stk::math::Vector3d(0.25, 0.0, 0.0), x);
}

}
