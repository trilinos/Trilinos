// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include <gtest/gtest.h>
#include <Akri_ContourSubElement.hpp>

#include <cmath>

namespace krino {

// Most of ContourSubElement is bound to a live ContourElement owner (it queries
// the owner for distances, the parametric element library, and the length
// scale), so the end-to-end decomposition is exercised through the
// ContourElement tests. The two pieces that are pure and can be pinned directly
// are the static helpers is_more and find_quadratic_crossing.

//--------------------------------------------------------------------------
// is_more: a tolerance-based lexicographic ordering of points, comparing x
// first, then y, then z. Two equal points compare as "not more".
//--------------------------------------------------------------------------

TEST(ContourSubElement, isMore_equalPoints_returnsFalse)
{
  const stk::math::Vector3d p{1.0, 2.0, 3.0};
  EXPECT_FALSE(ContourSubElement::is_more(p, p));
}

TEST(ContourSubElement, isMore_largerX_dominates)
{
  // x decides regardless of the trailing components.
  EXPECT_TRUE(ContourSubElement::is_more({2.0, 0.0, 0.0}, {1.0, 9.0, 9.0}));
  EXPECT_FALSE(ContourSubElement::is_more({1.0, 9.0, 9.0}, {2.0, 0.0, 0.0}));
}

TEST(ContourSubElement, isMore_tiedX_breaksOnY)
{
  EXPECT_TRUE(ContourSubElement::is_more({1.0, 2.0, 0.0}, {1.0, 1.0, 9.0}));
  EXPECT_FALSE(ContourSubElement::is_more({1.0, 1.0, 9.0}, {1.0, 2.0, 0.0}));
}

TEST(ContourSubElement, isMore_tiedXandY_breaksOnZ)
{
  EXPECT_TRUE(ContourSubElement::is_more({1.0, 2.0, 3.0}, {1.0, 2.0, 2.0}));
  EXPECT_FALSE(ContourSubElement::is_more({1.0, 2.0, 2.0}, {1.0, 2.0, 3.0}));
}

//--------------------------------------------------------------------------
// find_quadratic_crossing: given the distances d0, d1 at the endpoints and d2
// at the mid-edge, return the parametric location in [0,1] where the quadratic
// profile crosses zero. Inputs must straddle zero exactly once.
//--------------------------------------------------------------------------

// Any endpoint/midside distance that is (near) zero is snapped to the
// corresponding parametric location.
TEST(ContourSubElement, findQuadraticCrossing_zeroEndpointsSnapToEnds)
{
  EXPECT_DOUBLE_EQ(0.0, ContourSubElement::find_quadratic_crossing(0.0, -1.0, -0.5));
  EXPECT_DOUBLE_EQ(1.0, ContourSubElement::find_quadratic_crossing(1.0, 0.0, 0.5));
  EXPECT_DOUBLE_EQ(0.5, ContourSubElement::find_quadratic_crossing(1.0, -1.0, 0.0));
}

// When the mid-edge distance lies on the straight line between the endpoints the
// quadratic term vanishes and the root reduces to the linear crossing
// d0/(d0-d1). Here d0=1, d1=-3, mid=-1 (the linear average) -> loc 0.25.
TEST(ContourSubElement, findQuadraticCrossing_linearProfile_reducesToLinearRoot)
{
  EXPECT_DOUBLE_EQ(0.25, ContourSubElement::find_quadratic_crossing(1.0, -3.0, -1.0));
}

// A genuinely curved profile: d0=1, d1=-1, mid=0.5 describes d(t)=1-2t^2, whose
// root in [0,1] is 1/sqrt(2).
TEST(ContourSubElement, findQuadraticCrossing_curvedProfile_solvesQuadratic)
{
  const double loc = ContourSubElement::find_quadratic_crossing(1.0, -1.0, 0.5);
  EXPECT_NEAR(1.0/std::sqrt(2.0), loc, 1.e-12);
}

}
