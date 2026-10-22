// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include <gtest/gtest.h>
#include <Akri_AdaptiveContourUtils.hpp>

#include <array>

namespace krino {

// AdaptiveContourUtils holds the pure predicates and helpers that drive the
// adaptive (recursive) contour refinement. None of them touch a mesh -- they
// work on plain coordinate/distance arrays -- so they can be pinned directly.

//--------------------------------------------------------------------------
// have_possibly_cut_edge: an edge might be cut when either endpoint's distance
// magnitude is within half the edge length (dist^2 <= 0.25*edgeLengthSquared).
//--------------------------------------------------------------------------

// Distances far larger than every edge length: nothing can be cut.
TEST(AdaptiveContourUtils, havePossiblyCutEdge_allDistancesFar_returnsFalse)
{
  const std::array<stk::math::Vector3d,3> coords = {{ {0.0, 0.0, 0.0}, {1.0, 0.0, 0.0}, {0.0, 1.0, 0.0} }};
  const std::array<double,3> dist = {{ 10.0, 10.0, 10.0 }};
  EXPECT_FALSE(have_possibly_cut_edge(coords, dist));
}

// One node close (relative to its incident edge of length 1) flags a possible cut.
TEST(AdaptiveContourUtils, havePossiblyCutEdge_oneNodeClose_returnsTrue)
{
  const std::array<stk::math::Vector3d,3> coords = {{ {0.0, 0.0, 0.0}, {1.0, 0.0, 0.0}, {0.0, 1.0, 0.0} }};
  const std::array<double,3> dist = {{ 0.1, 10.0, 10.0 }};
  EXPECT_TRUE(have_possibly_cut_edge(coords, dist));
}

// A node exactly on the interface (distance 0) is always a possible cut.
TEST(AdaptiveContourUtils, havePossiblyCutEdge_nodeOnInterface_returnsTrue)
{
  const std::array<stk::math::Vector3d,4> coords = {{ {0.0, 0.0, 0.0}, {1.0, 0.0, 0.0}, {0.0, 1.0, 0.0}, {0.0, 0.0, 1.0} }};
  const std::array<double,4> dist = {{ 10.0, 10.0, 10.0, 0.0 }};
  EXPECT_TRUE(have_possibly_cut_edge(coords, dist));
}

// An edge whose endpoints have opposite signs is always a possible cut, even
// when both distance magnitudes exceed the narrow band (which happens for a
// non-normalized / non-unit-gradient field). Without the sign-change guard the
// crossing would be silently dropped.
TEST(AdaptiveContourUtils, havePossiblyCutEdge_signChangeOutsideBand_returnsTrue)
{
  const std::array<stk::math::Vector3d,3> coords = {{ {0.0, 0.0, 0.0}, {1.0, 0.0, 0.0}, {0.0, 1.0, 0.0} }};
  // Edge 0-1 has length 1 -> band threshold is |dist| <= 0.5. Both endpoints
  // sit outside the band (|+/-0.9| > 0.5) yet straddle zero.
  const std::array<double,3> dist = {{ -0.9, 0.9, 10.0 }};
  EXPECT_TRUE(have_possibly_cut_edge(coords, dist));
}

// The subElemVertexIndices overload evaluates the same predicate over a chosen
// subset of the parent nodes. Here a tri6 parent is subset down to its three
// corner vertices.
TEST(AdaptiveContourUtils, havePossiblyCutEdge_subElementIndices_evaluatesSubset)
{
  const std::array<stk::math::Vector3d,6> coords = {{ {0.0, 0.0, 0.0}, {1.0, 0.0, 0.0}, {0.0, 1.0, 0.0},
                                                      {0.5, 0.0, 0.0}, {0.5, 0.5, 0.0}, {0.0, 0.5, 0.0} }};
  const std::array<int,3> corners = {{ 0, 1, 2 }};

  const std::array<double,6> distFar = {{ 10.0, 10.0, 10.0, 10.0, 10.0, 10.0 }};
  EXPECT_FALSE(have_possibly_cut_edge(coords, distFar, corners));

  const std::array<double,6> distClose = {{ 0.1, 10.0, 10.0, 10.0, 10.0, 10.0 }};
  EXPECT_TRUE(have_possibly_cut_edge(coords, distClose, corners));
}

//--------------------------------------------------------------------------
// clip_midedge_distance: clamp the mid-edge distance into a band around the
// linear midpoint, [0.75*d0+0.25*d1, 0.25*d0+0.75*d1] (ordered by d0<d1).
//--------------------------------------------------------------------------
TEST(AdaptiveContourUtils, clipMidedgeDistance_withinBand_returnedUnchanged)
{
  // d0=0, d1=4 -> band is [1, 3]; d2=2 is inside.
  EXPECT_DOUBLE_EQ(2.0, clip_midedge_distance(0.0, 4.0, 2.0));
}

TEST(AdaptiveContourUtils, clipMidedgeDistance_aboveBand_clampedToUpper)
{
  // d0=0, d1=4 -> upper bound 3.
  EXPECT_DOUBLE_EQ(3.0, clip_midedge_distance(0.0, 4.0, 5.0));
}

TEST(AdaptiveContourUtils, clipMidedgeDistance_belowBand_clampedToLower)
{
  // d0=0, d1=4 -> lower bound 1.
  EXPECT_DOUBLE_EQ(1.0, clip_midedge_distance(0.0, 4.0, -1.0));
}

TEST(AdaptiveContourUtils, clipMidedgeDistance_reversedEndpoints_clampsSymmetrically)
{
  // d0=4, d1=0 -> band [1, 3] as well; d2=5 clamps to 3, d2=-1 to 1.
  EXPECT_DOUBLE_EQ(3.0, clip_midedge_distance(4.0, 0.0, 5.0));
  EXPECT_DOUBLE_EQ(1.0, clip_midedge_distance(4.0, 0.0, -1.0));
}

//--------------------------------------------------------------------------
// is_edge_converged: the mid-edge distance is within tol of the linear average.
//--------------------------------------------------------------------------
TEST(AdaptiveContourUtils, isEdgeConverged_midNearLinearAverage_isConverged)
{
  // average of 0 and 4 is 2; d2=2 is exactly linear.
  EXPECT_TRUE(is_edge_converged(0.0, 4.0, 2.0, 0.1));
}

TEST(AdaptiveContourUtils, isEdgeConverged_midFarFromLinearAverage_isNotConverged)
{
  // d2=2.5 deviates 0.5 from the linear average, exceeding tol 0.1.
  EXPECT_FALSE(is_edge_converged(0.0, 4.0, 2.5, 0.1));
}

//--------------------------------------------------------------------------
// compute_node_sign: strict three-way sign, with exact zero mapping to 0.
//--------------------------------------------------------------------------
TEST(AdaptiveContourUtils, computeNodeSign_returnsThreeWaySign)
{
  EXPECT_EQ(0, compute_node_sign(0.0));
  EXPECT_EQ(-1, compute_node_sign(-3.0));
  EXPECT_EQ(1, compute_node_sign(2.5));
}

//--------------------------------------------------------------------------
// snap_distance: values with magnitude below the tol are snapped to zero;
// larger ones are left untouched.
//--------------------------------------------------------------------------
TEST(AdaptiveContourUtils, snapDistance_snapsOnlySubTolValuesToZero)
{
  std::array<double,3> dist = {{ 0.001, -0.001, 0.5 }};
  snap_distance(dist, 0.01);
  EXPECT_DOUBLE_EQ(0.0, dist[0]);
  EXPECT_DOUBLE_EQ(0.0, dist[1]);
  EXPECT_DOUBLE_EQ(0.5, dist[2]);
}

//--------------------------------------------------------------------------
// subarray: gather selected entries in the requested order.
//--------------------------------------------------------------------------
TEST(AdaptiveContourUtils, subarray_gathersSelectedEntriesInOrder)
{
  const std::array<double,4> a = {{ 10.0, 20.0, 30.0, 40.0 }};
  const std::array<int,2> indices = {{ 2, 0 }};
  const std::array<double,2> sub = subarray(a, indices);
  EXPECT_DOUBLE_EQ(30.0, sub[0]);
  EXPECT_DOUBLE_EQ(10.0, sub[1]);
}

}
