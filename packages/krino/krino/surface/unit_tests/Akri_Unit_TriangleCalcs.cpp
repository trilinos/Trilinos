/*
 * Akri_Unit_TriangleCalcs.cpp
 *
 *  Created on: Sep 7, 2023
 *      Author: drnoble
 */

#include <gtest/gtest.h>

#include <cmath>
#include <vector>

#include <stk_math/StkVector.hpp>
#include <Akri_UnitTestUtils.hpp>
#include <Akri_FacetedSurfaceCalcs.hpp>
#include <Akri_Facet.hpp>

namespace
{

std::array<stk::math::Vector3d,3> triAtXPlus1
{{
    {1,0,-1},
    {1,1,1},
    {1,-1,1}
}};
std::array<stk::math::Vector3d,3> triAtXMinus1
{{
    {-1,0,-1},
    {-1,1,1},
    {-1,-1,1}
}};

class TriAtXEqualsOne : public ::testing::Test
{
protected:
    std::array<stk::math::Vector3d,3> triNodeLocations{triAtXPlus1};
};

void expect_tri_edge_intersection(const stk::math::Vector3d& edgePt1,
         const stk::math::Vector3d& edgePt2,
         const std::array<stk::math::Vector3d, 3> &triLocs,
         const std::pair<bool,double> & goldIntersection)
{
  krino::Facet3d facet(triLocs[0], triLocs[1], triLocs[2]);
  const auto [haveIntersection, intersectionLoc] = krino::compute_facet_edge_intersection(facet, edgePt1, edgePt2);
  EXPECT_EQ(goldIntersection.first, haveIntersection);
  if (haveIntersection)
  {
    EXPECT_NEAR(goldIntersection.second, intersectionLoc, 1.e-6);
  }
}

TEST_F(TriAtXEqualsOne,WhenChordDoesNotIntersect_noIntersections)
{
    expect_tri_edge_intersection({10,0,0}, {3,0,0}, triNodeLocations, {false, -1.});
}

TEST_F(TriAtXEqualsOne,WhenChordIntersectsCenter_oneIntersections)
{
    expect_tri_edge_intersection({0,0,0}, {3,0,0}, triNodeLocations, {true, 1./3.});
}

TEST_F(TriAtXEqualsOne,WhenChordDoesNotIntersectsCenter_NoIntersections)
{
    expect_tri_edge_intersection({0,5,0}, {3,5,0}, triNodeLocations, {false, -1.});
}

TEST_F(TriAtXEqualsOne,WhenChordDoesIntersectsExactlyAtNode_oneIntersection)
{
    expect_tri_edge_intersection({0,0,-1}, {2,0,-1}, triNodeLocations, {true, 0.5});
}

class TriNonAxisAligned : public ::testing::Test
{
protected:
    std::array<stk::math::Vector3d,3> triNodeLocations
    {{
        {1,0,0},
        {0,1,0},
        {0,0,1}
    }};
};

TEST_F(TriNonAxisAligned,WhenChordIntersectsCenter_oneIntersections)
{
    expect_tri_edge_intersection({0,0,0}, {2,2,2}, triNodeLocations,{true, 1./6.});
}

TEST_F(TriNonAxisAligned,WhenChordIntersectsJustInsideNode_oneIntersections)
{
    expect_tri_edge_intersection({0, .01, .01}, {1, .01, .01}, triNodeLocations,{true, 0.98});
}

TEST_F(TriNonAxisAligned,WhenChordIntersectsJustOutsideNode_noIntersections)
{
    expect_tri_edge_intersection({0, -.01, .01}, {1, -.01, .01}, triNodeLocations,{false, -1.});
}

class FourRightTriangles : public ::testing::Test
{
protected:
    const stk::math::Vector3d center{0.0, 0.0, 0.0};
    const stk::math::Vector3d posX{1.0, 0.0, 0.0};
    const stk::math::Vector3d negX{-1.0, 0.0, 0.0};
    const stk::math::Vector3d posY{0.0, 1.0, 0.0};
    const stk::math::Vector3d negY{0.0, -1.0, 0.0};
    std::vector<std::array<stk::math::Vector3d,3>> trisNodeLocations
    { {{ center, posX, posY }},
      {{ center, posY, negX }},
      {{ center, negX, negY }},
      {{ center, negY, posX }},
    };
};

void expect_have_any_tri_segment_intersections(const std::vector<std::array<stk::math::Vector3d,2>> & segments,
         const std::vector<std::array<stk::math::Vector3d,3>> &trisNodeLocations)
{
    size_t numIntersections = 0;
    for (auto & triLocs : trisNodeLocations)
    {
      krino::Facet3d facet(triLocs[0], triLocs[1], triLocs[2]);
      for (auto && segment : segments)
      {
        const auto [haveIntersection, intersectionLoc] = krino::compute_facet_edge_intersection(facet, segment[0], segment[1]);
        if (haveIntersection)
          ++numIntersections;
      }
    }
    EXPECT_TRUE(numIntersections > 0);
}

TEST_F(FourRightTriangles, SingleChordThatPassesThroughNode_GetAtLeastOneIntersection)
{
    expect_have_any_tri_segment_intersections({ {{{0,0,-1}, {0,0,1}}} }, trisNodeLocations);
}

TEST_F(FourRightTriangles, DoubleChordWithNodeOnTriNode_GetAtLeastOneIntersection)
{
    expect_have_any_tri_segment_intersections({ {{{0,0,-1}, {0,0,0}}}, {{{0,0,0}, {0,0,1}}} }, trisNodeLocations);
}

TEST_F(FourRightTriangles, SingleChordThatPassesThroughEdge_GetAtLeastOneIntersection)
{
    expect_have_any_tri_segment_intersections({ {{{0,0.5,-1}, {0,0.5,1}}} }, trisNodeLocations);

}

TEST_F(FourRightTriangles, DoubleChordWithNodeOnTriEdge_GetAtLeastOneIntersection)
{
    expect_have_any_tri_segment_intersections({ {{{0,0.5,-1}, {0,0.5,0}}}, {{{0,0.5,0}, {0,0.5,1}}} }, trisNodeLocations);
}

TEST(SixTetsInHaloAroundCentralEdge, whenCalculatingIntersectionsForStraightCurvePassingThroughAndPerpendicularToCentralEdgeAlongPlaneOfTwoTris_GetFourIntersections)
{
    const std::array<stk::math::Vector3d,2> segment = {{ {-5.000000, -5.000000, -5.000000},{-5.000000, 5.000000, -5.000000} }};
    std::array<stk::math::Vector3d,8> nodeLocs
    {{
      { -5.1000000000000000000000000, -0.9000000000000000000000000, -5.1000000000000000000000000 },
      { -4.8000000000000007105427358, -1.2000000000000001776356839, -5.4000000000000003552713679 },
      { -4.5000000000000000000000000, -1.5000000000000000000000000, -5.1000000000000005329070518 },
      { -4.8000000000000000000000000, -1.8000000000000000000000000, -4.8000000000000000000000000 },
      { -5.1000000000000005329070518, -1.5000000000000000000000000, -4.5000000000000000000000000 },
      { -5.4000000000000003552713679, -1.2000000000000001776356839, -4.8000000000000007105427358 },
      { -4.8000000000000007105427358, -1.2000000000000001776356839, -4.8000000000000007105427358 },
      { -5.1000000000000005329070518, -1.5000000000000000000000000, -5.1000000000000005329070518 }
   }};
    size_t numIntersections = 0;
    for (size_t tri{0}; tri<6 ; tri++)
    {
        krino::Facet3d facet(nodeLocs[7],nodeLocs[tri],nodeLocs[6]);
        const auto [haveIntersection, intersectionLoc] = krino::compute_facet_edge_intersection(facet, segment[0], segment[1]);
        if (haveIntersection)
          ++numIntersections;
    }
    EXPECT_EQ(4u, numIntersections);
}

void expect_closest_point(const stk::math::Vector3d& queryPt,
    const std::vector<const krino::Facet3d*> & nearestFacets,
    const stk::math::Vector3d& goldClosestPt)
{
  const stk::math::Vector3d closestPt = krino::compute_closest_point(queryPt, nearestFacets);
  krino::expect_near(goldClosestPt, closestPt);
}

void expect_same_closest_point_regardless_of_facet_order(const stk::math::Vector3d& queryPt,
    const krino::Facet3d & facetA,
    const krino::Facet3d & facetB)
{
  const std::vector<const krino::Facet3d*> facetsAB{&facetA, &facetB};
  const std::vector<const krino::Facet3d*> facetsBA{&facetB, &facetA};

  const stk::math::Vector3d closestPtAB = krino::compute_closest_point(queryPt, facetsAB);
  const stk::math::Vector3d closestPtBA = krino::compute_closest_point(queryPt, facetsBA);
  krino::expect_near(closestPtAB, closestPtBA);
}

void expect_same_signed_distance_regardless_of_facet_order(const stk::math::Vector3d& queryPt,
    const krino::Facet3d & facetA,
    const krino::Facet3d & facetB)
{
  const std::vector<const krino::Facet3d*> facetsAB{&facetA, &facetB};
  const std::vector<const krino::Facet3d*> facetsBA{&facetB, &facetA};

  const double signedDistanceAB = krino::compute_point_to_facets_distance_by_average_normal(queryPt, facetsAB);
  const double signedDistanceBA = krino::compute_point_to_facets_distance_by_average_normal(queryPt, facetsBA);
  EXPECT_EQ(signedDistanceAB, signedDistanceBA);
}

TEST(closest_facets, evalClosestPoint_getCorrectAnswerRegardlessOfFacetOrder)
{
  krino::Facet3d facetAtXMinus1(triAtXMinus1[0], triAtXMinus1[1], triAtXMinus1[2]);
  krino::Facet3d facetAtXPlus1(triAtXPlus1[0], triAtXPlus1[1], triAtXPlus1[2]);
  std::vector<const krino::Facet3d*> facets{&facetAtXMinus1, &facetAtXPlus1};

  expect_closest_point(stk::math::Vector3d(-0.1,0,0), facets, stk::math::Vector3d(-1,0,0));
  expect_closest_point(stk::math::Vector3d(0.1,0,0), facets, stk::math::Vector3d(1,0,0));

  const stk::math::Vector3d midPt(0,0,0);
  expect_same_closest_point_regardless_of_facet_order(midPt, facetAtXMinus1, facetAtXPlus1);
}

TEST(closest_facets, evalClosestPointAndSignedDistance_facetsThatWereGivingOppositeSignsForDistanceInParallel_getSameAnswerRegardlessOfFacetOrder)
{
  krino::Facet3d facet0(stk::math::Vector3d(0.588142308404064562,0.161857691595935382,0.488142308404064584),
      stk::math::Vector3d(0.588142308404064562,0.161857691595935382,0.511857691595935416),
      stk::math::Vector3d(0.599999999999999978,0.164656927944478504,0.5));
  krino::Facet3d facet1(stk::math::Vector3d(0.623499950831233773,0.173499950831233707,0.476500049168766204),
      stk::math::Vector3d(0.599999999999999867,0.164656927944478504,0.5),
      stk::math::Vector3d(0.623499950831233773,0.173499950831233707,0.523499950831233796));

  const stk::math::Vector3d queryPt(0.599999999999999978, 0.164656927944478532, 0.5);
  expect_same_closest_point_regardless_of_facet_order(queryPt, facet0, facet1);
  expect_same_signed_distance_regardless_of_facet_order(queryPt, facet0, facet1);
}

void expect_distance_sign(const krino::Facet3d & facet, const stk::math::Vector3d& queryPt, const int goldNegPosOrZero)
{
  const double signedDistance = krino::compute_point_to_facets_distance_by_average_normal(queryPt, std::vector<const krino::Facet3d*>{&facet});
  const bool match = (goldNegPosOrZero == 0) ? (signedDistance==0.) : ((goldNegPosOrZero<0) ? (signedDistance < 0.) : (signedDistance > 0.));
  EXPECT_TRUE(match) << "Expected sign " << goldNegPosOrZero << " and got distance " << signedDistance;
}

TEST(point_distance, evalSignedDistance_queryPointsThatAreJustAboveAndBelowOnFacetTolerance_getZeroIfBelowOnFacetTolerance)
{
  const krino::Facet3d facet(stk::math::Vector3d(1,0,0),stk::math::Vector3d(0,1,0), stk::math::Vector3d(0,0,1));
  const stk::math::Vector3d facetCentroid = facet.centroid();
  const stk::math::Vector3d facetNormal = facet.facet_normal();
  const double magGreaterThanOnFacetTol = 1.e-8;
  const double magLessThanOnFacetTol = 1.e-10;

  expect_distance_sign(facet, facetCentroid-magGreaterThanOnFacetTol*facetNormal, -1);
  expect_distance_sign(facet, facetCentroid+magGreaterThanOnFacetTol*facetNormal, +1);

  expect_distance_sign(facet, facetCentroid-magLessThanOnFacetTol*facetNormal, 0);
  expect_distance_sign(facet, facetCentroid+magLessThanOnFacetTol*facetNormal, 0);
}

//--------------------------------------------------------------------------
// (#1) Zero / indeterminate pseudo-normal branch.
//
// A triangle and its reverse-wound twin share every vertex but carry +z and -z
// normals.  A query point whose closest point is on the shared edge/vertex (so
// the average-normal path is taken) makes the two normals cancel to a
// zero-length pseudo-normal.  Both distance and normal free functions then take
// their documented indeterminate-branch fallbacks: distance is forced positive
// (= the true closest distance) and the returned normal degrades to the closest
// facet's own normal.  compute_pseudo_normal itself returns the zero vector.
//--------------------------------------------------------------------------

namespace {
// A query point off the +x edge of the triangle so that its closest point lands
// on that edge (an edge/vertex projection, not the face interior), forcing the
// average-normal / pseudo-normal path.
const stk::math::Vector3d offEdgeQueryPt{2,0,0};

// A triangle and its reverse-wound twin: identical vertices, +z and -z normals.
// Built as a local struct (not a static-init global) so its Facet3d members and
// the pointer vector referencing them are constructed inside each test.
struct OpposingCoincidentTriangles
{
  const krino::Facet3d plusZ{stk::math::Vector3d(0,0,0), stk::math::Vector3d(1,0,0), stk::math::Vector3d(0,1,0)};
  const krino::Facet3d minusZ{stk::math::Vector3d(0,0,0), stk::math::Vector3d(0,1,0), stk::math::Vector3d(1,0,0)};
  const std::vector<const krino::Facet3d*> facets{&plusZ, &minusZ};
};
}

TEST(zero_pseudo_normal, opposingCoincidentPair_pseudoNormalCancelsToZero)
{
  const OpposingCoincidentTriangles tris;

  const stk::math::Vector3d pseudoNormal = krino::compute_pseudo_normal(offEdgeQueryPt, tris.facets);
  EXPECT_EQ(0., pseudoNormal.length_squared());
}

TEST(zero_pseudo_normal, opposingCoincidentPair_distanceFallsBackToArbitraryPositive)
{
  const OpposingCoincidentTriangles tris;

  // The closest point is the vertex (1,0,0); true distance is 1.  With the
  // pseudo-normal zero, the sign is arbitrarily positive.
  const double signedDistance = krino::compute_point_to_facets_distance_by_average_normal(offEdgeQueryPt, tris.facets);
  EXPECT_NEAR(1., signedDistance, 1.e-9);
}

TEST(zero_pseudo_normal, opposingCoincidentPair_normalFallsBackToClosestFacetNormal)
{
  const OpposingCoincidentTriangles tris;

  // Indeterminate pseudo-normal -> normal degrades to the closest facet's own
  // (unit) normal, whichever facet won the closest-facet tie-break.  Either way
  // the magnitude is 1 and the vector is axis-aligned in z.
  const stk::math::Vector3d normal = krino::compute_closest_point_normal(offEdgeQueryPt, tris.facets);
  EXPECT_NEAR(1., normal.length(), 1.e-9);
  EXPECT_NEAR(1., std::abs(normal[2]), 1.e-9);
}

//--------------------------------------------------------------------------
// (#2) compute_closest_point_normal across its branches.
//
// compute_closest_point_normal returns a *consistent surface normal* (it is not
// query-side relative): interior projections return the closest facet's own
// normal directly; edge/vertex projections return the pseudo-normal-oriented
// unit vector, which for these coplanar configurations is the shared surface
// normal regardless of which side the query sits on.
//--------------------------------------------------------------------------

namespace {
// Two coplanar triangles in z=0 sharing the x-axis edge (0,0,0)->(1,0,0):
// halfA covers +y, halfB covers -y, both with unit normal +z.  A query off the
// shared edge takes the edge/vertex (pseudo-normal) branch with a clean +z
// average normal.  Built as a local struct rather than static-init globals.
struct CoplanarTrianglePair
{
  const krino::Facet3d halfA{stk::math::Vector3d(0,0,0), stk::math::Vector3d(1,0,0), stk::math::Vector3d(0,1,0)};
  const krino::Facet3d halfB{stk::math::Vector3d(0,0,0), stk::math::Vector3d(0,-1,0), stk::math::Vector3d(1,0,0)};
  const std::vector<const krino::Facet3d*> facets{&halfA, &halfB};
};
}

TEST(closest_point_normal, interiorProjection_returnsClosestFacetNormal)
{
  // Query projects onto the triangle interior (all closest-point weights > tol),
  // so the closest facet's own normal (+z) is returned directly.
  const krino::Facet3d facet(stk::math::Vector3d(0,0,0), stk::math::Vector3d(1,0,0), stk::math::Vector3d(0,1,0));
  const std::vector<const krino::Facet3d*> facets{&facet};

  krino::expect_near(stk::math::Vector3d(0,0,1),
      krino::compute_closest_point_normal(stk::math::Vector3d(0.25,0.25,0.5), facets));
  krino::expect_near(stk::math::Vector3d(0,0,1),
      krino::compute_closest_point_normal(stk::math::Vector3d(0.25,0.25,-0.5), facets));
}

TEST(closest_point_normal, edgeProjectionAboveTolerance_returnsQueryRelativeUnitVector)
{
  // Query off the shared x-axis edge at z=1.0 -- far above the on-surface
  // tolerance (sqrTolSoCloseThatItIsOnSurface = 1e-12) -- so the code takes the
  // `x - closestPoint` unitized branch, oriented by the +z pseudo-normal.  For
  // this coplanar surface that direction is the +z surface normal, and because
  // the orientation follows the surface (not the query side) it stays +z whether
  // the query is above or below.
  const CoplanarTrianglePair pair;

  krino::expect_near(stk::math::Vector3d(0,0,1),
      krino::compute_closest_point_normal(stk::math::Vector3d(0.5,0,1.0), pair.facets));
  krino::expect_near(stk::math::Vector3d(0,0,1),
      krino::compute_closest_point_normal(stk::math::Vector3d(0.5,0,-1.0), pair.facets));
}

TEST(closest_point_normal, edgeProjectionWithinOnSurfaceTolerance_returnsUnitPseudoNormal)
{
  // Companion to the test above: query at z=1e-8 gives a squared distance
  // (~1e-16) below sqrTolSoCloseThatItIsOnSurface = 1e-12, so the code takes the
  // *other* branch and returns pseudoNormal.unit_vector() directly rather than
  // the query-relative direction.  Result is still the unit +z surface normal.
  const CoplanarTrianglePair pair;

  const stk::math::Vector3d normal = krino::compute_closest_point_normal(stk::math::Vector3d(0.5,0,1.e-8), pair.facets);
  EXPECT_NEAR(1., normal.length(), 1.e-9);
  krino::expect_near(stk::math::Vector3d(0,0,1), normal);
}

//--------------------------------------------------------------------------
// (#3) compute_pseudo_normal angle-weighting.
//
// In 3D with more than two facets close to the query point, compute_pseudo_normal
// returns the angle-weighted sum of facet normals rather than the plain average.
// With exactly two (or fewer) close facets it returns the plain average.  Both
// paths are exercised, and the sign the pseudo-normal imparts to the average-
// normal signed distance is checked.
//--------------------------------------------------------------------------

TEST(pseudo_normal, twoFacetsAtEdge_usesPlainAverageNormal)
{
  // Two coplanar facets sharing the x-axis edge.  Only two facets are close, so
  // the plain average of the (unit) normals is returned: (0,0,1)+(0,0,1).
  const CoplanarTrianglePair pair;

  const stk::math::Vector3d pseudoNormal = krino::compute_pseudo_normal(stk::math::Vector3d(0.5,0,1.0), pair.facets);
  krino::expect_near(stk::math::Vector3d(0,0,1), pseudoNormal.unit_vector());
}

TEST(pseudo_normal, threeFacetsAtVertex_usesAngleWeightedNormal)
{
  // Three triangles meeting at the origin vertex, all lying in the z=0 plane so
  // every facet normal is +z.  With more than two close facets the angle-weighted
  // branch is taken; since all normals agree the weighted sum still points +z,
  // and (unlike the plain average) its magnitude equals the summed corner angles.
  const krino::Facet3d f0(stk::math::Vector3d(0,0,0), stk::math::Vector3d(1,0,0), stk::math::Vector3d(0,1,0));
  const krino::Facet3d f1(stk::math::Vector3d(0,0,0), stk::math::Vector3d(0,1,0), stk::math::Vector3d(-1,0,0));
  const krino::Facet3d f2(stk::math::Vector3d(0,0,0), stk::math::Vector3d(-1,0,0), stk::math::Vector3d(0,-1,0));
  const std::vector<const krino::Facet3d*> facets{&f0, &f1, &f2};

  // Query just above the shared origin vertex; closest point is the vertex.
  const stk::math::Vector3d pseudoNormal = krino::compute_pseudo_normal(stk::math::Vector3d(0,0,1.e-3), facets);

  ASSERT_GT(pseudoNormal.length(), 0.);
  krino::expect_near(stk::math::Vector3d(0,0,1), pseudoNormal.unit_vector());

  // Each right-isoceles corner at the origin subtends pi/2; three of them sum to
  // 3*pi/2, which is the magnitude of the angle-weighted (+z) pseudo-normal.
  EXPECT_NEAR(3.*M_PI/2., pseudoNormal.length(), 1.e-6);
}

TEST(point_distance, edgeProjectionOntoCoplanarPair_orientsSignedDistanceBySide)
{
  // Edge projection onto the flat two-facet surface: the +z pseudo-normal orients
  // the average-normal signed distance so points above are positive and points
  // below are negative, with magnitude equal to the true distance to the edge.
  const CoplanarTrianglePair pair;

  EXPECT_NEAR( 1., krino::compute_point_to_facets_distance_by_average_normal(stk::math::Vector3d(0.5,0, 1.0), pair.facets), 1.e-9);
  EXPECT_NEAR(-1., krino::compute_point_to_facets_distance_by_average_normal(stk::math::Vector3d(0.5,0,-1.0), pair.facets), 1.e-9);
}

//--------------------------------------------------------------------------
// 2D coverage.
//
// The same free functions are instantiated for Facet2d, but the pseudo-normal
// path differs by dimension: compute_pseudo_normal's angle-weighting is guarded
// by `if constexpr (3 == FACET::DIM)`, so in 2D the plain average_normal is
// always returned (no angle weighting is ever computed).  The 2D segment normal
// follows the crossZ convention: for a segment x0->x1 the normal is
// crossZ(x1-x0).unit_vector(), so (0,0,0)->(1,0,0) has normal (0,-1,0) and its
// reverse has (0,1,0) -- an opposing coincident pair that cancels to zero.
//--------------------------------------------------------------------------

TEST(zero_pseudo_normal_2d, opposingCoincidentSegments_pseudoNormalCancelsToZero)
{
  const krino::Facet2d segForward(stk::math::Vector3d(0,0,0), stk::math::Vector3d(1,0,0));  // normal (0,-1,0)
  const krino::Facet2d segReverse(stk::math::Vector3d(1,0,0), stk::math::Vector3d(0,0,0));  // normal (0, 1,0)
  const std::vector<const krino::Facet2d*> facets{&segForward, &segReverse};

  // Query off the +x end projects onto the shared vertex (1,0,0), forcing the
  // average-normal path; the opposing normals cancel.
  const stk::math::Vector3d pseudoNormal = krino::compute_pseudo_normal(stk::math::Vector3d(2,0,0), facets);
  EXPECT_EQ(0., pseudoNormal.length_squared());
}

TEST(zero_pseudo_normal_2d, opposingCoincidentSegments_distanceFallsBackToArbitraryPositive)
{
  const krino::Facet2d segForward(stk::math::Vector3d(0,0,0), stk::math::Vector3d(1,0,0));
  const krino::Facet2d segReverse(stk::math::Vector3d(1,0,0), stk::math::Vector3d(0,0,0));
  const std::vector<const krino::Facet2d*> facets{&segForward, &segReverse};

  // Closest point is the vertex (1,0,0); true distance is 1, sign arbitrarily +.
  EXPECT_NEAR(1., krino::compute_point_to_facets_distance_by_average_normal(stk::math::Vector3d(2,0,0), facets), 1.e-9);
}

TEST(zero_pseudo_normal_2d, opposingCoincidentSegments_normalFallsBackToClosestFacetNormal)
{
  const krino::Facet2d segForward(stk::math::Vector3d(0,0,0), stk::math::Vector3d(1,0,0));
  const krino::Facet2d segReverse(stk::math::Vector3d(1,0,0), stk::math::Vector3d(0,0,0));
  const std::vector<const krino::Facet2d*> facets{&segForward, &segReverse};

  // Indeterminate pseudo-normal -> degrades to the closest facet's own unit
  // normal, which is ±y aligned whichever segment won the tie-break.
  const stk::math::Vector3d normal = krino::compute_closest_point_normal(stk::math::Vector3d(2,0,0), facets);
  EXPECT_NEAR(1., normal.length(), 1.e-9);
  EXPECT_NEAR(1., std::abs(normal[1]), 1.e-9);
}

TEST(closest_point_normal_2d, interiorProjection_returnsClosestFacetNormal)
{
  // Query projects onto the segment interior (both weights > tol), so the
  // segment's own normal (0,-1,0) is returned directly.
  const krino::Facet2d seg(stk::math::Vector3d(0,0,0), stk::math::Vector3d(2,0,0));
  const std::vector<const krino::Facet2d*> facets{&seg};

  krino::expect_near(stk::math::Vector3d(0,-1,0),
      krino::compute_closest_point_normal(stk::math::Vector3d(1,-0.5,0), facets));
  krino::expect_near(stk::math::Vector3d(0,-1,0),
      krino::compute_closest_point_normal(stk::math::Vector3d(1, 0.5,0), facets));
}

TEST(closest_point_normal_2d, vertexProjection_returnsUnitSurfaceNormalIndependentOfSide)
{
  // Two collinear segments sharing the vertex (1,0,0), both with normal (0,-1,0).
  // A query off that vertex takes the edge/vertex (pseudo-normal) branch and
  // returns the shared (0,-1,0) surface normal regardless of side.
  const krino::Facet2d segA(stk::math::Vector3d(0,0,0), stk::math::Vector3d(1,0,0));
  const krino::Facet2d segB(stk::math::Vector3d(1,0,0), stk::math::Vector3d(2,0,0));
  const std::vector<const krino::Facet2d*> facets{&segA, &segB};

  krino::expect_near(stk::math::Vector3d(0,-1,0),
      krino::compute_closest_point_normal(stk::math::Vector3d(1,-1,0), facets));
  krino::expect_near(stk::math::Vector3d(0,-1,0),
      krino::compute_closest_point_normal(stk::math::Vector3d(1, 1,0), facets));
}

TEST(pseudo_normal_2d, alwaysUsesPlainAverageNormalRegardlessOfCloseCount)
{
  // Two collinear same-normal segments sharing a vertex: 2D never applies angle
  // weighting, so the result is the plain average of the (0,-1,0) normals.
  const krino::Facet2d segA(stk::math::Vector3d(0,0,0), stk::math::Vector3d(1,0,0));
  const krino::Facet2d segB(stk::math::Vector3d(1,0,0), stk::math::Vector3d(2,0,0));
  const std::vector<const krino::Facet2d*> facets{&segA, &segB};

  const stk::math::Vector3d pseudoNormal = krino::compute_pseudo_normal(stk::math::Vector3d(1,-1,0), facets);
  krino::expect_near(stk::math::Vector3d(0,-1,0), pseudoNormal.unit_vector());
}

}





