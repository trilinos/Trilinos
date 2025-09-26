/*
 * Akri_Unit_TriangleCalcs.cpp
 *
 *  Created on: Sep 7, 2023
 *      Author: drnoble
 */

#include <gtest/gtest.h>

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

}





