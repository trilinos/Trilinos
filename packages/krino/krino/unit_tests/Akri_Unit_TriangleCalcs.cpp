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

class TriAtXEqualsOne : public ::testing::Test
{
protected:
    std::array<stk::math::Vector3d,3> triNodeLocations
    {{
        {1,0,-1},
        {1,1,1},
        {1,-1,1}
    }};
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

}





