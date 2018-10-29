/*
 * UnitTestStkPlane.cpp
 *
 *  Created on: Oct 18, 2018
 *      Author: mlstate
 */

#include <stk_math/StkPlane.hpp>
#include <gtest/gtest.h>
#include <type_traits>
#include <cmath>
#include <string>
#include <vector>
#include "UnitTestStkVectorUtils.hpp"

namespace
{
TEST(stk_math_stkPlane, whenCallingDefaultConstructor_planeIsNotValid)
{
    stk::math::Plane3d plane;
    EXPECT_FALSE(plane.is_valid());
}

double tol = 1.e-10;

const stk::math::Vector3d unitNormal(0.,1.,0.);
const stk::math::Vector3d nonUnitNormal(0.,10.,0.);
const stk::math::Vector3d originPoint = stk::math::Vector3d::ZERO;
const stk::math::Vector3d inPlaneNonOriginPoint(10.4,0., -500.);
const double nonOriginPointX{13.3};
const double nonOriginPointY{17.7};
const double nonOriginPointY_minus1{nonOriginPointY-1};
const double nonOriginPointY_plus1{nonOriginPointY+1};
const double nonOriginPointY_minus3{nonOriginPointY-3};
const double nonOriginPointY_plus3{nonOriginPointY+3};
const double nonOriginPointZ{-100.};
const stk::math::Vector3d nonOriginPoint(nonOriginPointX, nonOriginPointY, nonOriginPointZ);
const double pointOffsetInNormalDirection = 17.;
const stk::math::Vector3d planePointOffsetInNormalDir{originPoint+ pointOffsetInNormalDirection*unitNormal};

std::vector<stk::math::Vector3d> create_3_points_in_plane
(
    const stk::math::Vector3d& planeNormal,
    const stk::math::Vector3d& planePoint
)
{
    stk::math::Vector3d vectorInPlaneWithPlaneNormal = planeNormal + stk::math::Vector3d{0.,0., 1.};
    stk::math::Vector3d vector1InPlane = stk::math::Cross(planeNormal, vectorInPlaneWithPlaneNormal);
    vector1InPlane.unitize();
    stk::math::Vector3d vector2InPlane = stk::math::Cross(planeNormal, vector1InPlane);
    vector2InPlane.unitize();
    stk::math::Vector3d planePoint2 = planePoint + vector1InPlane;
    stk::math::Vector3d planePoint3 = planePoint + vector2InPlane;
    return {{planePoint, planePoint2, planePoint3}};
}

void test_plane_normal_is_correct(const stk::math::Vector3d& goldPlaneUnitNormal, stk::math::Plane3d plane)
{
    stk::math::unitTestStkVectorUtils::expect_equal(goldPlaneUnitNormal, plane.normal(), tol);
    stk::math::unitTestStkVectorUtils::expect_equal(goldPlaneUnitNormal, ((const stk::math::Plane3d) (plane)).normal(), tol);
}

void test_plane_created_from_normal_and_point_has_correct_normal
(
    const stk::math::Vector3d &planeNormal,
    const stk::math::Vector3d& planePoint,
    const stk::math::Vector3d &goldPlaneUnitNormal
)
{
    stk::math::Plane3d plane(planeNormal, planePoint);
    ASSERT_TRUE(plane.is_valid());
    test_plane_normal_is_correct(goldPlaneUnitNormal, plane);
}

void test_plane_created_from_3_points_has_correct_normal
(
    const stk::math::Vector3d &point1,
    const stk::math::Vector3d &point2,
    const stk::math::Vector3d &point3,
    const stk::math::Vector3d &goldPlaneUnitNormal
)
{
    stk::math::Plane3d plane(point1, point2, point3);
    ASSERT_TRUE(plane.is_valid());
    test_plane_normal_is_correct(goldPlaneUnitNormal, plane);
}

void test_plane_created_has_correct_normal
(
    const stk::math::Vector3d& planeNormal,
    const stk::math::Vector3d& planePoint,
    const stk::math::Vector3d &goldPlaneUnitNormal
)
{
    test_plane_created_from_normal_and_point_has_correct_normal(planeNormal, planePoint, goldPlaneUnitNormal);
    std::vector<stk::math::Vector3d> planePoints = create_3_points_in_plane(planeNormal, planePoint);
    EXPECT_EQ(3u, planePoints.size());
    test_plane_created_from_3_points_has_correct_normal(planePoints[0], planePoints[1], planePoints[2], goldPlaneUnitNormal);
}

TEST(stk_math_stkPlane, whenConstructingPlaneFromUnitNormalWithOriginPoint_makeSureNormalIsCorrect)
{
    test_plane_created_has_correct_normal(unitNormal, originPoint, unitNormal);
}

TEST(stk_math_stkPlane, whenConstructingPlaneFromUnitNormalWithNonOriginPoint_makeSureNormalIsCorrect)
{
    test_plane_created_has_correct_normal(unitNormal, nonOriginPoint, unitNormal);
}

TEST(stk_math_stkPlane, whenConstructingPlaneFromNonUnitNormalWithOriginPoint_makeSureNormalIsCorrect)
{
    test_plane_created_has_correct_normal(nonUnitNormal, originPoint, unitNormal);
}

TEST(stk_math_stkPlane, whenConstructingPlaneFromNonUnitNormalWithNonOriginPoint_makeSureNormalIsCorrect)
{
    test_plane_created_has_correct_normal(nonUnitNormal, nonOriginPoint, unitNormal);
}

void test_constant_is_correct
(
    const stk::math::Plane3d &plane,
    const double goldPlaneConstant
)
{
    EXPECT_NEAR(goldPlaneConstant, plane.constant(), tol);
    EXPECT_NEAR(goldPlaneConstant, ((const stk::math::Plane3d )plane).constant(), tol);
}

void test_plane_created_from_3_points_has_correct_constant
(
    const stk::math::Vector3d &point1,
    const stk::math::Vector3d &point2,
    const stk::math::Vector3d &point3,
    const double goldPlaneConstant
)
{
    stk::math::Plane3d plane(point1, point2, point3);
    ASSERT_TRUE(plane.is_valid());
    test_constant_is_correct(plane, goldPlaneConstant);
}

void test_plane_created_from_normal_and_point_has_correct_constant
(
    const stk::math::Vector3d& planeNormal,
    const stk::math::Vector3d& planePoint,
    const double goldPlaneConstant
)
{
    stk::math::Plane3d plane(planeNormal, planePoint);
    ASSERT_TRUE(plane.is_valid());
    test_constant_is_correct(plane, goldPlaneConstant);
}

void test_plane_created_has_correct_constant
(
    const stk::math::Vector3d& planeNormal,
    const stk::math::Vector3d& planePoint,
    const double goldPlaneConstant
)
{
    test_plane_created_from_normal_and_point_has_correct_constant(planeNormal, planePoint, goldPlaneConstant);
    std::vector<stk::math::Vector3d> planePoints = create_3_points_in_plane(planeNormal, planePoint);
    EXPECT_EQ(3u, planePoints.size());
    test_plane_created_from_3_points_has_correct_constant(planePoints[0], planePoints[1], planePoints[2], goldPlaneConstant);
}

TEST(stk_math_stkPlane, whenConstructingPlaneFromUnitNormalWithOriginPoint_makeSureConstantIsCorrect)
{
    test_plane_created_has_correct_constant(unitNormal, originPoint, 0.0);
}

TEST(stk_math_stkPlane, whenConstructingPlaneFromNonUnitNormalWithOriginPoint_makeSureConstantIsCorrect)
{
    test_plane_created_has_correct_constant(nonUnitNormal, originPoint, 0.0);
}

TEST(stk_math_stkPlane, whenConstructingPlaneFromUnitNormalWithNonOriginaPointButOriginStillInPlane_makeSureConstantIsCorrect)
{
    test_plane_created_has_correct_constant(unitNormal, inPlaneNonOriginPoint, 0.0);
}

TEST(stk_math_stkPlane, whenConstructingPlaneFromUnitNormalWithNonOriginPointAndOriginNotInPlane_makeSureConstantIsCorrect)
{
    test_plane_created_has_correct_constant(unitNormal, planePointOffsetInNormalDir, pointOffsetInNormalDirection);
}

TEST(stk_math_stkPlane, whenConstructingPlaneFromNonUnitNormalWithNonOriginPointAndOriginNotInPlane_makeSureConstantIsCorrect)
{
    test_plane_created_has_correct_constant(nonUnitNormal, planePointOffsetInNormalDir, pointOffsetInNormalDirection);
}

TEST(stk_math_stkPlane, whenConstructingPlaneWith3IdenticalPoints_vectorMustBeInvalid)
{
    stk::math::Vector3d point(17., 18., -19.);
    EXPECT_FALSE(stk::math::Plane3d(point,point,point).is_valid());
}

TEST(stk_math_stkPlane, whenConstructingPlaneWith2IdenticalPointsAndOneOther_vectorMustBeInvalid)
{
    stk::math::Vector3d point1(17., 18., -19.);
    stk::math::Vector3d point2(17., 18., -30.);
    EXPECT_FALSE(stk::math::Plane3d(point1,point1,point2).is_valid());
}

TEST(stk_math_stkPlane, whenConstructingPlaneWith3ColinearPoints_vectorMustBeInvalid)
{
    stk::math::Vector3d point1(10., 0., 0.);
    stk::math::Vector3d point2(10., -17., 0.);
    stk::math::Vector3d point3(10., 34., 0.);
    EXPECT_FALSE(stk::math::Plane3d(point1,point2,point3).is_valid());
}

TEST(stk_math_stkPlane, whenConstructingPlaneWith3DistinctNonColinearPoints_vectorMustBeValid)
{
    stk::math::Vector3d point1(10., 0., 0.);
    stk::math::Vector3d point2(10., 1., 0.);
    stk::math::Vector3d point3(10., 1., 1.);
    EXPECT_TRUE(stk::math::Plane3d(point1,point2,point3).is_valid());
}

void test_plane_segment_intersection
(
    bool goldDoesIntersect,
    double goldIntersectionLocation,
    const stk::math::Vector3d &planeNormal,
    const stk::math::Vector3d &planePoint,
    const stk::math::Vector3d &segmentPoint1,
    const stk::math::Vector3d &segmentPoint2
)
{
    stk::math::Plane3d plane(planeNormal, planePoint);
    double location;
    EXPECT_EQ(goldDoesIntersect,
              plane.intersects_segment(std::array<stk::math::Vector3d,2>{segmentPoint1,segmentPoint2}, location));
    if ( std::isnan(goldIntersectionLocation))
        EXPECT_TRUE(std::isnan(location));
    else
        EXPECT_EQ(goldIntersectionLocation, location);
}

const double nanDoubleValue{std::nan("")};

TEST(stk_math_stkPlane, whenIntersectingSegmentParallelToPlane_noIntersectionFound)
{
    test_plane_segment_intersection(false, nanDoubleValue, nonUnitNormal, nonOriginPoint,
                                    stk::math::Vector3d{1,nonOriginPointY_minus1,0},
                                    stk::math::Vector3d{3,nonOriginPointY_minus1,0} );
}

TEST(stk_math_stkPlane, whenIntersectingSegmentParallelToAndInPlane_noIntersectionFound)
{
    test_plane_segment_intersection(false, nanDoubleValue, nonUnitNormal, nonOriginPoint,
                                    stk::math::Vector3d{1,nonOriginPointY,0},
                                    stk::math::Vector3d{3,nonOriginPointY,0} );
}

TEST(stk_math_stkPlane, whenIntersectingNonParallelSegmentWithBothPointsAbovePlane_noIntersectionFound)
{
    test_plane_segment_intersection(false, nanDoubleValue, nonUnitNormal, nonOriginPoint,
                                    stk::math::Vector3d{1,nonOriginPointY_plus1,0},
                                    stk::math::Vector3d{3,nonOriginPointY_plus3,0} );
}

TEST(stk_math_stkPlane, whenIntersectingNonParallelSegmentWithBothPointsBelowPlane_noIntersectionFound)
{
    test_plane_segment_intersection(false, nanDoubleValue, nonUnitNormal, nonOriginPoint,
                                    stk::math::Vector3d{1,nonOriginPointY_minus1,0},
                                    stk::math::Vector3d{3,nonOriginPointY_minus3,0} );
}

TEST(stk_math_stkPlane, whenIntersectingNonParallelSegmentWithPointsMirroredByPlane_intersectionFound)
{
    test_plane_segment_intersection(true, 0.5, nonUnitNormal, nonOriginPoint,
                                    stk::math::Vector3d{1,nonOriginPointY_plus3,0},
                                    stk::math::Vector3d{3,nonOriginPointY_minus3,0} );
}

TEST(stk_math_stkPlane, whenIntersectingNonParallelSegmentWithPointsOnOppositeSidesOfPlane_intersectionFound)
{
    test_plane_segment_intersection(true, 0.25, nonUnitNormal, nonOriginPoint,
                                    stk::math::Vector3d{1,nonOriginPointY_plus3,0},
                                    stk::math::Vector3d{3,nonOriginPointY_minus1,0} );
}

TEST(stk_math_stkPlane, whenIntersectingSegmentWithNANPlane_noIntersectionFound)
{
    test_plane_segment_intersection(false, nanDoubleValue, {0,0,0}, nonOriginPoint,
                                    stk::math::Vector3d{1,1,1},
                                    stk::math::Vector3d{2,2,2} );
}

void test_signed_distance_from_plane(const double goldDistanceFromPlane)
{
    const stk::math::Vector3d point1 {1., 1., 1.};
    const stk::math::Vector3d point2 {2., -17., 21.};
    const stk::math::Vector3d point3 {4., 17., 34.};
    const stk::math::Plane3d plane(point1, point2, point3);
    const stk::math::Vector3d pointToTest1 = point1 + goldDistanceFromPlane * plane.normal();
    const stk::math::Vector3d pointToTest2 = point2 + goldDistanceFromPlane * plane.normal();
    const stk::math::Vector3d pointToTest3 = point3 + goldDistanceFromPlane * plane.normal();
    EXPECT_NEAR(goldDistanceFromPlane, plane.signed_distance(pointToTest1), tol);
    EXPECT_NEAR(goldDistanceFromPlane, plane.signed_distance(pointToTest2), tol);
    EXPECT_NEAR(goldDistanceFromPlane, plane.signed_distance(pointToTest3), tol);
}

TEST(stk_math_stkPlane, whenComputingSignedDistanceToPointInPlane_return0)
{
    test_signed_distance_from_plane(0.0);
}

TEST(stk_math_stkPlane, whenComputingSignedDistanceToPointOnPositiveSideOfPlane_returnPositiveDistance)
{
    test_signed_distance_from_plane(17.);
}

TEST(stk_math_stkPlane, whenComputingSignedDistanceToPointOnNegativeSideOfPlane_returnNegativeDistance)
{
    test_signed_distance_from_plane(-3.14);
}

TEST(stk_math_stkPlane, whenComputingSignedDistanceToPointFromInvalidPlane_returnNAN)
{
    const stk::math::Vector3d point{1., 1., 1.};
    const stk::math::Plane3d plane(point, point, point);
    EXPECT_FALSE(plane.is_valid());
    EXPECT_TRUE(std::isnan(plane.signed_distance({2., 2., 2.})));
}
}
