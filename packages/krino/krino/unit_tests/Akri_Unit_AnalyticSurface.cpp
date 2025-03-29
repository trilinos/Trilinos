#include <gtest/gtest.h>

#include <stk_math/StkVector.hpp>
#include <Akri_UnitTestUtils.hpp>
#include <Akri_AnalyticSurf.hpp>

namespace krino {

TEST(sphere, intersectWithEdge_locationCorrect)
{
  Sphere sphere{stk::math::Vector3d::ZERO, 0.35};
  const stk::math::Vector3d pt0 = stk::math::Vector3d::ZERO;
  const stk::math::Vector3d pt1(0.5,0.0,0.0);
  const double edgeCrossingTol = 1.e-5;
  const auto & [crossing, loc] = sphere.compute_intersection_with_segment(pt0, pt1, edgeCrossingTol);
  EXPECT_EQ(+1, crossing);
  EXPECT_NEAR(0.7, loc, edgeCrossingTol);
}

TEST(cuboid, unrotatedCuboid_correctSignedDistanceForEveryFaceEdgeAndVertex)
{
  stk::math::Vector3d center = stk::math::Vector3d::ZERO;
  stk::math::Vector3d dimensions(1.,2.,0.5);
  Cuboid cuboid{center, dimensions};

  const double tol = 1.e-5;
  EXPECT_NEAR(-0.25, cuboid.point_signed_distance(center), tol);
  // face distances
  EXPECT_NEAR(0.2, cuboid.point_signed_distance(stk::math::Vector3d(0.7,0.,0.)), tol);
  EXPECT_NEAR(-0.1, cuboid.point_signed_distance(stk::math::Vector3d(0.4,0.,0.)), tol);
  EXPECT_NEAR(0.2, cuboid.point_signed_distance(stk::math::Vector3d(-0.7,0.,0.)), tol);
  EXPECT_NEAR(-0.1, cuboid.point_signed_distance(stk::math::Vector3d(-0.4,0.,0.)), tol);
  EXPECT_NEAR(0.2, cuboid.point_signed_distance(stk::math::Vector3d(0.,1.2,0.)), tol);
  EXPECT_NEAR(-0.1, cuboid.point_signed_distance(stk::math::Vector3d(0.,0.9,0.)), tol);
  EXPECT_NEAR(0.2, cuboid.point_signed_distance(stk::math::Vector3d(0.,0.,0.45)), tol);
  EXPECT_NEAR(-0.1, cuboid.point_signed_distance(stk::math::Vector3d(0.,0.,0.15)), tol);
  EXPECT_NEAR(0.2, cuboid.point_signed_distance(stk::math::Vector3d(0.,0.,-0.45)), tol);
  EXPECT_NEAR(-0.1, cuboid.point_signed_distance(stk::math::Vector3d(0.,0.,-0.15)), tol);
  // edge distances
  EXPECT_NEAR(std::sqrt(2.), cuboid.point_signed_distance(0.5*(cuboid.vertex_location(0)+cuboid.vertex_location(1))+stk::math::Vector3d(0.,-1.,-1.)), tol);
  EXPECT_NEAR(std::sqrt(2.), cuboid.point_signed_distance(0.5*(cuboid.vertex_location(1)+cuboid.vertex_location(2))+stk::math::Vector3d(+1.,0.,-1.)), tol);
  EXPECT_NEAR(std::sqrt(2.), cuboid.point_signed_distance(0.5*(cuboid.vertex_location(2)+cuboid.vertex_location(3))+stk::math::Vector3d(0.,+1.,-1.)), tol);
  EXPECT_NEAR(std::sqrt(2.), cuboid.point_signed_distance(0.5*(cuboid.vertex_location(3)+cuboid.vertex_location(0))+stk::math::Vector3d(-1.,0.,-1.)), tol);
  EXPECT_NEAR(std::sqrt(2.), cuboid.point_signed_distance(0.5*(cuboid.vertex_location(4)+cuboid.vertex_location(5))+stk::math::Vector3d(0.,-1.,+1.)), tol);
  EXPECT_NEAR(std::sqrt(2.), cuboid.point_signed_distance(0.5*(cuboid.vertex_location(5)+cuboid.vertex_location(6))+stk::math::Vector3d(+1.,0.,+1.)), tol);
  EXPECT_NEAR(std::sqrt(2.), cuboid.point_signed_distance(0.5*(cuboid.vertex_location(6)+cuboid.vertex_location(7))+stk::math::Vector3d(0.,+1.,+1.)), tol);
  EXPECT_NEAR(std::sqrt(2.), cuboid.point_signed_distance(0.5*(cuboid.vertex_location(7)+cuboid.vertex_location(4))+stk::math::Vector3d(-1.,0.,+1.)), tol);
  EXPECT_NEAR(std::sqrt(2.), cuboid.point_signed_distance(0.5*(cuboid.vertex_location(0)+cuboid.vertex_location(4))+stk::math::Vector3d(-1.,-1.,0.)), tol);
  EXPECT_NEAR(std::sqrt(2.), cuboid.point_signed_distance(0.5*(cuboid.vertex_location(1)+cuboid.vertex_location(5))+stk::math::Vector3d(+1.,-1.,0.)), tol);
  EXPECT_NEAR(std::sqrt(2.), cuboid.point_signed_distance(0.5*(cuboid.vertex_location(2)+cuboid.vertex_location(6))+stk::math::Vector3d(+1.,+1.,0.)), tol);
  EXPECT_NEAR(std::sqrt(2.), cuboid.point_signed_distance(0.5*(cuboid.vertex_location(3)+cuboid.vertex_location(7))+stk::math::Vector3d(-1.,+1.,0.)), tol);
  //corner distances
  EXPECT_NEAR(std::sqrt(3.), cuboid.point_signed_distance(cuboid.vertex_location(0)+stk::math::Vector3d(-1.,-1.,-1.)), tol);
  EXPECT_NEAR(std::sqrt(3.), cuboid.point_signed_distance(cuboid.vertex_location(1)+stk::math::Vector3d(+1.,-1.,-1.)), tol);
  EXPECT_NEAR(std::sqrt(3.), cuboid.point_signed_distance(cuboid.vertex_location(2)+stk::math::Vector3d(+1.,+1.,-1.)), tol);
  EXPECT_NEAR(std::sqrt(3.), cuboid.point_signed_distance(cuboid.vertex_location(3)+stk::math::Vector3d(-1.,+1.,-1.)), tol);
  EXPECT_NEAR(std::sqrt(3.), cuboid.point_signed_distance(cuboid.vertex_location(4)+stk::math::Vector3d(-1.,-1.,+1.)), tol);
  EXPECT_NEAR(std::sqrt(3.), cuboid.point_signed_distance(cuboid.vertex_location(5)+stk::math::Vector3d(+1.,-1.,+1.)), tol);
  EXPECT_NEAR(std::sqrt(3.), cuboid.point_signed_distance(cuboid.vertex_location(6)+stk::math::Vector3d(+1.,+1.,+1.)), tol);
  EXPECT_NEAR(std::sqrt(3.), cuboid.point_signed_distance(cuboid.vertex_location(7)+stk::math::Vector3d(-1.,+1.,+1.)), tol);
}

TEST(cuboid, rotatedAndTranslatedCuboid_correctSignedDistanceFromCenter)
{
  stk::math::Vector3d center(1.,-2.,3.);
  stk::math::Vector3d dimensions(1.,2.,0.5);
  stk::math::Vector3d rotation(0.5,0.25,0.1);
  Cuboid cuboid{center, dimensions, rotation};

  const double tol = 1.e-5;
  EXPECT_NEAR(-0.25, cuboid.point_signed_distance(center), tol);
}

}
