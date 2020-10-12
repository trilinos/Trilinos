#include "gtest/gtest.h"
#include "stk_math/SideGeometry.hpp"
#include "stk_math/StkVector.hpp"

namespace {

stk::math::QuadGeometry
small_quad(const stk::math::Vector3d & center, double eps)
{
  stk::math::Vector3d epsX(eps, 0.0, 0.0);
  stk::math::Vector3d epsY(0.0, eps, 0.0);

  return stk::math::QuadGeometry(center - epsX - epsY,
                                              center + epsX - epsY,
                                              center + epsX + epsY,
                                              center - epsX + epsY);
}

bool vector_compare_near(const stk::math::Vector3d & expected, const stk::math::Vector3d & computed, double tol)
{
  return (std::abs(computed[0] - expected[0]) < tol) &&
         (std::abs(computed[1] - expected[1]) < tol) &&
         (std::abs(computed[2] - expected[2]) < tol);
}

#define EXPECT_VECTOR_NEAR(v1, v2) EXPECT_TRUE(vector_compare_near(v1, v2, 1.e-12))


TEST(ParticleGeometry, constructWithVectorAndQueryNode)
{
  stk::math::PointGeometry particle({1, 2, 3});

  EXPECT_EQ(stk::math::Vector3d(1, 2, 3), particle.node(0));
}

TEST(ParticleGeometry, centroid)
{
  stk::math::PointGeometry particle({1, 2, 3});

  EXPECT_EQ(stk::math::Vector3d(1, 2, 3), particle.centroid());
}

TEST(ParticleGeometry, projectionToPointOnFace)
{
  stk::math::PointGeometry particle({1.0, 2.0, 3.0});

  EXPECT_VECTOR_NEAR(stk::math::Vector3d(1.0, 2.0, 3.0),
                     particle.closest_proj_on_face(stk::math::Vector3d(1.0, 2.0, 3.0)));
  EXPECT_VECTOR_NEAR(stk::math::Vector3d(1.0, 2.0, 3.0),
                     particle.closest_proj_on_face(stk::math::Vector3d(0.0, 0.0, 0.0)));
  EXPECT_VECTOR_NEAR(stk::math::Vector3d(1.0, 2.0, 3.0),
                     particle.closest_proj_on_face(stk::math::Vector3d(-1.0, -2.0, -3.0)));
}

TEST(ParticleGeometry, distanceToPoint)
{
  stk::math::PointGeometry particle({1.0, 2.0, 3.0});

  EXPECT_DOUBLE_EQ(0.0, particle.min_distance_to_point(stk::math::Vector3d(1.0, 2.0, 3.0)));
  EXPECT_DOUBLE_EQ(1.0, particle.min_distance_to_point(stk::math::Vector3d(0.0, 2.0, 3.0)));
  EXPECT_DOUBLE_EQ(1.0, particle.min_distance_to_point(stk::math::Vector3d(1.0, 3.0, 3.0)));
  EXPECT_DOUBLE_EQ(1.0, particle.min_distance_to_point(stk::math::Vector3d(1.0, 2.0, 2.0)));
  EXPECT_DOUBLE_EQ(5.0, particle.min_distance_to_point(stk::math::Vector3d(4.0, 6.0, 3.0)));
}

TEST(ParticleGeometry, isNodeCloseToSide)
{
  stk::math::Vector3d notClose(2.0, 2.0, 0.0);

  stk::math::Vector3d node(0.0, 0.0, 0.0);

  stk::math::PointGeometry particle(node);

  stk::math::QuadGeometry quadNode = small_quad(node, 0.01);
  stk::math::QuadGeometry quadNotClose = small_quad(notClose, 0.01);

  EXPECT_TRUE(particle.are_nodes_close_to_side(quadNode, 0.01));
  EXPECT_FALSE(particle.are_nodes_close_to_side(quadNotClose, 0.01));
}

TEST(LineGeometry, constructWithVectorsAndQueryNodes)
{
  stk::math::LineGeometry line({1, 2, 3}, {4, 5, 6});

  EXPECT_EQ(stk::math::Vector3d(1, 2, 3), line.node(0));
  EXPECT_EQ(stk::math::Vector3d(4, 5, 6), line.node(1));
}

TEST(LineGeometry, centroid)
{
  stk::math::LineGeometry line({0.0, 0.0, 0.0}, {1.0, 2.0, 3.0});

  EXPECT_EQ(stk::math::Vector3d(0.5, 1.0, 1.5), line.centroid());
}

TEST(LineGeometry, projectionToPointOnLine)
{
  stk::math::LineGeometry line({0.0, 0.0, 0.0}, {1.0, 2.0, 3.0});

  // At the nodes
  EXPECT_VECTOR_NEAR(stk::math::Vector3d(0.0, 0.0, 0.0), line.closest_proj_on_face(stk::math::Vector3d(0.0, 0.0, 0.0)));
  EXPECT_VECTOR_NEAR(stk::math::Vector3d(1.0, 2.0, 3.0), line.closest_proj_on_face(stk::math::Vector3d(1.0, 2.0, 3.0)));

  // In the middle
  EXPECT_VECTOR_NEAR(stk::math::Vector3d(0.1, 0.2, 0.3), line.closest_proj_on_face(stk::math::Vector3d(0.1, 0.2, 0.3)));

  // Offline, projects to middle
  EXPECT_VECTOR_NEAR(stk::math::Vector3d(0.1, 0.2, 0.3), line.closest_proj_on_face(stk::math::Vector3d(0.4,-0.1, 0.4)));

  // Offline, projects to first endpoint
  EXPECT_VECTOR_NEAR(stk::math::Vector3d(0.0, 0.0, 0.0), line.closest_proj_on_face(stk::math::Vector3d(-1.0,-3.0,-2.0)));
  EXPECT_VECTOR_NEAR(stk::math::Vector3d(0.0, 0.0, 0.0), line.closest_proj_on_face(stk::math::Vector3d(-3.0,-2.0,-1.0)));

  // Offline, projects to second endpoint
  EXPECT_VECTOR_NEAR(stk::math::Vector3d(1.0, 2.0, 3.0), line.closest_proj_on_face(stk::math::Vector3d(4.0, 5.0, 6.0)));
  EXPECT_VECTOR_NEAR(stk::math::Vector3d(1.0, 2.0, 3.0), line.closest_proj_on_face(stk::math::Vector3d(6.0, 7.0, 4.0)));
}

TEST(LineGeometry, projectionBetweenTwoLines)
{
  stk::math::LineGeometry line1({0.0, 1.0, 0.0}, {0.0, 0.0, 0.0});
  stk::math::LineGeometry line2({4.0, 0.0, 0.0}, {4.0, 1.0, 0.0});

  EXPECT_VECTOR_NEAR(stk::math::Vector3d(0.0, 0.0, 0.0), line1.closest_proj_on_face(stk::math::Vector3d(4.0, 0.0, 0.0)));
  EXPECT_VECTOR_NEAR(stk::math::Vector3d(0.0, 1.0, 0.0), line1.closest_proj_on_face(stk::math::Vector3d(4.0, 1.0, 0.0)));
  EXPECT_VECTOR_NEAR(stk::math::Vector3d(0.0, 0.5, 0.0), line1.closest_proj_on_face(stk::math::Vector3d(4.0, 0.5, 0.0)));

  EXPECT_VECTOR_NEAR(stk::math::Vector3d(4.0, 1.0, 0.0), line2.closest_proj_on_face(stk::math::Vector3d(0.0, 1.0, 0.0)));
  EXPECT_VECTOR_NEAR(stk::math::Vector3d(4.0, 0.0, 0.0), line2.closest_proj_on_face(stk::math::Vector3d(0.0, 0.0, 0.0)));
  EXPECT_VECTOR_NEAR(stk::math::Vector3d(4.0, 0.5, 0.0), line2.closest_proj_on_face(stk::math::Vector3d(0.0, 0.5, 0.0)));
}

TEST(LineGeometry, distanceBetweenTwoLines)
{
  stk::math::LineGeometry line1({0.0, 1.0, 0.0}, {0.0, 0.0, 0.0});
  stk::math::LineGeometry line2({4.0, 0.0, 0.0}, {4.0, 1.0, 0.0});

  EXPECT_DOUBLE_EQ(4.0, line1.min_distance_to_point(stk::math::Vector3d(4.0, 0.0, 0.0)));
  EXPECT_DOUBLE_EQ(4.0, line1.min_distance_to_point(stk::math::Vector3d(4.0, 1.0, 0.0)));
  EXPECT_DOUBLE_EQ(4.0, line1.min_distance_to_point(stk::math::Vector3d(4.0, 0.5, 0.0)));

  EXPECT_DOUBLE_EQ(4.0, line2.min_distance_to_point(stk::math::Vector3d(0.0, 1.0, 0.0)));
  EXPECT_DOUBLE_EQ(4.0, line2.min_distance_to_point(stk::math::Vector3d(0.0, 0.0, 0.0)));
  EXPECT_DOUBLE_EQ(4.0, line2.min_distance_to_point(stk::math::Vector3d(0.0, 0.5, 0.0)));

  EXPECT_FALSE(line1.are_nodes_close_to_side(line2, 1.1));
  EXPECT_FALSE(line2.are_nodes_close_to_side(line1, 1.1));
}

TEST(Tri3dGeometry, constructWithVectorsAndQueryNodes)
{
  stk::math::TriGeometry tri({1, 2, 3}, {4, 5, 6}, {7, 8, 9});

  EXPECT_EQ(stk::math::Vector3d(1, 2, 3), tri.node(0));
  EXPECT_EQ(stk::math::Vector3d(4, 5, 6), tri.node(1));
  EXPECT_EQ(stk::math::Vector3d(7, 8, 9), tri.node(2));
}

TEST(Tri3dGeometry, centroid)
{
  stk::math::TriGeometry tri({0.0, 0.0, 0.0}, {1.0, 2.0, 3.0}, {0.0, 2.0, 0.0});

  EXPECT_EQ(stk::math::Vector3d(1.0/3.0, 4.0/3.0, 1.0), tri.centroid());
}

TEST(Tri3dGeometry, projectionToPointOnFace)
{
  stk::math::TriGeometry tri({0.0, 0.0, 0.0}, {1.0, 0.0, 0.0}, {0.0, 1.0, 0.0});
  stk::math::TriGeometry triYZ({0.0, 0.0, 0.0}, {0.0, 1.0, 0.0}, {0.0, 0.0, 1.0});

  // At the corners
  EXPECT_VECTOR_NEAR(stk::math::Vector3d(0.0, 0.0, 0.0), tri.closest_proj_on_face(stk::math::Vector3d(0.0, 0.0, 0.0)));
  EXPECT_VECTOR_NEAR(stk::math::Vector3d(1.0, 0.0, 0.0), tri.closest_proj_on_face(stk::math::Vector3d(1.0, 0.0, 0.0)));
  EXPECT_VECTOR_NEAR(stk::math::Vector3d(0.0, 1.0, 0.0), tri.closest_proj_on_face(stk::math::Vector3d(0.0, 1.0, 0.0)));
  EXPECT_VECTOR_NEAR(stk::math::Vector3d(0.0, 0.0, 1.0), triYZ.closest_proj_on_face(stk::math::Vector3d(0.0, 0.0, 1.0)));

  // On the edges
  EXPECT_VECTOR_NEAR(stk::math::Vector3d(0.5, 0.0, 0.0), tri.closest_proj_on_face(stk::math::Vector3d(0.5, 0.0, 0.0)));
  EXPECT_VECTOR_NEAR(stk::math::Vector3d(0.0, 0.5, 0.0), tri.closest_proj_on_face(stk::math::Vector3d(0.0, 0.5, 0.0)));
  EXPECT_VECTOR_NEAR(stk::math::Vector3d(0.5, 0.5, 0.0), tri.closest_proj_on_face(stk::math::Vector3d(0.5, 0.5, 0.0)));
  EXPECT_VECTOR_NEAR(stk::math::Vector3d(0.0, 0.5, 0.5), triYZ.closest_proj_on_face(stk::math::Vector3d(0.0, 0.5, 0.5)));

  // In the middle
  EXPECT_VECTOR_NEAR(stk::math::Vector3d(0.1, 0.1, 0.0), tri.closest_proj_on_face(stk::math::Vector3d(0.1, 0.1, 0.0)));
  EXPECT_VECTOR_NEAR(stk::math::Vector3d(0.8, 0.1, 0.0), tri.closest_proj_on_face(stk::math::Vector3d(0.8, 0.1, 0.0)));
  EXPECT_VECTOR_NEAR(stk::math::Vector3d(0.1, 0.8, 0.0), tri.closest_proj_on_face(stk::math::Vector3d(0.1, 0.8, 0.0)));
  EXPECT_VECTOR_NEAR(stk::math::Vector3d(0.0, 0.1, 0.8), triYZ.closest_proj_on_face(stk::math::Vector3d(0.0, 0.1, 0.8)));
}

TEST(Tri3dGeometry, projectionToPointAboveFace)
{
  stk::math::TriGeometry tri({0.0, 0.0, 0.0}, {1.0, 0.0, 0.0}, {0.0, 1.0, 0.0});
  stk::math::TriGeometry triYZ({0.0, 0.0, 0.0}, {0.0, 1.0, 0.0}, {0.0, 0.0, 1.0});

  // Above the corners
  EXPECT_VECTOR_NEAR(stk::math::Vector3d(0.0, 0.0, 0.0), tri.closest_proj_on_face(stk::math::Vector3d(0.0, 0.0, 1.0)));
  EXPECT_VECTOR_NEAR(stk::math::Vector3d(1.0, 0.0, 0.0), tri.closest_proj_on_face(stk::math::Vector3d(1.0, 0.0, 1.0)));
  EXPECT_VECTOR_NEAR(stk::math::Vector3d(0.0, 1.0, 0.0), tri.closest_proj_on_face(stk::math::Vector3d(0.0, 1.0, 1.0)));
  EXPECT_VECTOR_NEAR(stk::math::Vector3d(0.0, 0.0, 1.0), triYZ.closest_proj_on_face(stk::math::Vector3d(1.0, 0.0, 1.0)));

  // Above the edges
  EXPECT_VECTOR_NEAR(stk::math::Vector3d(0.5, 0.0, 0.0), tri.closest_proj_on_face(stk::math::Vector3d(0.5, 0.0, 1.0)));
  EXPECT_VECTOR_NEAR(stk::math::Vector3d(0.0, 0.5, 0.0), tri.closest_proj_on_face(stk::math::Vector3d(0.0, 0.5, 1.0)));
  EXPECT_VECTOR_NEAR(stk::math::Vector3d(0.5, 0.5, 0.0), tri.closest_proj_on_face(stk::math::Vector3d(0.5, 0.5, 1.0)));
  EXPECT_VECTOR_NEAR(stk::math::Vector3d(0.0, 0.5, 0.5), triYZ.closest_proj_on_face(stk::math::Vector3d(1.0, 0.5, 0.5)));

  // Above the middle
  EXPECT_VECTOR_NEAR(stk::math::Vector3d(0.1, 0.1, 0.0), tri.closest_proj_on_face(stk::math::Vector3d(0.1, 0.1, 1.0)));
  EXPECT_VECTOR_NEAR(stk::math::Vector3d(0.8, 0.1, 0.0), tri.closest_proj_on_face(stk::math::Vector3d(0.8, 0.1, 1.0)));
  EXPECT_VECTOR_NEAR(stk::math::Vector3d(0.1, 0.8, 0.0), tri.closest_proj_on_face(stk::math::Vector3d(0.1, 0.8, 1.0)));
  EXPECT_VECTOR_NEAR(stk::math::Vector3d(0.0, 0.1, 0.8), triYZ.closest_proj_on_face(stk::math::Vector3d(1.0, 0.1, 0.8)));
}

TEST(Tri3dGeometry, projectionToPointOutsideFace_inPlane)
{
  stk::math::TriGeometry tri({0.0, 0.0, 0.0}, {1.0, 0.0, 0.0}, {0.0, 1.0, 0.0});
  stk::math::TriGeometry triYZ({0.0, 0.0, 0.0}, {0.0, 1.0, 0.0}, {0.0, 0.0, 1.0});

  // Along rays projected from left edge
  EXPECT_VECTOR_NEAR(stk::math::Vector3d(0.0, 1.0, 0.0), tri.closest_proj_on_face(stk::math::Vector3d( 0.0, 2.0, 0.0)));
  EXPECT_VECTOR_NEAR(stk::math::Vector3d(0.0, 0.0, 0.0), tri.closest_proj_on_face(stk::math::Vector3d( 0.0,-1.0, 0.0)));
  EXPECT_VECTOR_NEAR(stk::math::Vector3d(0.0, 0.0, 0.0), triYZ.closest_proj_on_face(stk::math::Vector3d( 0.0, 0.0,-1.0)));

  // Along rays projected from bottom edge
  EXPECT_VECTOR_NEAR(stk::math::Vector3d(0.0, 0.0, 0.0), tri.closest_proj_on_face(stk::math::Vector3d(-1.0, 0.0, 0.0)));
  EXPECT_VECTOR_NEAR(stk::math::Vector3d(1.0, 0.0, 0.0), tri.closest_proj_on_face(stk::math::Vector3d( 2.0, 0.0, 0.0)));
  EXPECT_VECTOR_NEAR(stk::math::Vector3d(0.0, 1.0, 0.0), triYZ.closest_proj_on_face(stk::math::Vector3d( 0.0, 2.0, 0.0)));

  // Along rays projected from hypotenuse
  EXPECT_VECTOR_NEAR(stk::math::Vector3d(0.0, 1.0, 0.0), tri.closest_proj_on_face(stk::math::Vector3d(-1.0, 2.0, 0.0)));
  EXPECT_VECTOR_NEAR(stk::math::Vector3d(1.0, 0.0, 0.0), tri.closest_proj_on_face(stk::math::Vector3d( 2.0,-1.0, 0.0)));
  EXPECT_VECTOR_NEAR(stk::math::Vector3d(0.0, 1.0, 0.0), triYZ.closest_proj_on_face(stk::math::Vector3d( 0.0, 2.0,-1.0)));

  // Outside the corners
  EXPECT_VECTOR_NEAR(stk::math::Vector3d(0.0, 1.0, 0.0), tri.closest_proj_on_face(stk::math::Vector3d(-0.5, 2.0, 0.0)));
  EXPECT_VECTOR_NEAR(stk::math::Vector3d(1.0, 0.0, 0.0), tri.closest_proj_on_face(stk::math::Vector3d( 2.0,-0.5, 0.0)));
  EXPECT_VECTOR_NEAR(stk::math::Vector3d(0.0, 0.0, 0.0), tri.closest_proj_on_face(stk::math::Vector3d(-1.0,-1.0, 0.0)));
  EXPECT_VECTOR_NEAR(stk::math::Vector3d(0.0, 0.0, 0.0), triYZ.closest_proj_on_face(stk::math::Vector3d( 0.0,-1.0,-1.0)));

  // Outside the edges
  EXPECT_VECTOR_NEAR(stk::math::Vector3d(0.0, 0.5, 0.0), tri.closest_proj_on_face(stk::math::Vector3d(-1.0, 0.5, 0.0)));
  EXPECT_VECTOR_NEAR(stk::math::Vector3d(0.5, 0.0, 0.0), tri.closest_proj_on_face(stk::math::Vector3d( 0.5,-1.0, 0.0)));
  EXPECT_VECTOR_NEAR(stk::math::Vector3d(0.5, 0.5, 0.0), tri.closest_proj_on_face(stk::math::Vector3d( 1.0, 1.0, 0.0)));
  EXPECT_VECTOR_NEAR(stk::math::Vector3d(0.0, 0.5, 0.5), triYZ.closest_proj_on_face(stk::math::Vector3d( 0.0, 1.0, 1.0)));
}

TEST(Tri3dGeometry, projectionToPointOutsideFace_abovePlane)
{
  stk::math::TriGeometry tri({0.0, 0.0, 0.0}, {1.0, 0.0, 0.0}, {0.0, 1.0, 0.0});

  // Along rays projected from left edge
  EXPECT_VECTOR_NEAR(stk::math::Vector3d(0.0, 1.0, 0.0), tri.closest_proj_on_face(stk::math::Vector3d( 0.0, 2.0, 1.0)));
  EXPECT_VECTOR_NEAR(stk::math::Vector3d(0.0, 0.0, 0.0), tri.closest_proj_on_face(stk::math::Vector3d( 0.0,-1.0, 1.0)));

  // Along rays projected from bottom edge
  EXPECT_VECTOR_NEAR(stk::math::Vector3d(0.0, 0.0, 0.0), tri.closest_proj_on_face(stk::math::Vector3d(-1.0, 0.0, 1.0)));
  EXPECT_VECTOR_NEAR(stk::math::Vector3d(1.0, 0.0, 0.0), tri.closest_proj_on_face(stk::math::Vector3d( 2.0, 0.0, 1.0)));

  // Along rays projected from hypotenuse
  EXPECT_VECTOR_NEAR(stk::math::Vector3d(0.0, 1.0, 0.0), tri.closest_proj_on_face(stk::math::Vector3d(-1.0, 2.0, 1.0)));
  EXPECT_VECTOR_NEAR(stk::math::Vector3d(1.0, 0.0, 0.0), tri.closest_proj_on_face(stk::math::Vector3d( 2.0,-1.0, 1.0)));

  // Outside the corners
  EXPECT_VECTOR_NEAR(stk::math::Vector3d(0.0, 1.0, 0.0), tri.closest_proj_on_face(stk::math::Vector3d(-0.5, 2.0, 1.0)));
  EXPECT_VECTOR_NEAR(stk::math::Vector3d(1.0, 0.0, 0.0), tri.closest_proj_on_face(stk::math::Vector3d( 2.0,-0.5, 1.0)));
  EXPECT_VECTOR_NEAR(stk::math::Vector3d(0.0, 0.0, 0.0), tri.closest_proj_on_face(stk::math::Vector3d(-1.0,-1.0, 1.0)));

  // Outside the edges
  EXPECT_VECTOR_NEAR(stk::math::Vector3d(0.0, 0.5, 0.0), tri.closest_proj_on_face(stk::math::Vector3d(-1.0, 0.5, 1.0)));
  EXPECT_VECTOR_NEAR(stk::math::Vector3d(0.5, 0.0, 0.0), tri.closest_proj_on_face(stk::math::Vector3d( 0.5,-1.0, 1.0)));
  EXPECT_VECTOR_NEAR(stk::math::Vector3d(0.5, 0.5, 0.0), tri.closest_proj_on_face(stk::math::Vector3d( 1.0, 1.0, 1.0)));
}

TEST(Tri3dGeometry, distanceToPoint)
{
  stk::math::TriGeometry tri({0.0, 0.0, 0.0}, {1.0, 0.0, 0.0}, {0.0, 1.0, 0.0});

  EXPECT_DOUBLE_EQ(0.0                , tri.min_distance_to_point(stk::math::Vector3d( 0.0, 0.0, 0.0)));
  EXPECT_DOUBLE_EQ(0.0                , tri.min_distance_to_point(stk::math::Vector3d( 0.1, 0.1, 0.0)));
  EXPECT_DOUBLE_EQ(std::sqrt(2.0)/2.0 , tri.min_distance_to_point(stk::math::Vector3d( 1.0, 1.0, 0.0)));
  EXPECT_DOUBLE_EQ(1.0                , tri.min_distance_to_point(stk::math::Vector3d(-1.0, 0.0, 0.0)));
  EXPECT_DOUBLE_EQ(std::sqrt(1.0+0.25), tri.min_distance_to_point(stk::math::Vector3d(-0.5, 2.0, 0.0)));

  EXPECT_DOUBLE_EQ(1.0                , tri.min_distance_to_point(stk::math::Vector3d( 0.0, 0.0, 1.0)));
  EXPECT_DOUBLE_EQ(1.0                , tri.min_distance_to_point(stk::math::Vector3d( 0.1, 0.1, 1.0)));
  EXPECT_DOUBLE_EQ(std::sqrt(0.5+1.0) , tri.min_distance_to_point(stk::math::Vector3d( 1.0, 1.0, 1.0)));
  EXPECT_DOUBLE_EQ(std::sqrt(1.0+1.0) , tri.min_distance_to_point(stk::math::Vector3d(-1.0, 0.0, 1.0)));
  EXPECT_DOUBLE_EQ(std::sqrt(1.25+1.0), tri.min_distance_to_point(stk::math::Vector3d(-0.5, 2.0, 1.0)));
}

TEST(Tri3dGeometry, areNodesCloseToSide)
{
  stk::math::Vector3d notClose(2.0, 2.0, 0.0);

  stk::math::Vector3d node1(0.0, 0.0, 0.0);
  stk::math::Vector3d node2(1.0, 0.0, 0.0);
  stk::math::Vector3d node3(0.0, 1.0, 0.0);

  stk::math::TriGeometry tri(node1, node2, node3);

  stk::math::Vector3d centroid = tri.centroid();

  stk::math::QuadGeometry quadNode1 = small_quad(node1, 0.01);
  stk::math::QuadGeometry quadNode2 = small_quad(node2, 0.01);
  stk::math::QuadGeometry quadNode3 = small_quad(node3, 0.01);
  stk::math::QuadGeometry quadCentroid = small_quad(centroid, 0.01);
  stk::math::QuadGeometry quadNotClose = small_quad(notClose, 0.01);

  EXPECT_TRUE(tri.are_nodes_close_to_side(quadNode1, 0.01));
  EXPECT_TRUE(tri.are_nodes_close_to_side(quadNode2, 0.01));
  EXPECT_TRUE(tri.are_nodes_close_to_side(quadNode3, 0.01));
  EXPECT_TRUE(tri.are_nodes_close_to_side(quadCentroid, 0.01));
  EXPECT_FALSE(tri.are_nodes_close_to_side(quadNotClose, 0.01));
}

TEST(Quad3dGeometry, constructWithVectorsAndQueryNodes)
{
  stk::math::QuadGeometry quad({1, 2, 3}, {4, 5, 6}, {7, 8, 9}, {10, 11, 12});

  EXPECT_EQ(stk::math::Vector3d( 1,  2,  3), quad.node(0));
  EXPECT_EQ(stk::math::Vector3d( 4,  5,  6), quad.node(1));
  EXPECT_EQ(stk::math::Vector3d( 7,  8,  9), quad.node(2));
  EXPECT_EQ(stk::math::Vector3d(10, 11, 12), quad.node(3));
}

TEST(Quad3dGeometry, centroid)
{
  stk::math::QuadGeometry quad({0.0, 0.0, 0.0}, {1.0, 0.0, 0.0}, {1.0, 1.5, 2.0}, {0.0, 1.5, 2.0});

  EXPECT_EQ(stk::math::Vector3d(0.5, 0.75, 1.0), quad.centroid());
}

TEST(Quad3dGeometry, projectionToPointOnFace)
{
  stk::math::QuadGeometry quad({0.0, 0.0, 0.0}, {1.0, 0.0, 0.0}, {1.0, 1.0, 0.0}, {0.0, 1.0, 0.0});
  stk::math::QuadGeometry quadYZ({0.0, 0.0, 0.0}, {0.0, 1.0, 0.0}, {0.0, 1.0, 1.0}, {0.0, 0.0, 1.0});

  // At the corners
  EXPECT_VECTOR_NEAR(stk::math::Vector3d(0.0, 0.0, 0.0), quad.closest_proj_on_face(stk::math::Vector3d(0.0, 0.0, 0.0)));
  EXPECT_VECTOR_NEAR(stk::math::Vector3d(1.0, 0.0, 0.0), quad.closest_proj_on_face(stk::math::Vector3d(1.0, 0.0, 0.0)));
  EXPECT_VECTOR_NEAR(stk::math::Vector3d(1.0, 1.0, 0.0), quad.closest_proj_on_face(stk::math::Vector3d(1.0, 1.0, 0.0)));
  EXPECT_VECTOR_NEAR(stk::math::Vector3d(0.0, 1.0, 0.0), quad.closest_proj_on_face(stk::math::Vector3d(0.0, 1.0, 0.0)));
  EXPECT_VECTOR_NEAR(stk::math::Vector3d(0.0, 0.0, 1.0), quadYZ.closest_proj_on_face(stk::math::Vector3d(0.0, 0.0, 1.0)));

  // On the edges
  EXPECT_VECTOR_NEAR(stk::math::Vector3d(0.5, 0.0, 0.0), quad.closest_proj_on_face(stk::math::Vector3d(0.5, 0.0, 0.0)));
  EXPECT_VECTOR_NEAR(stk::math::Vector3d(1.0, 0.5, 0.0), quad.closest_proj_on_face(stk::math::Vector3d(1.0, 0.5, 0.0)));
  EXPECT_VECTOR_NEAR(stk::math::Vector3d(0.5, 1.0, 0.0), quad.closest_proj_on_face(stk::math::Vector3d(0.5, 1.0, 0.0)));
  EXPECT_VECTOR_NEAR(stk::math::Vector3d(0.0, 0.5, 0.0), quad.closest_proj_on_face(stk::math::Vector3d(0.0, 0.5, 0.0)));
  EXPECT_VECTOR_NEAR(stk::math::Vector3d(0.0, 0.0, 0.5), quadYZ.closest_proj_on_face(stk::math::Vector3d(0.0, 0.0, 0.5)));

  // In the middle
  EXPECT_VECTOR_NEAR(stk::math::Vector3d(0.1, 0.1, 0.0), quad.closest_proj_on_face(stk::math::Vector3d(0.1, 0.1, 0.0)));
  EXPECT_VECTOR_NEAR(stk::math::Vector3d(0.9, 0.1, 0.0), quad.closest_proj_on_face(stk::math::Vector3d(0.9, 0.1, 0.0)));
  EXPECT_VECTOR_NEAR(stk::math::Vector3d(0.9, 0.9, 0.0), quad.closest_proj_on_face(stk::math::Vector3d(0.9, 0.9, 0.0)));
  EXPECT_VECTOR_NEAR(stk::math::Vector3d(0.1, 0.9, 0.0), quad.closest_proj_on_face(stk::math::Vector3d(0.1, 0.9, 0.0)));
  EXPECT_VECTOR_NEAR(stk::math::Vector3d(0.5, 0.5, 0.0), quad.closest_proj_on_face(stk::math::Vector3d(0.5, 0.5, 0.0)));
  EXPECT_VECTOR_NEAR(stk::math::Vector3d(0.0, 0.5, 0.5), quadYZ.closest_proj_on_face(stk::math::Vector3d(0.0, 0.5, 0.5)));
}

TEST(Quad3dGeometry, projectionToPointAboveFace)
{
  stk::math::QuadGeometry quad({0.0, 0.0, 0.0}, {1.0, 0.0, 0.0}, {1.0, 1.0, 0.0}, {0.0, 1.0, 0.0});

  // At the corners
  EXPECT_VECTOR_NEAR(stk::math::Vector3d(0.0, 0.0, 0.0), quad.closest_proj_on_face(stk::math::Vector3d(0.0, 0.0, 1.0)));
  EXPECT_VECTOR_NEAR(stk::math::Vector3d(1.0, 0.0, 0.0), quad.closest_proj_on_face(stk::math::Vector3d(1.0, 0.0, 1.0)));
  EXPECT_VECTOR_NEAR(stk::math::Vector3d(1.0, 1.0, 0.0), quad.closest_proj_on_face(stk::math::Vector3d(1.0, 1.0, 1.0)));
  EXPECT_VECTOR_NEAR(stk::math::Vector3d(0.0, 1.0, 0.0), quad.closest_proj_on_face(stk::math::Vector3d(0.0, 1.0, 1.0)));

  // On the edges
  EXPECT_VECTOR_NEAR(stk::math::Vector3d(0.5, 0.0, 0.0), quad.closest_proj_on_face(stk::math::Vector3d(0.5, 0.0, 1.0)));
  EXPECT_VECTOR_NEAR(stk::math::Vector3d(1.0, 0.5, 0.0), quad.closest_proj_on_face(stk::math::Vector3d(1.0, 0.5, 1.0)));
  EXPECT_VECTOR_NEAR(stk::math::Vector3d(0.5, 1.0, 0.0), quad.closest_proj_on_face(stk::math::Vector3d(0.5, 1.0, 1.0)));
  EXPECT_VECTOR_NEAR(stk::math::Vector3d(0.0, 0.5, 0.0), quad.closest_proj_on_face(stk::math::Vector3d(0.0, 0.5, 1.0)));

  // In the middle
  EXPECT_VECTOR_NEAR(stk::math::Vector3d(0.1, 0.1, 0.0), quad.closest_proj_on_face(stk::math::Vector3d(0.1, 0.1, 1.0)));
  EXPECT_VECTOR_NEAR(stk::math::Vector3d(0.9, 0.1, 0.0), quad.closest_proj_on_face(stk::math::Vector3d(0.9, 0.1, 1.0)));
  EXPECT_VECTOR_NEAR(stk::math::Vector3d(0.9, 0.9, 0.0), quad.closest_proj_on_face(stk::math::Vector3d(0.9, 0.9, 1.0)));
  EXPECT_VECTOR_NEAR(stk::math::Vector3d(0.1, 0.9, 0.0), quad.closest_proj_on_face(stk::math::Vector3d(0.1, 0.9, 1.0)));
  EXPECT_VECTOR_NEAR(stk::math::Vector3d(0.5, 0.5, 0.0), quad.closest_proj_on_face(stk::math::Vector3d(0.5, 0.5, 1.0)));
}

TEST(Quad3dGeometry, projectionToPointOutsideFace_inPlane)
{
  stk::math::QuadGeometry quad({0.0, 0.0, 0.0}, {1.0, 0.0, 0.0}, {1.0, 1.0, 0.0}, {0.0, 1.0, 0.0});

  // Along rays projected from bottom edge
  EXPECT_VECTOR_NEAR(stk::math::Vector3d(0.0, 0.0, 0.0), quad.closest_proj_on_face(stk::math::Vector3d(-1.0, 0.0, 0.0)));
  EXPECT_VECTOR_NEAR(stk::math::Vector3d(1.0, 0.0, 0.0), quad.closest_proj_on_face(stk::math::Vector3d( 2.0, 0.0, 0.0)));

  // Along rays projected from right edge
  EXPECT_VECTOR_NEAR(stk::math::Vector3d(1.0, 0.0, 0.0), quad.closest_proj_on_face(stk::math::Vector3d( 1.0,-1.0, 0.0)));
  EXPECT_VECTOR_NEAR(stk::math::Vector3d(1.0, 1.0, 0.0), quad.closest_proj_on_face(stk::math::Vector3d( 1.0, 2.0, 0.0)));

  // Along rays projected from top edge
  EXPECT_VECTOR_NEAR(stk::math::Vector3d(0.0, 1.0, 0.0), quad.closest_proj_on_face(stk::math::Vector3d(-1.0, 1.0, 0.0)));
  EXPECT_VECTOR_NEAR(stk::math::Vector3d(1.0, 1.0, 0.0), quad.closest_proj_on_face(stk::math::Vector3d( 2.0, 1.0, 0.0)));

  // Along rays projected from left edge
  EXPECT_VECTOR_NEAR(stk::math::Vector3d(0.0, 0.0, 0.0), quad.closest_proj_on_face(stk::math::Vector3d( 0.0,-1.0, 0.0)));
  EXPECT_VECTOR_NEAR(stk::math::Vector3d(0.0, 1.0, 0.0), quad.closest_proj_on_face(stk::math::Vector3d( 0.0, 2.0, 0.0)));

  // Outside the corners
  EXPECT_VECTOR_NEAR(stk::math::Vector3d(0.0, 0.0, 0.0), quad.closest_proj_on_face(stk::math::Vector3d(-1.0,-1.0, 0.0)));
  EXPECT_VECTOR_NEAR(stk::math::Vector3d(1.0, 0.0, 0.0), quad.closest_proj_on_face(stk::math::Vector3d( 2.0,-1.0, 0.0)));
  EXPECT_VECTOR_NEAR(stk::math::Vector3d(1.0, 1.0, 0.0), quad.closest_proj_on_face(stk::math::Vector3d( 2.0, 2.0, 0.0)));
  EXPECT_VECTOR_NEAR(stk::math::Vector3d(0.0, 1.0, 0.0), quad.closest_proj_on_face(stk::math::Vector3d(-1.0, 2.0, 0.0)));

  // Outside the edges
  EXPECT_VECTOR_NEAR(stk::math::Vector3d(0.5, 0.0, 0.0), quad.closest_proj_on_face(stk::math::Vector3d( 0.5,-1.0, 0.0)));
  EXPECT_VECTOR_NEAR(stk::math::Vector3d(1.0, 0.5, 0.0), quad.closest_proj_on_face(stk::math::Vector3d( 2.0, 0.5, 0.0)));
  EXPECT_VECTOR_NEAR(stk::math::Vector3d(0.5, 1.0, 0.0), quad.closest_proj_on_face(stk::math::Vector3d( 0.5, 2.0, 0.0)));
  EXPECT_VECTOR_NEAR(stk::math::Vector3d(0.0, 0.5, 0.0), quad.closest_proj_on_face(stk::math::Vector3d(-1.0, 0.5, 0.0)));
}

TEST(Quad3dGeometry, projectionToPointOutsideFace_abovePlane)
{
  stk::math::QuadGeometry quad({0.0, 0.0, 0.0}, {1.0, 0.0, 0.0}, {1.0, 1.0, 0.0}, {0.0, 1.0, 0.0});

  // Along rays projected from bottom edge
  EXPECT_VECTOR_NEAR(stk::math::Vector3d(0.0, 0.0, 0.0), quad.closest_proj_on_face(stk::math::Vector3d(-1.0, 0.0, 1.0)));
  EXPECT_VECTOR_NEAR(stk::math::Vector3d(1.0, 0.0, 0.0), quad.closest_proj_on_face(stk::math::Vector3d( 2.0, 0.0, 1.0)));

  // Along rays projected from right edge
  EXPECT_VECTOR_NEAR(stk::math::Vector3d(1.0, 0.0, 0.0), quad.closest_proj_on_face(stk::math::Vector3d( 1.0,-1.0, 1.0)));
  EXPECT_VECTOR_NEAR(stk::math::Vector3d(1.0, 1.0, 0.0), quad.closest_proj_on_face(stk::math::Vector3d( 1.0, 2.0, 1.0)));

  // Along rays projected from top edge
  EXPECT_VECTOR_NEAR(stk::math::Vector3d(0.0, 1.0, 0.0), quad.closest_proj_on_face(stk::math::Vector3d(-1.0, 1.0, 1.0)));
  EXPECT_VECTOR_NEAR(stk::math::Vector3d(1.0, 1.0, 0.0), quad.closest_proj_on_face(stk::math::Vector3d( 2.0, 1.0, 1.0)));

  // Along rays projected from left edge
  EXPECT_VECTOR_NEAR(stk::math::Vector3d(0.0, 0.0, 0.0), quad.closest_proj_on_face(stk::math::Vector3d( 0.0,-1.0, 1.0)));
  EXPECT_VECTOR_NEAR(stk::math::Vector3d(0.0, 1.0, 0.0), quad.closest_proj_on_face(stk::math::Vector3d( 0.0, 2.0, 1.0)));

  // Outside the corners
  EXPECT_VECTOR_NEAR(stk::math::Vector3d(0.0, 0.0, 0.0), quad.closest_proj_on_face(stk::math::Vector3d(-1.0,-1.0, 1.0)));
  EXPECT_VECTOR_NEAR(stk::math::Vector3d(1.0, 0.0, 0.0), quad.closest_proj_on_face(stk::math::Vector3d( 2.0,-1.0, 1.0)));
  EXPECT_VECTOR_NEAR(stk::math::Vector3d(1.0, 1.0, 0.0), quad.closest_proj_on_face(stk::math::Vector3d( 2.0, 2.0, 1.0)));
  EXPECT_VECTOR_NEAR(stk::math::Vector3d(0.0, 1.0, 0.0), quad.closest_proj_on_face(stk::math::Vector3d(-1.0, 2.0, 1.0)));

  // Outside the edges
  EXPECT_VECTOR_NEAR(stk::math::Vector3d(0.5, 0.0, 0.0), quad.closest_proj_on_face(stk::math::Vector3d( 0.5,-1.0, 1.0)));
  EXPECT_VECTOR_NEAR(stk::math::Vector3d(1.0, 0.5, 0.0), quad.closest_proj_on_face(stk::math::Vector3d( 2.0, 0.5, 1.0)));
  EXPECT_VECTOR_NEAR(stk::math::Vector3d(0.5, 1.0, 0.0), quad.closest_proj_on_face(stk::math::Vector3d( 0.5, 2.0, 1.0)));
  EXPECT_VECTOR_NEAR(stk::math::Vector3d(0.0, 0.5, 0.0), quad.closest_proj_on_face(stk::math::Vector3d(-1.0, 0.5, 1.0)));
}

TEST(Quad3dGeometry, distanceToPoint)
{
  stk::math::QuadGeometry quad({0.0, 0.0, 0.0}, {1.0, 0.0, 0.0}, {1.0, 1.0, 0.0}, {0.0, 2.0, 0.0});

  EXPECT_DOUBLE_EQ(0.0                , quad.min_distance_to_point(stk::math::Vector3d( 0.0, 0.0, 0.0)));
  EXPECT_DOUBLE_EQ(0.0                , quad.min_distance_to_point(stk::math::Vector3d( 0.0, 2.0, 0.0)));
  EXPECT_DOUBLE_EQ(0.0                , quad.min_distance_to_point(stk::math::Vector3d( 0.5, 0.5, 0.0)));
  EXPECT_DOUBLE_EQ(0.0                , quad.min_distance_to_point(stk::math::Vector3d( 0.3, 1.3, 0.0)));
  EXPECT_DOUBLE_EQ(1.0                , quad.min_distance_to_point(stk::math::Vector3d(-1.0, 1.0, 0.0)));
  EXPECT_DOUBLE_EQ(std::sqrt(2.0)/2.0 , quad.min_distance_to_point(stk::math::Vector3d( 1.0, 2.0, 0.0)));
  EXPECT_DOUBLE_EQ(1.0                , quad.min_distance_to_point(stk::math::Vector3d( 2.0, 0.0, 0.0)));

  EXPECT_DOUBLE_EQ(1.0                , quad.min_distance_to_point(stk::math::Vector3d( 0.0, 0.0, 1.0)));
  EXPECT_DOUBLE_EQ(1.0                , quad.min_distance_to_point(stk::math::Vector3d( 0.0, 2.0, 1.0)));
  EXPECT_DOUBLE_EQ(1.0                , quad.min_distance_to_point(stk::math::Vector3d( 0.5, 0.5, 1.0)));
  EXPECT_DOUBLE_EQ(1.0                , quad.min_distance_to_point(stk::math::Vector3d( 0.3, 1.3, 1.0)));
  EXPECT_DOUBLE_EQ(std::sqrt(2.0)     , quad.min_distance_to_point(stk::math::Vector3d(-1.0, 1.0, 1.0)));
  EXPECT_DOUBLE_EQ(std::sqrt(1.0+0.5) , quad.min_distance_to_point(stk::math::Vector3d( 1.0, 2.0, 1.0)));
  EXPECT_DOUBLE_EQ(std::sqrt(2.0)     , quad.min_distance_to_point(stk::math::Vector3d( 2.0, 0.0, 1.0)));
}

TEST(Quad3dGeometry, areNodesCloseToSide)
{
  stk::math::Vector3d notClose(2.0, 2.0, 0.0);

  stk::math::Vector3d node1(0.0, 0.0, 0.0);
  stk::math::Vector3d node2(1.0, 0.0, 0.0);
  stk::math::Vector3d node3(1.0, 1.0, 0.0);
  stk::math::Vector3d node4(0.0, 1.0, 0.0);

  stk::math::QuadGeometry quad(node1, node2, node3, node4);

  stk::math::Vector3d centroid = quad.centroid();

  stk::math::QuadGeometry quadNode1 = small_quad(node1, 0.01);
  stk::math::QuadGeometry quadNode2 = small_quad(node2, 0.01);
  stk::math::QuadGeometry quadNode3 = small_quad(node3, 0.01);
  stk::math::QuadGeometry quadNode4 = small_quad(node4, 0.01);
  stk::math::QuadGeometry quadCentroid = small_quad(centroid, 0.01);
  stk::math::QuadGeometry quadNotClose = small_quad(notClose, 0.01);

  EXPECT_TRUE(quad.are_nodes_close_to_side(quadNode1, 0.01));
  EXPECT_TRUE(quad.are_nodes_close_to_side(quadNode2, 0.01));
  EXPECT_TRUE(quad.are_nodes_close_to_side(quadNode3, 0.01));
  EXPECT_TRUE(quad.are_nodes_close_to_side(quadNode4, 0.01));
  EXPECT_TRUE(quad.are_nodes_close_to_side(quadCentroid, 0.01));
  EXPECT_FALSE(quad.are_nodes_close_to_side(quadNotClose, 0.01));
}

}
