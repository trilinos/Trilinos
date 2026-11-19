// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include <gtest/gtest.h>

#include <Akri_AnalyticSurf.hpp>
#include <Akri_BoundingBox.hpp>
#include <Akri_MeshFacetedSurface.hpp>
#include <Akri_Transformation.hpp>
#include <Akri_TransformedSurface.hpp>
#include <Akri_UnitTestUtils.hpp>

#include <stk_math/StkVector.hpp>
#include <stk_util/parallel/Parallel.hpp>

#include <array>
#include <cmath>
#include <memory>
#include <vector>

namespace krino {

namespace {

const double quarterTurn = 0.5*M_PI;
const double transformedTol = 1.e-12;
const BoundingBox anyBox(stk::math::Vector3d(-10.,-10.,-10.), stk::math::Vector3d(10.,10.,10.));

void expect_sphere_centered_at(const Surface & sphere, const double radius, const stk::math::Vector3d & center)
{
  EXPECT_NEAR(-radius, sphere.point_signed_distance(center), transformedTol) << center;
}

// A motion with all of its parts: an initial rotation and displacement and velocities.
Transformation general_motion()
{
  Transformation motion;
  motion.set_reference_point(stk::math::Vector3d(0.5,-0.3,0.2));
  motion.set_initial_displacement(stk::math::Vector3d(0.1,0.2,-0.3));
  motion.set_initial_rotation(stk::math::Vector3d(0.3,-0.7,1.1));
  motion.set_translational_velocity(stk::math::Vector3d(-0.1,0.05,0.2));
  motion.set_rotational_velocity(stk::math::Vector3d(0.,0.,quarterTurn));
  return motion;
}

const std::vector<stk::math::Vector3d> samplePoints{{0.,0.,0.}, {1.,-2.,0.5}, {-0.7,0.4,1.3}, {2.,2.,-1.}, {0.3,0.1,-0.2}};

}

TEST(TransformedSurface, sphereWithInitialDisplacement_displacedAtTimeZero)
{
  const stk::math::Vector3d center(1.,0.,0.);
  const stk::math::Vector3d displacement(0.,2.,0.);
  Transformation motion;
  motion.set_initial_displacement(displacement);
  TransformedSurface sphere(std::make_shared<Sphere>(center, 0.5), motion);
  expect_sphere_centered_at(sphere, 0.5, center+displacement);
  sphere.prepare_to_compute(0., anyBox, 0.);
  expect_sphere_centered_at(sphere, 0.5, center+displacement);
}

TEST(TransformedSurface, sphereWithInitialRotationAboutReferencePoint_rotatedAtTimeZero)
{
  Transformation motion;
  motion.set_reference_point(stk::math::Vector3d(1.,1.,0.));
  motion.set_initial_rotation(stk::math::Vector3d(0.,0.,quarterTurn));
  TransformedSurface sphere(std::make_shared<Sphere>(stk::math::Vector3d(2.,1.,0.), 0.5), motion);
  sphere.prepare_to_compute(0., anyBox, 0.);
  expect_sphere_centered_at(sphere, 0.5, stk::math::Vector3d(1.,2.,0.));
}

TEST(TransformedSurface, preparedTwiceAtSameTimeOrBackInTime_whereTheMotionPutsItThen)
{
  const stk::math::Vector3d center(1.,0.,0.);
  const stk::math::Vector3d velocity(0.,1.,0.);
  Transformation motion;
  motion.set_translational_velocity(velocity);
  TransformedSurface sphere(std::make_shared<Sphere>(center, 0.5), motion);
  for (const double time : {1., 1., 3., 2.})
  {
    sphere.prepare_to_compute(time, anyBox, 0.);
    expect_sphere_centered_at(sphere, 0.5, center+time*velocity);
  }
}

TEST(TransformedSurface, twoSurfacesWithTheSameMotion_bothMoved)
{
  Transformation motion;
  motion.set_translational_velocity(stk::math::Vector3d(0.,1.,0.));
  TransformedSurface sphere0(std::make_shared<Sphere>(stk::math::Vector3d(1.,0.,0.), 0.5), motion);
  TransformedSurface sphere1(std::make_shared<Sphere>(stk::math::Vector3d(-1.,0.,0.), 0.5), motion);
  for (const double time : {1., 2.})
  {
    sphere0.prepare_to_compute(time, anyBox, 0.);
    sphere1.prepare_to_compute(time, anyBox, 0.);
  }
  expect_sphere_centered_at(sphere0, 0.5, stk::math::Vector3d(1.,2.,0.));
  expect_sphere_centered_at(sphere1, 0.5, stk::math::Vector3d(-1.,2.,0.));
}

TEST(TransformedSurface, sphereTurningAboutDisplacedReferencePoint_turnsAboutIt)
{
  Transformation motion;
  motion.set_reference_point(stk::math::Vector3d(1.,0.,0.));
  motion.set_initial_displacement(stk::math::Vector3d(0.,3.,0.));
  motion.set_rotational_velocity(stk::math::Vector3d(0.,0.,quarterTurn));
  TransformedSurface sphere(std::make_shared<Sphere>(stk::math::Vector3d(2.,0.,0.), 0.5), motion);
  sphere.prepare_to_compute(1., anyBox, 0.);
  // Displaced to (2,3,0), then a quarter turn about the displaced reference point (1,3,0).
  expect_sphere_centered_at(sphere, 0.5, stk::math::Vector3d(1.,4.,0.));
}

TEST(TransformedSurface, changesWithTime)
{
  const TransformedSurface sphere(std::make_shared<Sphere>(stk::math::Vector3d::ZERO, 0.5), Transformation());
  EXPECT_TRUE(sphere.changes_with_time());
}

TEST(TransformedSurface, innerSurface_notMoved)
{
  const stk::math::Vector3d center(1.,0.,0.);
  auto inner = std::make_shared<Sphere>(center, 0.5);
  Transformation motion;
  motion.set_translational_velocity(stk::math::Vector3d(0.,1.,0.));
  TransformedSurface sphere(inner, motion);
  sphere.prepare_to_compute(2., anyBox, 0.);
  expect_sphere_centered_at(*inner, 0.5, center);
}

// The queries of a moving cylinder are those of the cylinder with its end points where the motion
// puts them, including the normals, which turn with it.
TEST(TransformedSurface, cylinder_sameQueriesAsCylinderWhereTheMotionPutsIt)
{
  const stk::math::Vector3d end0(-0.5,0.,0.), end1(0.5,0.,0.);
  const double radius = 0.3;
  const Transformation motion = general_motion();
  TransformedSurface moving(std::make_shared<Cylinder>(end0, end1, radius), motion);
  const double time = 0.7;
  moving.prepare_to_compute(time, anyBox, 0.);
  const RigidTransformation placement = motion.rigid_transformation_at(time);
  Cylinder placed(placement.transform_point(end0), placement.transform_point(end1), radius);
  placed.prepare_to_compute(time, anyBox, 0.);

  EXPECT_TRUE(moving.can_approximate_closest_point_normal());
  for (auto & x : samplePoints)
  {
    EXPECT_NEAR(placed.point_signed_distance(x), moving.point_signed_distance(x), transformedTol) << x;
    expect_eq_absolute(placed.closest_point_normal(x), moving.closest_point_normal(x), transformedTol);
  }
  const stk::math::Vector3d pt0 = placement.transform_point(stk::math::Vector3d(0.,0.,0.));
  const stk::math::Vector3d pt1 = placement.transform_point(stk::math::Vector3d(0.,0.,1.));
  const auto [placedSign, placedLoc] = placed.compute_intersection_with_segment(pt0, pt1, 1.e-10);
  const auto [movingSign, movingLoc] = moving.compute_intersection_with_segment(pt0, pt1, 1.e-10);
  EXPECT_EQ(1, movingSign);
  EXPECT_EQ(placedSign, movingSign);
  EXPECT_NEAR(placedLoc, movingLoc, 1.e-9);
  EXPECT_NEAR(0.3, movingLoc, 1.e-9);
}

// The bounding box of a moving surface is around where the motion puts it, and it intersects only
// the boxes near there.
TEST(TransformedSurface, sphere_boundingBoxWhereTheMotionPutsIt)
{
  const stk::math::Vector3d center(1.,0.,0.);
  Transformation motion;
  motion.set_translational_velocity(stk::math::Vector3d(0.,5.,0.));
  TransformedSurface sphere(std::make_shared<Sphere>(center, 0.5), motion);
  sphere.prepare_to_compute(1., anyBox, 0.);

  BoundingBox bbox;
  sphere.insert_into(bbox);
  const double bboxTol = 1.e-5; // BoundingBox is single precision
  for (int i=0; i<3; ++i)
  {
    EXPECT_NEAR((stk::math::Vector3d(0.5,4.5,-0.5))[i], bbox.get_min()[i], bboxTol);
    EXPECT_NEAR((stk::math::Vector3d(1.5,5.5,0.5))[i], bbox.get_max()[i], bboxTol);
  }
  EXPECT_TRUE(sphere.does_intersect(BoundingBox(stk::math::Vector3d(0.9,4.9,-0.1), stk::math::Vector3d(1.1,5.1,0.1))));
  EXPECT_FALSE(sphere.does_intersect(BoundingBox(stk::math::Vector3d(0.9,-0.1,-0.1), stk::math::Vector3d(1.1,0.1,0.1))));
}

// A surface that is everywhere, such as a plane, stays everywhere.
TEST(TransformedSurface, plane_boundingBoxStaysTheEntireDomain)
{
  Transformation motion;
  motion.set_initial_rotation(stk::math::Vector3d(0.3,-0.7,1.1));
  TransformedSurface plane(std::make_shared<Plane>(stk::math::Vector3d(0.,0.,1.), 0.), motion);
  BoundingBox bbox;
  plane.insert_into(bbox);
  EXPECT_TRUE(bbox.contains(BoundingBox::ENTIRE_DOMAIN));
}

TEST(TransformedSurface, planeAtFractionOfExtent_throws)
{
  EXPECT_ANY_THROW(TransformedSurface(Plane::at_fraction_of_extent(stk::math::Vector3d(0.,0.,1.), 0.5, 1.), Transformation()));
}

// The sharp features of a moving cuboid are those of the cuboid where the motion puts it.
TEST(TransformedSurface, cuboid_sameSharpFeaturesAsCuboidWhereTheMotionPutsIt)
{
  const stk::math::Vector3d center(0.1,-0.2,0.3);
  const stk::math::Vector3d dimensions(1.,0.6,0.8);
  const stk::math::Vector3d initialRotation(0.3,-0.7,1.1);
  Transformation motion;
  motion.set_reference_point(stk::math::Vector3d(0.5,-0.3,0.2));
  motion.set_initial_displacement(stk::math::Vector3d(0.1,0.2,-0.3));
  motion.set_initial_rotation(initialRotation);
  TransformedSurface moving(std::make_shared<Cuboid>(center, dimensions), motion);
  moving.prepare_to_compute(0., anyBox, 0.);
  const Cuboid placed(motion.rigid_transformation_at(0.).transform_point(center), dimensions, initialRotation);

  EXPECT_TRUE(moving.has_edges());
  EXPECT_TRUE(moving.has_corners());
  for (auto & x : samplePoints)
  {
    EXPECT_NEAR(placed.point_signed_distance(x), moving.point_signed_distance(x), transformedTol) << x;
  }

  std::vector<stk::math::Vector3d> placedParamCoords, movingParamCoords;
  size_t numIntersections = 0;
  for (unsigned i=0; i<8; ++i)
  {
    const stk::math::Vector3d corner = placed.vertex_location(i);
    const std::array<stk::math::Vector3d,4> tet{{corner+stk::math::Vector3d(-0.05,-0.05,-0.05),
        corner+stk::math::Vector3d(0.15,-0.05,-0.05), corner+stk::math::Vector3d(-0.05,0.15,-0.05), corner+stk::math::Vector3d(-0.05,-0.05,0.15)}};
    placed.fill_tetrahedon_intersection_parametric_coordinates(tet, placedParamCoords);
    moving.fill_tetrahedon_intersection_parametric_coordinates(tet, movingParamCoords);
    ASSERT_EQ(placedParamCoords.size(), movingParamCoords.size());
    for (size_t j=0; j<placedParamCoords.size(); ++j)
      expect_eq_absolute(placedParamCoords[j], movingParamCoords[j], transformedTol);
    numIntersections += movingParamCoords.size();

    const std::array<stk::math::Vector3d,3> tri{{tet[0], tet[1], tet[2]}};
    placed.fill_triangle_intersection_parametric_coordinates(tri, placedParamCoords);
    moving.fill_triangle_intersection_parametric_coordinates(tri, movingParamCoords);
    ASSERT_EQ(placedParamCoords.size(), movingParamCoords.size());
    for (size_t j=0; j<placedParamCoords.size(); ++j)
      expect_eq_absolute(placedParamCoords[j], movingParamCoords[j], transformedTol);
    numIntersections += movingParamCoords.size();
  }
  EXPECT_LE(8u, numIntersections);
}

// The crossings of a segment with a moving cuboid are those with the cuboid where the motion puts
// it.
TEST(TransformedSurface, cuboid_sameSegmentCrossingsAsCuboidWhereTheMotionPutsIt)
{
  const stk::math::Vector3d center(0.1, -0.2, 0.3);
  const stk::math::Vector3d dimensions(1., 0.6, 0.8);
  const stk::math::Vector3d initialRotation(0.3, -0.7, 1.1);
  Transformation motion;
  motion.set_reference_point(stk::math::Vector3d(0.5, -0.3, 0.2));
  motion.set_initial_displacement(stk::math::Vector3d(0.1, 0.2, -0.3));
  motion.set_initial_rotation(initialRotation);
  TransformedSurface moving(std::make_shared<Cuboid>(center, dimensions), motion);
  moving.prepare_to_compute(0., anyBox, 0.);
  const Cuboid placed(
      motion.rigid_transformation_at(0.).transform_point(center), dimensions, initialRotation);

  std::vector<double> placedPositions, movingPositions;
  size_t numCrossings = 0;
  for (auto & x0 : samplePoints)
    for (auto & x1 : samplePoints)
    {
      placed.fill_segment_crossings(x0, x1, 1.e-10, placedPositions);
      moving.fill_segment_crossings(x0, x1, 1.e-10, movingPositions);
      ASSERT_EQ(placedPositions.size(), movingPositions.size()) << x0 << " to " << x1;
      for (size_t i = 0; i < placedPositions.size(); ++i)
      {
        EXPECT_NEAR(placedPositions[i], movingPositions[i], 1.e-12) << x0 << " to " << x1;
      }
      numCrossings += movingPositions.size();
    }
  EXPECT_LT(0u, numCrossings);
}

// A moving box of facets is the box of facets with its vertices where the motion puts them.
TEST(TransformedSurface, facetedBox_sameDistancesAsFacetsWhereTheMotionPutsThem)
{
  const std::vector<stk::math::Vector3d> vertices{
      {0, 0, 0}, {1, 0, 0}, {1, 1, 0}, {0, 1, 0}, {0, 0, 1}, {1, 0, 1}, {1, 1, 1}, {0, 1, 1}};
  const std::vector<std::array<unsigned, 3>> allFacets{{{0, 2, 1}}, {{0, 3, 2}}, {{4, 5, 6}}, {{4, 6, 7}},
      {{0, 1, 5}}, {{0, 5, 4}}, {{1, 2, 6}}, {{1, 6, 5}}, {{2, 3, 7}}, {{2, 7, 6}}, {{3, 0, 4}}, {{3, 4, 7}}};
  const int numProcs = stk::parallel_machine_size(MPI_COMM_WORLD);
  const int me = stk::parallel_machine_rank(MPI_COMM_WORLD);
  std::vector<std::array<unsigned, 3>> myFacets;
  for (size_t i = 0; i < allFacets.size(); ++i)
    if (int(i % numProcs) == me) myFacets.push_back(allFacets[i]);

  const Transformation motion = general_motion();
  const double time = 0.7;
  const RigidTransformation placement = motion.rigid_transformation_at(time);
  std::vector<stk::math::Vector3d> placedVertices;
  for (auto & x : vertices)
    placedVertices.push_back(placement.transform_point(x));

  TransformedSurface moving(std::make_shared<MeshFacetedSurface>(MPI_COMM_WORLD, 3, vertices, myFacets), motion);
  MeshFacetedSurface placed(MPI_COMM_WORLD, 3, placedVertices, myFacets);
  moving.prepare_to_compute(time, anyBox, 0.);
  placed.prepare_to_compute(time, anyBox, 0.);
  for (auto & x : samplePoints)
  {
    EXPECT_NEAR(placed.point_signed_distance(x), moving.point_signed_distance(x), 1.e-10) << x;
  }
}

} // namespace krino
