// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include <gtest/gtest.h>

#include <Akri_Cutting_Surface.hpp>

#include <Akri_LevelSet.hpp>
#include <Akri_Plane.hpp>
#include <Akri_Plane_Intersections.hpp>
#include <Akri_Segment.hpp>
#include <Akri_UnitTestUtils.hpp>

#include <stk_util/parallel/Parallel.hpp>
#include <random>

namespace krino
{

namespace {
  int num_random_cases() { return 10000; }
  double clip(const double in) { return std::floor(1000*in)/1000.; }
} // namespace

TEST(Plane_Cutting_Surface, random_edge_cuts)
{
  const bool debug_output = false;
  const Vector3d plane_pt0(0., 0., 0.);
  const Vector3d plane_pt1(1., 0., 0.);
  const Vector3d plane_pt2(0., 1., 0.);
  Plane_Cutting_Surface surf(plane_pt0, plane_pt1, plane_pt2);


  stk::ParallelMachine pm = MPI_COMM_WORLD;
  const int parallel_rank = stk::parallel_machine_rank(pm);
  std::mt19937 mt(std::mt19937::default_seed + parallel_rank);
  std::uniform_real_distribution<double> coord(-1., 1.);

  const int num_cases = num_random_cases();
  for (int icase=0; icase<num_cases; ++icase)
  {
    const Vector3d pt0(clip(coord(mt)), clip(coord(mt)), clip(coord(mt)));
    const Vector3d pt1(clip(coord(mt)), clip(coord(mt)), clip(coord(mt)));
    const Segment3d segment(pt0, pt1);

    if (debug_output) std::cout << "Case: " << icase << std::endl;
    if (debug_output) std::cout << " Pt0 " << pt0 << std::endl;
    if (debug_output) std::cout << " Pt1 " << pt1 << std::endl;

    const double sign0 = surf.sign_at_position(pt0);
    const double sign1 = surf.sign_at_position(pt1);

    if (debug_output) std::cout << " sign " << sign0  << " " << sign1 << std::endl;

    if (sign0 != sign1)
    {
      const double position = surf.interface_crossing_position(segment);
      EXPECT_NEAR(position, (pt0[2] / (pt0[2] - pt1[2])), std::numeric_limits<double>::epsilon());
    }
    else
    {
      EXPECT_ANY_THROW(surf.interface_crossing_position(segment));
    }
  }
}

TEST(Intersecting_Planes_Cutting_Surface, planar_with_random_edge_cuts)
{
  const Vector3d plane_pt0(0., 0., 0.);
  const Vector3d plane_pt1(1., 0., 0.);
  const Vector3d plane_pt2(1., 1., 0.);
  const Vector3d plane_pt3(0., 1., 0.);
  Intersecting_Planes_Cutting_Surface surf(plane_pt0, plane_pt1, plane_pt2, plane_pt3);


  stk::ParallelMachine pm = MPI_COMM_WORLD;
  const int parallel_rank = stk::parallel_machine_rank(pm);
  std::mt19937 mt(std::mt19937::default_seed + parallel_rank);
  std::uniform_real_distribution<double> coord(-1., 1.);

  const int num_cases = num_random_cases();
  for (int icase=0; icase<num_cases; ++icase)
  {
    const Vector3d pt0(clip(coord(mt)), clip(coord(mt)), clip(coord(mt)));
    const Vector3d pt1(clip(coord(mt)), clip(coord(mt)), clip(coord(mt)));
    const Segment3d segment(pt0, pt1);

    const double sign0 = surf.sign_at_position(pt0);
    const double sign1 = surf.sign_at_position(pt1);

    if (sign0 != sign1)
    {
      const double position = surf.interface_crossing_position(segment);
      EXPECT_NEAR(position, (pt0[2] / (pt0[2] - pt1[2])), std::numeric_limits<double>::epsilon());
    }
    else
    {
      EXPECT_ANY_THROW(surf.interface_crossing_position(segment));
    }
  }
}

TEST(Intersecting_Planes_Cutting_Surface, positive_dihedral_with_random_edge_cuts)
{
  const bool debug_output = false;
  const Vector3d plane_pt0(0., 0., 0.);
  const Vector3d plane_pt1(0., 1., 1.);
  const Vector3d plane_pt2(1., 0., 0.);
  const Vector3d plane_pt3(0., 1., 0.);
  Intersecting_Planes_Cutting_Surface surf(plane_pt0, plane_pt1, plane_pt2, plane_pt3);


  stk::ParallelMachine pm = MPI_COMM_WORLD;
  const int parallel_rank = stk::parallel_machine_rank(pm);
  std::mt19937 mt(std::mt19937::default_seed + parallel_rank);
  std::uniform_real_distribution<double> coord(-1., 1.);

  const int num_cases = num_random_cases();
  for (int icase=0; icase<num_cases; ++icase)
  {
    const Vector3d pt0(clip(coord(mt)), clip(coord(mt)), clip(coord(mt)));
    const Vector3d pt1(clip(coord(mt)), clip(coord(mt)), clip(coord(mt)));
    const Segment3d segment(pt0, pt1);

    if (debug_output) std::cout << "Case: " << icase << std::endl;
    if (debug_output) std::cout << " Pt0 " << pt0 << std::endl;
    if (debug_output) std::cout << " Pt1 " << pt1 << std::endl;

    const double sign0 = surf.sign_at_position(pt0);
    const double sign1 = surf.sign_at_position(pt1);

    if (debug_output) std::cout << " sign " << sign0  << " " << sign1 << std::endl;

    if (sign0 != sign1)
    {
      EXPECT_FALSE((pt0[1] < 0.) && (pt1[1] < 0.));
      EXPECT_FALSE((pt0[2] < 0.) && (pt1[2] < 0.));
      const double position = surf.interface_crossing_position(segment);
      const double pos0 = (pt0[2]-pt0[1]) / ((pt0[2]-pt0[1]) - (pt1[2]-pt1[1]));
      const double pos1 = pt0[2] / (pt0[2] - pt1[2]);

      // Only checking cases where both points are on same side of dividing plane
      const int case_id = ((pt0[2] > 0.5*pt0[1]) ? 1 : 0) + ((pt1[2] > 0.5*pt1[1]) ? 2 : 0);
      switch (case_id)
      {
        case 0: EXPECT_NEAR(position, pos1, std::numeric_limits<float>::epsilon()); break;
        case 3: EXPECT_NEAR(position, pos0, std::numeric_limits<float>::epsilon()); break;
      }
    }
    else
    {
      EXPECT_ANY_THROW(surf.interface_crossing_position(segment));
    }
  }
}

TEST(Intersecting_Planes_Cutting_Surface, negative_dihedral_with_random_edge_cuts)
{
  const bool debug_output = false;
  const Vector3d plane_pt0(0., 0., 0.);
  const Vector3d plane_pt1(0., 1., 0.);
  const Vector3d plane_pt2(1., 0., 0.);
  const Vector3d plane_pt3(0., 1., 1.);
  Intersecting_Planes_Cutting_Surface surf(plane_pt0, plane_pt1, plane_pt2, plane_pt3);


  stk::ParallelMachine pm = MPI_COMM_WORLD;
  const int parallel_rank = stk::parallel_machine_rank(pm);
  std::mt19937 mt(std::mt19937::default_seed + parallel_rank);
  std::uniform_real_distribution<double> coord(-1., 1.);

  const int num_cases = num_random_cases();
  for (int icase=0; icase<num_cases; ++icase)
  {
    const Vector3d pt0(clip(coord(mt)), clip(coord(mt)), clip(coord(mt)));
    const Vector3d pt1(clip(coord(mt)), clip(coord(mt)), coord(mt));
    const Segment3d segment(pt0, pt1);

    if (debug_output) std::cout << "Case: " << icase << std::endl;
    if (debug_output) std::cout << " Pt0 " << pt0 << std::endl;
    if (debug_output) std::cout << " Pt1 " << pt1 << std::endl;

    const double sign0 = surf.sign_at_position(pt0);
    const double sign1 = surf.sign_at_position(pt1);

    if (debug_output) std::cout << " sign " << sign0  << " " << sign1 << std::endl;

    if (sign0 != sign1)
    {
      EXPECT_FALSE((pt0[1] <= 0.) && (pt1[1] <= 0.));
      EXPECT_FALSE((pt0[2] <= 0.) && (pt1[2] <= 0.));
      const double position = surf.interface_crossing_position(segment);
      const double pos0 = (pt0[2]-pt0[1]) / ((pt0[2]-pt0[1]) - (pt1[2]-pt1[1]));
      const double pos1 = pt0[2] / (pt0[2] - pt1[2]);
      // Only checking cases where both points are on same side of dividing plane
      const int case_id = ((pt0[2] > 0.5*pt0[1]) ? 1 : 0) + ((pt1[2] > 0.5*pt1[1]) ? 2 : 0);
      switch (case_id)
      {
        case 0: EXPECT_NEAR(position, pos1, std::numeric_limits<float>::epsilon()); break;
        case 3: EXPECT_NEAR(position, pos0, std::numeric_limits<float>::epsilon()); break;
      }
    }
    else
    {
      EXPECT_ANY_THROW(surf.interface_crossing_position(segment));
    }
  }
}

TEST(Intersecting_Planes_Cutting_Surface, infinitesimal_triangle_that_requires_robust_dihedral_angle)
{
  const double goldPosition = 0.5181038869168293;
  const double otherCrossing = 0.4818961330726974;
  const Vector3d plane_pt0(otherCrossing, 0., 1.-otherCrossing);
  const Vector3d plane_pt1(0., 0., goldPosition);
  const Vector3d plane_pt2(0., (1.-1.e-10), 0.);
  const Vector3d plane_pt3(1.e-10, (1.-1.e-10), 0.);
  Intersecting_Planes_Cutting_Surface surf(plane_pt0, plane_pt1, plane_pt2, plane_pt3);

  const Vector3d node0(0., 0., 0.);
  const Vector3d node3(0., 0., 1.);
  const double position = surf.interface_crossing_position(Segment3d(node0,node3));
  EXPECT_DOUBLE_EQ(position, goldPosition);
}

void expect_intersection(const Plane3d & plane0, const Plane3d & plane1, const Plane3d & plane2, const Vector3d & goldIntersection)
{
  Vector3d intersectionPoint;
  EXPECT_TRUE(find_intersection_of_three_planes(plane0, plane1, plane2, intersectionPoint));
  expect_eq(goldIntersection, intersectionPoint);
}

void expect_intersection_with_side_of_tet(const int side, const Plane3d & plane0, const Plane3d & plane1, const Vector3d & goldIntersection)
{
  Vector3d intersectionPoint;
  EXPECT_TRUE(find_intersection_of_two_planes_and_side_of_tet(side, plane0, plane1, intersectionPoint));
  expect_eq(goldIntersection, intersectionPoint);
}

void expect_no_intersection(const Plane3d & plane0, const Plane3d & plane1, const Plane3d & plane2)
{
  Vector3d intersectionPoint;
  EXPECT_FALSE(find_intersection_of_three_planes(plane0, plane1, plane2, intersectionPoint));
}

void expect_no_intersection_within_tet(const Plane3d & plane0, const Plane3d & plane1, const Plane3d & plane2)
{
  Vector3d intersectionPoint;
  EXPECT_FALSE(find_intersection_of_three_planes_within_tet(plane0, plane1, plane2, intersectionPoint));
}

void expect_no_intersection_with_side_of_tet(const int side, const Plane3d & plane0, const Plane3d & plane1)
{
  Vector3d intersectionPoint;
  EXPECT_FALSE(find_intersection_of_two_planes_and_side_of_tet(side, plane0, plane1, intersectionPoint));
}

struct IntersectPlanes : public ::testing::Test
{
  const Vector3d pt0{0., 0., 0.};
  const Vector3d pt1{1., 0., 0.};
  const Vector3d pt2{0., 1., 0.};
  const Vector3d pt3{0., 0., 1.};
  const Plane3d plane0{pt0, pt1, pt3};
  const Plane3d plane1{pt1, pt2, pt3};
  const Plane3d plane2{pt2, pt0, pt3};
  const Plane3d plane3{pt0, pt2, pt1};

  const Vector3d offsetPt0{-0.1, -0.1, -0.1};
  const Vector3d offsetPt1{1.1, -0.1, -0.1};
  const Vector3d offsetPt2{-0.1, 1.1, -0.1};
  const Vector3d offsetPt3{-0.1, -0.1, 1.1};
  const Plane3d offsetPlane0{offsetPt0, offsetPt1, offsetPt3};
  const Plane3d offsetPlane1{offsetPt1, offsetPt2, offsetPt3};
  const Plane3d offsetPlane2{offsetPt2, offsetPt0, offsetPt3};
  const Plane3d offsetPlane3{offsetPt0, offsetPt2, offsetPt1};

};

TEST_F(IntersectPlanes, given3Planes_FindCorrectIntersection)
{
  expect_intersection(plane0, plane1, plane2, Vector3d{0., 0., 1.});
  expect_intersection(plane0, plane1, plane3, Vector3d{1., 0., 0.});
  expect_intersection(plane1, plane2, plane3, Vector3d{0., 1., 0.});
  expect_intersection(plane0, plane2, plane3, Vector3d{0., 0., 0.});

  expect_no_intersection(plane0, plane0, plane1);
}

TEST_F(IntersectPlanes, given3PlanesThatDoNotIntersectWithinTetBounds_FindNoIntersection)
{
  expect_no_intersection_within_tet(offsetPlane0, offsetPlane1, offsetPlane2);
  expect_no_intersection_within_tet(offsetPlane0, offsetPlane1, offsetPlane3);
  expect_no_intersection_within_tet(offsetPlane1, offsetPlane2, offsetPlane3);
  expect_no_intersection_within_tet(offsetPlane0, offsetPlane2, offsetPlane3);
}

TEST_F(IntersectPlanes, given2PlanesThatIntersectSideOfTet_FindCorrectIntersection)
{
  expect_intersection_with_side_of_tet(0, plane1, plane2, Vector3d{0., 0., 1.});
  expect_intersection_with_side_of_tet(0, plane1, plane3, Vector3d{1., 0., 0.});
  expect_intersection_with_side_of_tet(0, plane2, plane3, Vector3d{0., 0., 0.});

  expect_intersection_with_side_of_tet(1, plane2, plane3, Vector3d{0., 1., 0.});
  expect_intersection_with_side_of_tet(1, plane0, plane2, Vector3d{0., 0., 1.});
  expect_intersection_with_side_of_tet(1, plane0, plane3, Vector3d{1., 0., 0.});

  expect_intersection_with_side_of_tet(2, plane0, plane1, Vector3d{0., 0., 1.});
  expect_intersection_with_side_of_tet(2, plane0, plane3, Vector3d{0., 0., 0.});
  expect_intersection_with_side_of_tet(2, plane1, plane3, Vector3d{0., 1., 0.});

  expect_intersection_with_side_of_tet(3, plane0, plane1, Vector3d{1., 0., 0.});
  expect_intersection_with_side_of_tet(3, plane0, plane2, Vector3d{0., 0., 0.});
  expect_intersection_with_side_of_tet(3, plane1, plane2, Vector3d{0., 1., 0.});
}

TEST_F(IntersectPlanes, given2PlanesThatDoNotIntersectSideOfTetAtAll_FindNoIntersection)
{
  expect_no_intersection_with_side_of_tet(0, plane0, plane2);
  expect_no_intersection_with_side_of_tet(1, plane1, plane2);
  expect_no_intersection_with_side_of_tet(2, plane2, plane1);
  expect_no_intersection_with_side_of_tet(3, plane3, plane2);
}

TEST_F(IntersectPlanes, given2PlanesThatDoNotIntersectSideOfTetWithinTetBounds_FindNoIntersection)
{
  expect_no_intersection_with_side_of_tet(0, offsetPlane1, offsetPlane2);
  expect_no_intersection_with_side_of_tet(0, offsetPlane1, offsetPlane3);
  expect_no_intersection_with_side_of_tet(0, offsetPlane2, offsetPlane3);

  expect_no_intersection_with_side_of_tet(1, offsetPlane2, offsetPlane3);
  expect_no_intersection_with_side_of_tet(1, offsetPlane0, offsetPlane2);
  expect_no_intersection_with_side_of_tet(1, offsetPlane0, offsetPlane3);

  expect_no_intersection_with_side_of_tet(2, offsetPlane0, offsetPlane1);
  expect_no_intersection_with_side_of_tet(2, offsetPlane0, offsetPlane3);
  expect_no_intersection_with_side_of_tet(2, offsetPlane1, offsetPlane3);

  expect_no_intersection_with_side_of_tet(3, offsetPlane0, offsetPlane1);
  expect_no_intersection_with_side_of_tet(3, offsetPlane0, offsetPlane2);
  expect_no_intersection_with_side_of_tet(3, offsetPlane1, offsetPlane2);
}

static Vector3d compute_2d_plane_direction(const Vector3d & pt0, const Vector3d & pt1)
{
  return Vector3d(pt1[1]-pt0[1],pt0[0]-pt1[0],0.);
}

static Plane3d build_2d_plane(const Vector3d & pt0, const Vector3d & pt1)
{
  return Plane3d(compute_2d_plane_direction(pt0,pt1),pt0);
}

struct Intersect2DPlanes : public ::testing::Test
{
  Intersect2DPlanes()
  : plane0{build_2d_plane(pt0, pt1)},
    plane1{build_2d_plane(pt1, pt2)},
    plane2{build_2d_plane(pt0, pt2)},
    offsetPlane0{build_2d_plane(offsetPt0, offsetPt1)},
    offsetPlane1{build_2d_plane(offsetPt1, offsetPt2)},
    offsetPlane2{build_2d_plane(offsetPt0, offsetPt2)}
  {
  }

  const Vector3d pt0{0., 0., 0.};
  const Vector3d pt1{1., 0., 0.};
  const Vector3d pt2{0., 1., 0.};

  const Vector3d offsetPt0{-0.1, -0.1, 0.};
  const Vector3d offsetPt1{1.1, -0.1, 0.};
  const Vector3d offsetPt2{-0.1, 1.1, 0.};

  Plane3d plane0;
  Plane3d plane1;
  Plane3d plane2;

  Plane3d offsetPlane0;
  Plane3d offsetPlane1;
  Plane3d offsetPlane2;
};

void expect_2d_intersection(const Plane3d & plane0, const Plane3d & plane1, const Vector3d & goldIntersection)
{
  Vector3d intersectionPoint;
  EXPECT_TRUE(find_intersection_of_two_2D_planes(plane0, plane1, intersectionPoint));
  expect_eq(goldIntersection, intersectionPoint);
}

void expect_no_2d_intersection_within_tri(const Plane3d & plane0, const Plane3d & plane1)
{
  Vector3d intersectionPoint;
  EXPECT_FALSE(find_intersection_of_two_2D_planes_within_tri(plane0, plane1, intersectionPoint));
}

void expect_no_2d_intersection(const Plane3d & plane0, const Plane3d & plane1)
{
  Vector3d intersectionPoint;
  EXPECT_FALSE(find_intersection_of_two_2D_planes(plane0, plane1, intersectionPoint));
}

TEST_F(Intersect2DPlanes, given2Plane2Ds_FindCorrectIntersection)
{
  expect_2d_intersection(plane0, plane1, Vector3d{1., 0., 0.});
  expect_2d_intersection(plane1, plane2, Vector3d{0., 1., 0.});
  expect_2d_intersection(plane0, plane2, Vector3d{0., 0., 0.});

  expect_no_2d_intersection(plane0, plane0);
}

TEST_F(Intersect2DPlanes, given2Plane2DsThatDoNotIntersectWithinTriBounds_FindNoIntersection)
{
  expect_no_2d_intersection_within_tri(offsetPlane0, offsetPlane1);
  expect_no_2d_intersection_within_tri(offsetPlane1, offsetPlane2);
  expect_no_2d_intersection_within_tri(offsetPlane0, offsetPlane2);
}

TEST(ProjectionOf3DPointsIntoTriangle, givenTriangleCheckProjectionOfPointOnAndOffPlane)
{
  const Vector3d pt0{0., 0., 0.};
  const Vector3d pt1{1., 0., 0.};
  const Vector3d pt2{0., 1., 1.};
  const std::array<Vector3d,3> triCoords{{ pt0, pt1, pt2 }};

  expect_eq(Vector3d{0.,0.,0.}, triangle_parametric_coordinates_of_projected_point(triCoords, pt0));
  expect_eq(Vector3d{1.,0.,0.}, triangle_parametric_coordinates_of_projected_point(triCoords, pt1));
  expect_eq(Vector3d{0.,1.,0.}, triangle_parametric_coordinates_of_projected_point(triCoords, pt2));

  const Vector3d normal = Cross(pt1-pt0, pt2-pt0).unit_vector();
  expect_eq(Vector3d{0.,0.,0.}, triangle_parametric_coordinates_of_projected_point(triCoords, pt0+normal));
  expect_eq(Vector3d{1.,0.,0.}, triangle_parametric_coordinates_of_projected_point(triCoords, pt1+normal));
  expect_eq(Vector3d{0.,1.,0.}, triangle_parametric_coordinates_of_projected_point(triCoords, pt2+normal));

  const double oneThird = 1./3.;
  const Vector3d midPt = oneThird*(pt0+pt1+pt2);
  expect_eq(Vector3d{oneThird,oneThird,0.}, triangle_parametric_coordinates_of_projected_point(triCoords, midPt+normal));
}

} // namespace krino
