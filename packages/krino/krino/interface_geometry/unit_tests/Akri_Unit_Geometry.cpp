// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include <gtest/gtest.h>

#include <Akri_Cutting_Surface.hpp>

#include <stk_math/StkVector.hpp>
#include <stk_util/parallel/Parallel.hpp>

#include <array>
#include <cmath>
#include <limits>
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
  const stk::math::Vector3d plane_pt0(0., 0., 0.);
  const stk::math::Vector3d plane_pt1(1., 0., 0.);
  const stk::math::Vector3d plane_pt2(0., 1., 0.);
  Plane_Cutting_Surface surf(plane_pt0, plane_pt1, plane_pt2);


  stk::ParallelMachine pm = MPI_COMM_WORLD;
  const int parallel_rank = stk::parallel_machine_rank(pm);
  std::mt19937 mt(std::mt19937::default_seed + parallel_rank);
  std::uniform_real_distribution<double> coord(-1., 1.);

  const int num_cases = num_random_cases();
  for (int icase=0; icase<num_cases; ++icase)
  {
    const stk::math::Vector3d pt0(clip(coord(mt)), clip(coord(mt)), clip(coord(mt)));
    const stk::math::Vector3d pt1(clip(coord(mt)), clip(coord(mt)), clip(coord(mt)));
    const std::array<stk::math::Vector3d,2> segment{pt0, pt1};

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
  const stk::math::Vector3d plane_pt0(0., 0., 0.);
  const stk::math::Vector3d plane_pt1(1., 0., 0.);
  const stk::math::Vector3d plane_pt2(1., 1., 0.);
  const stk::math::Vector3d plane_pt3(0., 1., 0.);
  Intersecting_Planes_Cutting_Surface surf(plane_pt0, plane_pt1, plane_pt2, plane_pt3);


  stk::ParallelMachine pm = MPI_COMM_WORLD;
  const int parallel_rank = stk::parallel_machine_rank(pm);
  std::mt19937 mt(std::mt19937::default_seed + parallel_rank);
  std::uniform_real_distribution<double> coord(-1., 1.);

  const int num_cases = num_random_cases();
  for (int icase=0; icase<num_cases; ++icase)
  {
    const stk::math::Vector3d pt0(clip(coord(mt)), clip(coord(mt)), clip(coord(mt)));
    const stk::math::Vector3d pt1(clip(coord(mt)), clip(coord(mt)), clip(coord(mt)));
    const std::array<stk::math::Vector3d,2> segment{pt0, pt1};

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
  const stk::math::Vector3d plane_pt0(0., 0., 0.);
  const stk::math::Vector3d plane_pt1(0., 1., 1.);
  const stk::math::Vector3d plane_pt2(1., 0., 0.);
  const stk::math::Vector3d plane_pt3(0., 1., 0.);
  Intersecting_Planes_Cutting_Surface surf(plane_pt0, plane_pt1, plane_pt2, plane_pt3);


  stk::ParallelMachine pm = MPI_COMM_WORLD;
  const int parallel_rank = stk::parallel_machine_rank(pm);
  std::mt19937 mt(std::mt19937::default_seed + parallel_rank);
  std::uniform_real_distribution<double> coord(-1., 1.);

  const int num_cases = num_random_cases();
  for (int icase=0; icase<num_cases; ++icase)
  {
    const stk::math::Vector3d pt0(clip(coord(mt)), clip(coord(mt)), clip(coord(mt)));
    const stk::math::Vector3d pt1(clip(coord(mt)), clip(coord(mt)), clip(coord(mt)));
    const std::array<stk::math::Vector3d,2> segment{pt0, pt1};

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
  const stk::math::Vector3d plane_pt0(0., 0., 0.);
  const stk::math::Vector3d plane_pt1(0., 1., 0.);
  const stk::math::Vector3d plane_pt2(1., 0., 0.);
  const stk::math::Vector3d plane_pt3(0., 1., 1.);
  Intersecting_Planes_Cutting_Surface surf(plane_pt0, plane_pt1, plane_pt2, plane_pt3);


  stk::ParallelMachine pm = MPI_COMM_WORLD;
  const int parallel_rank = stk::parallel_machine_rank(pm);
  std::mt19937 mt(std::mt19937::default_seed + parallel_rank);
  std::uniform_real_distribution<double> coord(-1., 1.);

  const int num_cases = num_random_cases();
  for (int icase=0; icase<num_cases; ++icase)
  {
    const stk::math::Vector3d pt0(clip(coord(mt)), clip(coord(mt)), clip(coord(mt)));
    const stk::math::Vector3d pt1(clip(coord(mt)), clip(coord(mt)), coord(mt));
    const std::array<stk::math::Vector3d,2> edgeNodeCoords{pt0, pt1};

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
      const double position = surf.interface_crossing_position(edgeNodeCoords);
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
      EXPECT_ANY_THROW(surf.interface_crossing_position(edgeNodeCoords));
    }
  }
}

TEST(Intersecting_Planes_Cutting_Surface, infinitesimal_triangle_that_requires_robust_dihedral_angle)
{
  const double goldPosition = 0.5181038869168293;
  const double otherCrossing = 0.4818961330726974;
  const stk::math::Vector3d plane_pt0(otherCrossing, 0., 1.-otherCrossing);
  const stk::math::Vector3d plane_pt1(0., 0., goldPosition);
  const stk::math::Vector3d plane_pt2(0., (1.-1.e-10), 0.);
  const stk::math::Vector3d plane_pt3(1.e-10, (1.-1.e-10), 0.);
  Intersecting_Planes_Cutting_Surface surf(plane_pt0, plane_pt1, plane_pt2, plane_pt3);

  const stk::math::Vector3d node0(0., 0., 0.);
  const stk::math::Vector3d node3(0., 0., 1.);
  const std::array<stk::math::Vector3d,2> segment{node0,node3};
  const double position = surf.interface_crossing_position(segment);
  EXPECT_DOUBLE_EQ(position, goldPosition);
}

}
