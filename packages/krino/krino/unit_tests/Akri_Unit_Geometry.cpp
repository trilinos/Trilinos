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
#include <stk_math/StkPlane.hpp>
#include <Akri_Plane_Intersections.hpp>
#include <Akri_Segment.hpp>
#include <Akri_SegmentWithSensitivities.hpp>
#include <Akri_Triangle.hpp>
#include <Akri_TriangleWithSensitivities.hpp>
#include <Akri_UnitTestUtils.hpp>

#include <stk_util/parallel/Parallel.hpp>
#include <random>

namespace krino
{

namespace {
  int num_random_cases() { return 10000; }
  double clip(const double in) { return std::floor(1000*in)/1000.; }
} // namespace

namespace analytic_detail
{
  // Provided by Sean Hardesty

  void make_crmat(const stk::math::Vector3d &x, double *M) {
      M[1] += x[2];
      M[2] -= x[1];
      M[3] -= x[2];
      M[5] += x[0];
      M[6] += x[1];
      M[7] -= x[0];
  }
  void darea_mat(const stk::math::Vector3d v0, const stk::math::Vector3d v1, const stk::math::Vector3d v2, double *A) {
      std::fill_n(A, 3*9, 0.0);
      make_crmat(v2-v1, A);
      make_crmat(v0-v2, A+9);
      make_crmat(v1-v0, A+18);
  }

  double two_area_and_gradient(const stk::math::Vector3d v0, const stk::math::Vector3d v1, const stk::math::Vector3d v2, double *dtwoArea)
  {
    stk::math::Vector3d n = CalcTriangle3<double>::two_times_area_vector(v0,v1,v2);
    const double twoArea = n.length();
    n /= twoArea;
    double A[3*9];
    darea_mat(v0, v1, v2, A);

    //blas::gemhv(3, 9, 1.0, A, n.getdata(), dtwoArea);
    std::fill_n(dtwoArea,9,0.0);
    for(unsigned i=0;i<3;i++) {
        for(unsigned j=0;j<9;j++) {
            dtwoArea[j] += A[3*j+i]*n[i];
        }
    }

    return twoArea;
  }

  stk::math::Vector3d normal_and_gradient(const stk::math::Vector3d v0, const stk::math::Vector3d v1, const stk::math::Vector3d v2, double *dnormal)
  {
    stk::math::Vector3d n = CalcTriangle3<double>::two_times_area_vector(v0,v1,v2);
    const double twoArea = n.length();
    n /= twoArea;
    double A[3*9];
    darea_mat(v0, v1 ,v2, A);
    double B[3*3] = {1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0};

    //blas::ger(3, 3, -1.0, n.getdata(), 1, n.getdata(), 1, B, 3);
    for(unsigned i=0;i<3;i++) {
        for(unsigned j=0;j<3;j++) {
            B[3*j+i] -= n[i]*n[j];
        }
    }

    //blas::gemm(3, 3, 9, 1.0/twoArea, B, 3, A, 3, dnormal, 3);
    const double twoArea_ = 1.0/twoArea;
    std::fill_n(dnormal, 3*9, 0.0);
    for(unsigned i=0;i<3;i++) {
        for(unsigned j=0;j<9;j++) {
            for(unsigned k=0;k<3;k++) {
                dnormal[3*j+i] += B[3*k+i]*A[3*j+k]*twoArea_;
            }
        }
    }

    return n;
  }
}

namespace analytic_detail
{
stk::math::Vector3d difference_and_sensitivity(const stk::math::Vector3d v0, const stk::math::Vector3d v1, double *dDiff)
{
  for(unsigned n=0;n<2;n++)
    for(unsigned i=0;i<3;i++)
      for(unsigned j=0;j<3;j++)
        dDiff[6*i + (3*n+j)] = ((n==0)?(-1.):(1.)) * ((i==j)?(1.):(0.));

  return v1-v0;
}

double length_and_gradient(const stk::math::Vector3d v0, const stk::math::Vector3d v1, double *dLength)
{
  double dDiff[3*2*3];
  const stk::math::Vector3d diff = difference_and_sensitivity(v0, v1, dDiff);

  const double length = diff.length();

  std::fill_n(dLength, 3*2, 0.0);

  for(unsigned n=0;n<2;n++)
    for(unsigned i=0;i<3;i++)
      for(unsigned j=0;j<3;j++)
        dLength[3*n + j] += diff[i]*dDiff[6*i + (3*n+j)] / length;

  return length;
}

stk::math::Vector3d normal2d_and_gradient(const stk::math::Vector3d v0, const stk::math::Vector3d v1, double *dNormal)
{
  // use 3d calc
  const stk::math::Vector3d p0(v0[0], v0[1], 0.);
  const stk::math::Vector3d p1(v1[0], v1[1], 0.);
  const stk::math::Vector3d p2(v0[0], v0[1], 1.);

  double dNormalAnalytic[27];
  const stk::math::Vector3d normalAnalytic = analytic_detail::normal_and_gradient(p0, p1, p2, dNormalAnalytic);

  for(unsigned i=0;i<3;i++)
    for(unsigned j=0;j<6;j++)
      dNormal[j*3+i] = ((i==2) ? 0. : dNormalAnalytic[j*3+i]);

  return normalAnalytic;
}
}

TEST(Triangle_Calcs, Cartesian_normal_and_area)
{
  const double tol = 1.e-9;
  const stk::math::Vector3d pt0{0., 0., 0.};
  const stk::math::Vector3d ptx{1., 0., 0.};
  const stk::math::Vector3d pty{0., 1., 0.};
  const stk::math::Vector3d ptz{0., 0., 1.};

  expect_eq(ptz, CalcTriangle3<double>::normal(pt0, ptx, pty));
  expect_eq(ptx, CalcTriangle3<double>::normal(pt0, pty, ptz));
  expect_eq(pty, CalcTriangle3<double>::normal(pt0, ptz, ptx));
  expect_eq(-ptz, CalcTriangle3<double>::normal(pt0, pty, ptx));
  expect_eq(-ptx, CalcTriangle3<double>::normal(pt0, ptz, pty));
  expect_eq(-pty, CalcTriangle3<double>::normal(pt0, ptx, ptz));

  EXPECT_NEAR(0.5, CalcTriangle3<double>::area(pt0, ptx, pty), tol);
  EXPECT_NEAR(0.5, CalcTriangle3<double>::area(pt0, pty, ptz), tol);
  EXPECT_NEAR(0.5, CalcTriangle3<double>::area(pt0, ptz, ptx), tol);
}

void expect_matching_area_and_normal_values_and_sensitivities(const stk::math::Vector3d v0, const stk::math::Vector3d v1, const stk::math::Vector3d v2)
{
  const double tol = 1.e-9;

  double dNormal[27];
  const stk::math::Vector3d normal = TriangleWithSens::normal_and_optional_sensitivities(v0, v1, v2, dNormal);

  double dArea[9];
  const double area = TriangleWithSens::area_and_optional_sensitivities(v0, v1, v2, dArea);

  double dAreaVector[27];
  const stk::math::Vector3d areaVector = TriangleWithSens::area_vector_and_optional_sensitivities(v0, v1, v2, dAreaVector);

  double dNormalAnalytic[27];
  const stk::math::Vector3d normalAnalytic = analytic_detail::normal_and_gradient(v0, v1, v2, dNormalAnalytic);

  double dtwoAreaAnalytic[9];
  const double twoAreaAnalytic = analytic_detail::two_area_and_gradient(v0, v1, v2, dtwoAreaAnalytic);

  for (unsigned i=0; i<3; ++i)
  {
    EXPECT_NEAR(normalAnalytic[i], normal[i], tol);
    EXPECT_NEAR(0.5*twoAreaAnalytic*normalAnalytic[i], areaVector[i], tol);

    for (unsigned j=0; j<9; ++j)
    {
      EXPECT_NEAR(dNormalAnalytic[j*3+i], dNormal[j*3+i], tol);
      EXPECT_NEAR(0.5*twoAreaAnalytic*dNormalAnalytic[j*3+i] + 0.5*normalAnalytic[i]*dtwoAreaAnalytic[j], dAreaVector[j*3+i], tol);
    }
  }

  EXPECT_NEAR(0.5*twoAreaAnalytic, area, tol);
  for (unsigned j=0; j<9; ++j)
  {
    EXPECT_NEAR(0.5*dtwoAreaAnalytic[j], dArea[j], tol);
  }
}

TEST(Triangle_Calcs, normal_and_area_sens_check)
{
  stk::ParallelMachine pm = MPI_COMM_WORLD;
  const int parallel_rank = stk::parallel_machine_rank(pm);
  std::mt19937 mt(std::mt19937::default_seed + parallel_rank);
  std::uniform_real_distribution<double> coord(-1., 1.);

  const int num_cases = 100;
  for (int icase=0; icase<num_cases; ++icase)
  {
    const stk::math::Vector3d pt0(coord(mt), coord(mt), coord(mt));
    const stk::math::Vector3d pt1(coord(mt), coord(mt), coord(mt));
    const stk::math::Vector3d pt2(coord(mt), coord(mt), coord(mt));

    expect_matching_area_and_normal_values_and_sensitivities(pt0, pt1, pt2);
  }
}

TEST(Triangle_Calcs, normal_and_area_performance_check)
{
  stk::ParallelMachine pm = MPI_COMM_WORLD;
  const int parallel_rank = stk::parallel_machine_rank(pm);
  std::uniform_real_distribution<double> coord(-1., 1.);

  double dNormal[27];
  double dArea[9];

  const int numCases = 100000;
  stk::diag::Timer timer{"Calc Triangle With Sensitivities", sierra::Diag::sierraTimer()};

  {
    stk::diag::TimeBlock timer__(timer);
    std::mt19937 mt(std::mt19937::default_seed + parallel_rank);

    for (int icase=0; icase<numCases; ++icase)
    {
      const stk::math::Vector3d pt0(coord(mt), coord(mt), coord(mt));
      const stk::math::Vector3d pt1(coord(mt), coord(mt), coord(mt));
      const stk::math::Vector3d pt2(coord(mt), coord(mt), coord(mt));

      analytic_detail::normal_and_gradient(pt0, pt1, pt2, dNormal);
      analytic_detail::two_area_and_gradient(pt0, pt1, pt2, dArea);
    }
  }
  std::cout << "Using analytic sensitivities, " << numCases << " cases uses time = " << timer.getMetric<stk::diag::CPUTime>().getLap() << std::endl;

  {
    stk::diag::TimeBlock timer__(timer);
    std::mt19937 mt(std::mt19937::default_seed + parallel_rank);

    for (int icase=0; icase<numCases; ++icase)
    {
      const stk::math::Vector3d pt0(coord(mt), coord(mt), coord(mt));
      const stk::math::Vector3d pt1(coord(mt), coord(mt), coord(mt));
      const stk::math::Vector3d pt2(coord(mt), coord(mt), coord(mt));

      TriangleWithSens::normal_and_optional_sensitivities(pt0, pt1, pt2, dNormal);
      TriangleWithSens::area_and_optional_sensitivities(pt0, pt1, pt2, dArea);
    }
  }
  std::cout << "Using FAD, " << numCases << " cases uses time = " << timer.getMetric<stk::diag::CPUTime>().getLap() << std::endl;
}

void expect_matching_length_and_normal2d_values_and_sensitivities(const stk::math::Vector3d v0, const stk::math::Vector3d v1)
{
  const double tol = 1.e-9;

  double dNormal[18];
  const stk::math::Vector3d normal = SegmentWithSens::normal2d_and_optional_sensitivities(v0, v1, dNormal);

  double dLength[6];
  const double length = SegmentWithSens::length_and_optional_sensitivities(v0, v1, dLength);

  double dNormalAnalytic[18];
  const stk::math::Vector3d normalAnalytic = analytic_detail::normal2d_and_gradient(v0, v1, dNormalAnalytic);

  double dLengthAnalytic[9];
  const double lengthAnalytic = analytic_detail::length_and_gradient(v0, v1, dLengthAnalytic);

  for (unsigned i=0; i<3; ++i)
  {
    EXPECT_NEAR(normalAnalytic[i], normal[i], tol);

    for (unsigned j=0; j<6; ++j)
    {
      EXPECT_NEAR(dNormalAnalytic[j*3+i], dNormal[j*3+i], tol);
    }
  }

  EXPECT_NEAR(lengthAnalytic, length, tol);
  for (unsigned j=0; j<6; ++j)
  {
    EXPECT_NEAR(dLengthAnalytic[j], dLength[j], tol);
  }
}

TEST(Segment_Calcs, normal2d_and_length_sens_check)
{
  stk::ParallelMachine pm = MPI_COMM_WORLD;
  const int parallel_rank = stk::parallel_machine_rank(pm);
  std::mt19937 mt(std::mt19937::default_seed + parallel_rank);
  std::uniform_real_distribution<double> coord(-1., 1.);

  const int num_cases = 100;
  for (int icase=0; icase<num_cases; ++icase)
  {
    const stk::math::Vector3d pt0(coord(mt), coord(mt), 0.);
    const stk::math::Vector3d pt1(coord(mt), coord(mt), 0.);

    expect_matching_length_and_normal2d_values_and_sensitivities(pt0, pt1);
  }
}

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

void expect_intersection(const stk::math::Plane3d & plane0, const stk::math::Plane3d & plane1, const stk::math::Plane3d & plane2, const stk::math::Vector3d & goldIntersection)
{
  stk::math::Vector3d intersectionPoint;
  EXPECT_TRUE(find_intersection_of_three_planes(plane0, plane1, plane2, intersectionPoint));
  expect_eq(goldIntersection, intersectionPoint);
}

void expect_intersection_with_side_of_tet(const int side, const stk::math::Plane3d & plane0, const stk::math::Plane3d & plane1, const stk::math::Vector3d & goldIntersection)
{
  stk::math::Vector3d intersectionPoint;
  EXPECT_TRUE(find_intersection_of_two_planes_and_side_of_tet(side, plane0, plane1, intersectionPoint));
  expect_eq(goldIntersection, intersectionPoint);
}

void expect_no_intersection(const stk::math::Plane3d & plane0, const stk::math::Plane3d & plane1, const stk::math::Plane3d & plane2)
{
  stk::math::Vector3d intersectionPoint;
  EXPECT_FALSE(find_intersection_of_three_planes(plane0, plane1, plane2, intersectionPoint));
}

void expect_no_intersection_within_tet(const stk::math::Plane3d & plane0, const stk::math::Plane3d & plane1, const stk::math::Plane3d & plane2)
{
  stk::math::Vector3d intersectionPoint;
  EXPECT_FALSE(find_intersection_of_three_planes_within_tet(plane0, plane1, plane2, intersectionPoint));
}

void expect_no_intersection_with_side_of_tet(const int side, const stk::math::Plane3d & plane0, const stk::math::Plane3d & plane1)
{
  stk::math::Vector3d intersectionPoint;
  EXPECT_FALSE(find_intersection_of_two_planes_and_side_of_tet(side, plane0, plane1, intersectionPoint));
}

struct IntersectPlanes : public ::testing::Test
{
  const stk::math::Vector3d pt0{0., 0., 0.};
  const stk::math::Vector3d pt1{1., 0., 0.};
  const stk::math::Vector3d pt2{0., 1., 0.};
  const stk::math::Vector3d pt3{0., 0., 1.};
  const stk::math::Plane3d plane0{pt0, pt1, pt3};
  const stk::math::Plane3d plane1{pt1, pt2, pt3};
  const stk::math::Plane3d plane2{pt2, pt0, pt3};
  const stk::math::Plane3d plane3{pt0, pt2, pt1};

  const stk::math::Vector3d offsetPt0{-0.1, -0.1, -0.1};
  const stk::math::Vector3d offsetPt1{1.1, -0.1, -0.1};
  const stk::math::Vector3d offsetPt2{-0.1, 1.1, -0.1};
  const stk::math::Vector3d offsetPt3{-0.1, -0.1, 1.1};
  const stk::math::Plane3d offsetPlane0{offsetPt0, offsetPt1, offsetPt3};
  const stk::math::Plane3d offsetPlane1{offsetPt1, offsetPt2, offsetPt3};
  const stk::math::Plane3d offsetPlane2{offsetPt2, offsetPt0, offsetPt3};
  const stk::math::Plane3d offsetPlane3{offsetPt0, offsetPt2, offsetPt1};

};

TEST_F(IntersectPlanes, given3Planes_FindCorrectIntersection)
{
  expect_intersection(plane0, plane1, plane2, stk::math::Vector3d{0., 0., 1.});
  expect_intersection(plane0, plane1, plane3, stk::math::Vector3d{1., 0., 0.});
  expect_intersection(plane1, plane2, plane3, stk::math::Vector3d{0., 1., 0.});
  expect_intersection(plane0, plane2, plane3, stk::math::Vector3d{0., 0., 0.});

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
  expect_intersection_with_side_of_tet(0, plane1, plane2, stk::math::Vector3d{0., 0., 1.});
  expect_intersection_with_side_of_tet(0, plane1, plane3, stk::math::Vector3d{1., 0., 0.});
  expect_intersection_with_side_of_tet(0, plane2, plane3, stk::math::Vector3d{0., 0., 0.});

  expect_intersection_with_side_of_tet(1, plane2, plane3, stk::math::Vector3d{0., 1., 0.});
  expect_intersection_with_side_of_tet(1, plane0, plane2, stk::math::Vector3d{0., 0., 1.});
  expect_intersection_with_side_of_tet(1, plane0, plane3, stk::math::Vector3d{1., 0., 0.});

  expect_intersection_with_side_of_tet(2, plane0, plane1, stk::math::Vector3d{0., 0., 1.});
  expect_intersection_with_side_of_tet(2, plane0, plane3, stk::math::Vector3d{0., 0., 0.});
  expect_intersection_with_side_of_tet(2, plane1, plane3, stk::math::Vector3d{0., 1., 0.});

  expect_intersection_with_side_of_tet(3, plane0, plane1, stk::math::Vector3d{1., 0., 0.});
  expect_intersection_with_side_of_tet(3, plane0, plane2, stk::math::Vector3d{0., 0., 0.});
  expect_intersection_with_side_of_tet(3, plane1, plane2, stk::math::Vector3d{0., 1., 0.});
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

static stk::math::Vector3d compute_2d_plane_direction(const stk::math::Vector3d & pt0, const stk::math::Vector3d & pt1)
{
  return stk::math::Vector3d(pt1[1]-pt0[1],pt0[0]-pt1[0],0.);
}

static stk::math::Plane3d build_2d_plane(const stk::math::Vector3d & pt0, const stk::math::Vector3d & pt1)
{
  return stk::math::Plane3d(compute_2d_plane_direction(pt0,pt1),pt0);
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

  const stk::math::Vector3d pt0{0., 0., 0.};
  const stk::math::Vector3d pt1{1., 0., 0.};
  const stk::math::Vector3d pt2{0., 1., 0.};

  const stk::math::Vector3d offsetPt0{-0.1, -0.1, 0.};
  const stk::math::Vector3d offsetPt1{1.1, -0.1, 0.};
  const stk::math::Vector3d offsetPt2{-0.1, 1.1, 0.};

  stk::math::Plane3d plane0;
  stk::math::Plane3d plane1;
  stk::math::Plane3d plane2;

  stk::math::Plane3d offsetPlane0;
  stk::math::Plane3d offsetPlane1;
  stk::math::Plane3d offsetPlane2;
};

void expect_2d_intersection(const stk::math::Plane3d & plane0, const stk::math::Plane3d & plane1, const stk::math::Vector3d & goldIntersection)
{
  stk::math::Vector3d intersectionPoint;
  EXPECT_TRUE(find_intersection_of_two_2D_planes(plane0, plane1, intersectionPoint));
  expect_eq(goldIntersection, intersectionPoint);
}

void expect_no_2d_intersection_within_tri(const stk::math::Plane3d & plane0, const stk::math::Plane3d & plane1)
{
  stk::math::Vector3d intersectionPoint;
  EXPECT_FALSE(find_intersection_of_two_2D_planes_within_tri(plane0, plane1, intersectionPoint));
}

void expect_no_2d_intersection(const stk::math::Plane3d & plane0, const stk::math::Plane3d & plane1)
{
  stk::math::Vector3d intersectionPoint;
  EXPECT_FALSE(find_intersection_of_two_2D_planes(plane0, plane1, intersectionPoint));
}

TEST_F(Intersect2DPlanes, given2Plane2Ds_FindCorrectIntersection)
{
  expect_2d_intersection(plane0, plane1, stk::math::Vector3d{1., 0., 0.});
  expect_2d_intersection(plane1, plane2, stk::math::Vector3d{0., 1., 0.});
  expect_2d_intersection(plane0, plane2, stk::math::Vector3d{0., 0., 0.});

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
  const stk::math::Vector3d pt0{0., 0., 0.};
  const stk::math::Vector3d pt1{1., 0., 0.};
  const stk::math::Vector3d pt2{0., 1., 1.};
  const std::array<stk::math::Vector3d,3> triCoords{{ pt0, pt1, pt2 }};

  expect_eq(stk::math::Vector3d{0.,0.,0.}, triangle_parametric_coordinates_of_projected_point(triCoords, pt0));
  expect_eq(stk::math::Vector3d{1.,0.,0.}, triangle_parametric_coordinates_of_projected_point(triCoords, pt1));
  expect_eq(stk::math::Vector3d{0.,1.,0.}, triangle_parametric_coordinates_of_projected_point(triCoords, pt2));

  const stk::math::Vector3d normal = CalcTriangle3<double>::normal(pt0, pt1, pt2);
  expect_eq(stk::math::Vector3d{0.,0.,0.}, triangle_parametric_coordinates_of_projected_point(triCoords, pt0+normal));
  expect_eq(stk::math::Vector3d{1.,0.,0.}, triangle_parametric_coordinates_of_projected_point(triCoords, pt1+normal));
  expect_eq(stk::math::Vector3d{0.,1.,0.}, triangle_parametric_coordinates_of_projected_point(triCoords, pt2+normal));

  const double oneThird = 1./3.;
  const stk::math::Vector3d midPt = oneThird*(pt0+pt1+pt2);
  expect_eq(stk::math::Vector3d{oneThird,oneThird,0.}, triangle_parametric_coordinates_of_projected_point(triCoords, midPt+normal));
}

} // namespace krino
