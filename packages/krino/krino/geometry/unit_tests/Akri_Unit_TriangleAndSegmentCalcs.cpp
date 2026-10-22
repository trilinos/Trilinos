// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include <gtest/gtest.h>

#include <Akri_Segment.hpp>
#include <Akri_SegmentWithSensitivities.hpp>
#include <Akri_Triangle.hpp>
#include <Akri_TriangleWithSensitivities.hpp>
#include <Akri_UnitTestUtils.hpp>

#include <stk_math/StkVector.hpp>
#include <stk_util/diag/Timer.hpp>
#include <stk_util/parallel/Parallel.hpp>

#include <random>

namespace krino
{

// These tests exercise CalcTriangle3 / CalcSegment3 and their forward-mode
// (FAD) sensitivity counterparts TriangleWithSens / SegmentWithSens. The
// analytic_detail helpers below independently compute the area/normal/length
// gradients (contributed by Sean Hardesty) so the FAD sensitivities can be
// checked against a hand-rolled analytic reference. Moved here from the
// krino_lib Akri_Unit_Geometry monolith so the tests live with the geometry
// library whose classes they cover.

namespace analytic_detail
{
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

}
