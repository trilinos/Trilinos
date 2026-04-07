// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#ifndef KRINO_KRINO_KRINO_LIB_AKRI_CRAMERSRULESOLVER_HPP_
#define KRINO_KRINO_KRINO_LIB_AKRI_CRAMERSRULESOLVER_HPP_
#include <array>
#include <cmath>


namespace krino {

template<typename MAT>
bool compute_determinant3x3(const MAT & A)
{
  return A[0]*(A[1][1]*A[2][2] - A[1][2]*A[2][1]) + A[1]*(A[1][2]*A[2][0] - A[1][0]*A[2][2]) + A[1]*(A[1][0]*A[2][1] - A[1][1]*A[2][0]);
}

template<typename MAT>
bool invert3x3(const MAT & A, MAT & invA)
{
  const std::array<double,3> A1xA2 = {A[1][1]*A[2][2] - A[1][2]*A[2][1], A[1][2]*A[2][0] - A[1][0]*A[2][2], A[1][0]*A[2][1] - A[1][1]*A[2][0]};
  const double det = A[0][0]*A1xA2[0] + A[0][1]*A1xA2[1] + A[0][2]*A1xA2[2];

  if(std::abs(det)<1e-16) return false;
  const double invDet = 1./det;

  invA[0][0] = A1xA2[0]*invDet;
  invA[0][1] = (A[2][1]*A[0][2] - A[2][2]*A[0][1])*invDet;
  invA[0][2] = (A[0][1]*A[1][2] - A[0][2]*A[1][1])*invDet;

  invA[1][0] = A1xA2[1]*invDet;
  invA[1][1] = (A[2][2]*A[0][0] - A[2][0]*A[0][2])*invDet;
  invA[1][2] = (A[0][2]*A[1][0] - A[0][0]*A[1][2])*invDet;

  invA[2][0] = A1xA2[2]*invDet;
  invA[2][1] = (A[2][0]*A[0][1] - A[2][1]*A[0][0])*invDet;
  invA[2][2] = (A[0][0]*A[1][1] - A[0][1]*A[1][0])*invDet;

  return true;
}

template<typename MAT>
void multiply_3x3_3x3(const MAT & fullA, const MAT & fullB, MAT & result)
{
  result[0][0] = fullA[0][0] * fullB[0][0] + fullA[0][1] * fullB[1][0] + fullA[0][2] * fullB[2][0];
  result[0][1] = fullA[0][0] * fullB[0][1] + fullA[0][1] * fullB[1][1] + fullA[0][2] * fullB[2][1];
  result[0][2] = fullA[0][0] * fullB[0][2] + fullA[0][1] * fullB[1][2] + fullA[0][2] * fullB[2][2];
  result[1][0] = fullA[1][0] * fullB[0][0] + fullA[1][1] * fullB[1][0] + fullA[1][2] * fullB[2][0];
  result[1][1] = fullA[1][0] * fullB[0][1] + fullA[1][1] * fullB[1][1] + fullA[1][2] * fullB[2][1];
  result[1][2] = fullA[1][0] * fullB[0][2] + fullA[1][1] * fullB[1][2] + fullA[1][2] * fullB[2][2];
  result[2][0] = fullA[2][0] * fullB[0][0] + fullA[2][1] * fullB[1][0] + fullA[2][2] * fullB[2][0];
  result[2][1] = fullA[2][0] * fullB[0][1] + fullA[2][1] * fullB[1][1] + fullA[2][2] * fullB[2][1];
  result[2][2] = fullA[2][0] * fullB[0][2] + fullA[2][1] * fullB[1][2] + fullA[2][2] * fullB[2][2];
}

template<typename MAT, typename VEC>
void multiply_3x3_vec3(const MAT & fullA, const VEC & b, VEC & result)
{
  result[0] = fullA[0][0] * b[0] + fullA[0][1] * b[1] + fullA[0][2] * b[2];
  result[1] = fullA[1][0] * b[0] + fullA[1][1] * b[1] + fullA[1][2] * b[2];
  result[2] = fullA[2][0] * b[0] + fullA[2][1] * b[1] + fullA[2][2] * b[2];
}

namespace CramersRuleSolver {

std::array<double,3> solve3x3(
        double a11, double a12, double a13,
        double a21, double a22, double a23,
        double a31, double a32, double a33,
        double k1,  double k2,  double k3 );

std::array<double,3> solve3x3(const std::array<std::array<double,3>,3> &A, std::array<double,3> & b);

std::array<double,5> solve5x5(
        double a11, double a12, double a13, double a14, double a15,
        double a21, double a22, double a23, double a24, double a25,
        double a31, double a32, double a33, double a34, double a35,
        double a41, double a42, double a43, double a44, double a45,
        double a51, double a52, double a53, double a54, double a55,
        double k1,  double k2,  double k3,  double k4,  double k5 );

std::array<double,5> solve5x5(const std::array<std::array<double,5>,5> &A, std::array<double,5> & b);

double compute_determinant3x3(
        double a11, double a12, double a13,
        double a21, double a22, double a23,
        double a31, double a32, double a33 );

double compute_determinant4x4(
        double a11, double a12, double a13, double a14,
        double a21, double a22, double a23, double a24,
        double a31, double a32, double a33, double a34,
        double a41, double a42, double a43, double a44 );

double compute_determinant5x5(
        double a11, double a12, double a13, double a14, double a15,
        double a21, double a22, double a23, double a24, double a25,
        double a31, double a32, double a33, double a34, double a35,
        double a41, double a42, double a43, double a44, double a45,
        double a51, double a52, double a53, double a54, double a55 );

} // namespace CramersRuleSolver
} // namespace krino

#endif /* KRINO_KRINO_KRINO_LIB_AKRI_CRAMERSRULESOLVER_HPP_ */
