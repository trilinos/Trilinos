// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include "Akri_CramersRuleSolver.hpp"
#include <array>
#include <stk_util/util/ReportHandler.hpp>


namespace krino {
namespace CramersRuleSolver {

std::array<double,3> solve3x3(
        double a11, double a12, double a13,
        double a21, double a22, double a23,
        double a31, double a32, double a33,
        double k1,  double k2,  double k3 )
{
    double det = compute_determinant3x3( a11, a12, a13, a21, a22, a23, a31, a32, a33 );

    STK_ThrowRequireMsg(det != 0.0, "Bad determinant.  Are the points really unique?");

    std::array<double,3> answer;

    answer[0] = compute_determinant3x3( k1,  a12, a13,  k2, a22, a23,  k3, a32, a33 )/det;
    answer[1] = compute_determinant3x3( a11,  k1, a13, a21,  k2, a23, a31,  k3, a33 )/det;
    answer[2] = compute_determinant3x3( a11, a12,  k1, a21, a22,  k2, a31, a32,  k3 )/det;
    return( answer );
}

std::array<double,3> solve3x3(const std::array<std::array<double,3>,3> &A, std::array<double,3> & b)
{
    return solve3x3(A[0][0], A[0][1], A[0][2], A[1][0], A[1][1], A[1][2], A[2][0], A[2][1], A[2][2], b[0], b[1], b[2]);
}

std::array<double,5> solve5x5(
        double a11, double a12, double a13, double a14, double a15,
        double a21, double a22, double a23, double a24, double a25,
        double a31, double a32, double a33, double a34, double a35,
        double a41, double a42, double a43, double a44, double a45,
        double a51, double a52, double a53, double a54, double a55,
        double k1,  double k2,  double k3,  double k4,  double k5 )
{
    const double det = compute_determinant5x5(
        a11, a12, a13, a14, a15,
        a21, a22, a23, a24, a25,
        a31, a32, a33, a34, a35,
        a41, a42, a43, a44, a45,
        a51, a52, a53, a54, a55 );

    STK_ThrowRequireMsg(det != 0.0, "Bad determinant.  Are the points really unique?");

    std::array<double,5> answer;
    answer[0] = compute_determinant5x5(
        k1, a12, a13, a14, a15,
        k2, a22, a23, a24, a25,
        k3, a32, a33, a34, a35,
        k4, a42, a43, a44, a45,
        k5, a52, a53, a54, a55)/det;
    answer[1] = compute_determinant5x5(
        a11, k1, a13, a14, a15,
        a21, k2, a23, a24, a25,
        a31, k3, a33, a34, a35,
        a41, k4, a43, a44, a45,
        a51, k5, a53, a54, a55)/det;
    answer[2] = compute_determinant5x5(
        a11, a12, k1, a14, a15,
        a21, a22, k2, a24, a25,
        a31, a32, k3, a34, a35,
        a41, a42, k4, a44, a45,
        a51, a52, k5, a54, a55)/det;
    answer[3] = compute_determinant5x5(
        a11, a12, a13, k1, a15,
        a21, a22, a23, k2, a25,
        a31, a32, a33, k3, a35,
        a41, a42, a43, k4, a45,
        a51, a52, a53, k5, a55)/det;
    answer[4] = compute_determinant5x5(
        a11, a12, a13, a14, k1,
        a21, a22, a23, a24, k2,
        a31, a32, a33, a34, k3,
        a41, a42, a43, a44, k4,
        a51, a52, a53, a54, k5)/det;
    return( answer );
}

std::array<double,5> solve5x5(const std::array<std::array<double,5>,5> &A, std::array<double,5> & b)
{
    return solve5x5(
        A[0][0], A[0][1], A[0][2], A[0][3], A[0][4],
        A[1][0], A[1][1], A[1][2], A[1][3], A[1][4],
        A[2][0], A[2][1], A[2][2], A[2][3], A[2][4],
        A[3][0], A[3][1], A[3][2], A[3][3], A[3][4],
        A[4][0], A[4][1], A[4][2], A[4][3], A[4][4],
        b[0], b[1], b[2], b[3], b[4]);
}


double compute_determinant3x3(
        double a11, double a12, double a13,
        double a21, double a22, double a23,
        double a31, double a32, double a33 )
{
    return( a11*a22*a33 + a12*a23*a31 + a13*a21*a32 -
            a13*a22*a31 - a12*a21*a33 - a11*a23*a32 );
}

double compute_determinant4x4(
        double a11, double a12, double a13, double a14,
        double a21, double a22, double a23, double a24,
        double a31, double a32, double a33, double a34,
        double a41, double a42, double a43, double a44 )
{
    return( a11*compute_determinant3x3(a22, a23, a24, a32, a33, a34, a42, a43, a44) -
            a12*compute_determinant3x3(a21, a23, a24, a31, a33, a34, a41, a43, a44) +
            a13*compute_determinant3x3(a21, a22, a24, a31, a32, a34, a41, a42, a44) -
            a14*compute_determinant3x3(a21, a22, a23, a31, a32, a33, a41, a42, a43));
}

double compute_determinant5x5(
        double a11, double a12, double a13, double a14, double a15,
        double a21, double a22, double a23, double a24, double a25,
        double a31, double a32, double a33, double a34, double a35,
        double a41, double a42, double a43, double a44, double a45,
        double a51, double a52, double a53, double a54, double a55 )
{
    return( a11*compute_determinant4x4(a22, a23, a24, a25, a32, a33, a34, a35, a42, a43, a44, a45, a52, a53, a54, a55 ) -
            a12*compute_determinant4x4(a21, a23, a24, a25, a31, a33, a34, a35, a41, a43, a44, a45, a51, a53, a54, a55 ) +
            a13*compute_determinant4x4(a21, a22, a24, a25, a31, a32, a34, a35, a41, a42, a44, a45, a51, a52, a54, a55 ) -
            a14*compute_determinant4x4(a21, a22, a23, a25, a31, a32, a33, a35, a41, a42, a43, a45, a51, a52, a53, a55 ) +
            a15*compute_determinant4x4(a21, a22, a23, a24, a31, a32, a33, a34, a41, a42, a43, a44, a51, a52, a53, a54 ));
}

} // namespace CramersRuleSolver
} // namespace krino
