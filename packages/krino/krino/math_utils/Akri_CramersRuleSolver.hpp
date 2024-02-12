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


namespace krino {
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
