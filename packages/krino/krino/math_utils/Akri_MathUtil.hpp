// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#ifndef KRINO_INCLUDE_AKRI_MATHUTIL_H_
#define KRINO_INCLUDE_AKRI_MATHUTIL_H_
#include <stk_math/StkVector.hpp>
#include <vector>
#include <functional>

namespace krino {

double compute_parametric_square_distance(const stk::math::Vector3d childPCoords);

stk::math::Vector3d get_parametric_coordinates_of_point(const std::vector<stk::math::Vector3d> & nodeCoords, const stk::math::Vector3d & pt);

std::pair<bool, double> find_root( const std::function<double(const double)> & f,
    const double xa,
    const double xb,
    const double fa,
    const double fb,
    const unsigned maxIters = 100,
    const double tol = 1.e-4);

std::pair<bool, double> find_root_newton_raphson( const std::function<std::pair<double,double>(const double)> & f,
    const double guess,
    const unsigned maxIters = 100,
    const double fTol = 1.e-4);

double find_quadratic_crossing( const double d0, const double d1, const double d2 );

}


#endif /* KRINO_INCLUDE_AKRI_MATHUTIL_H_ */
