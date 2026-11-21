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
#include <functional>

namespace krino {

// A sum over processors rounded to float, so that it is the same for every distribution of its addends
// across processors (but for sums at a float rounding boundary): the last bits of a sum of doubles depend
// on the order of its addends once they span more than the double mantissa.  This suits sums of many
// addends of widely varying size whose result needs only about float precision, as the dot products of
// an optimizer.  It does not suit sums that cancel to the roundoff of their addends, whose roundoff
// depends on the order, nor results that need more than float precision, as an objective compared in a
// line search: rounding each addend to float instead keeps such a sum exact, and so independent of the
// order, while its addends span less than about 2^29 in magnitude, and otherwise order dependent only at
// the level of its smallest addends (as MeshObjectiveBase does).
inline double round_sum_for_distribution_independence(const double sum) { return float(sum); }

double compute_parametric_square_distance(const stk::math::Vector3d childPCoords);

template<typename NODECOORDCONTAINER>
stk::math::Vector3d get_parametric_coordinates_of_point(const NODECOORDCONTAINER & nodeCoords, const stk::math::Vector3d & pt);

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
