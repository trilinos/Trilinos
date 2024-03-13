// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#ifndef KRINO_UNIT_TESTS_INCLUDE_AKRI_UNITTESTUTILS_H_
#define KRINO_UNIT_TESTS_INCLUDE_AKRI_UNITTESTUTILS_H_
#include <stk_math/StkVector.hpp>

namespace krino {

class Facet2d;
class Facet3d;

void expect_eq(const stk::math::Vector3d & gold, const stk::math::Vector3d & result, const double relativeTol=1.e-6);
void expect_eq_absolute(const stk::math::Vector3d & gold, const stk::math::Vector3d & result, const double absoluteTol=1.e-6);
void expect_eq(const Facet2d & gold, const Facet2d & result, const double relativeTol=1.e-6);
void expect_eq(const Facet3d & gold, const Facet3d & result, const double relativeTol=1.e-6);
bool is_near_relative(const Facet2d & gold, const Facet2d & result, const double relativeTol=1.e-6);
bool is_near_relative(const Facet3d & gold, const Facet3d & result, const double relativeTol=1.e-6);

}

#endif /* KRINO_UNIT_TESTS_INCLUDE_AKRI_UNITTESTUTILS_H_ */
