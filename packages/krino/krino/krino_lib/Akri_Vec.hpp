// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#ifndef Akri_Vec_h
#define Akri_Vec_h

#include <stk_math/StkVector.hpp>

namespace krino {

using stk::math::MemberInit;
using stk::math::Vec;
typedef stk::math::Vec<double,3> Vector3d;
typedef stk::math::Vec<double,2> Vector2d;
typedef stk::math::Vec<float,3> Float3d;

bool is_less_than_in_x_then_y_then_z(const Vector3d& A, const Vector3d &B);

} // namespace krino

#endif // Akri_Vec_h
