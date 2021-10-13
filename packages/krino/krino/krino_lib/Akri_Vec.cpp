// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include <Akri_Vec.hpp>

namespace krino {

static bool float_less(double a, double b)
{
  return static_cast<float>(a) < static_cast<float>(b);
}

bool is_less_than_in_x_then_y_then_z(const Vector3d& A, const Vector3d &B)
{
    if (float_less(A[0], B[0]))
        return true;
    else if (float_less(B[0], A[0]))
        return false;

    if (float_less(A[1], B[1]))
        return true;
    else if (float_less(B[1], A[1]))
        return false;

    if (float_less(A[2], B[2]))
        return true;
    else if (float_less(B[2], A[2]))
        return false;

    return false;
}

}
