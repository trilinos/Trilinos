// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include <Akri_UnitTestUtils.hpp>
#include <Akri_Vec.hpp>
#include <gtest/gtest.h>

namespace krino {

void expect_eq(const Vector3d & gold, const Vector3d & result, const double relativeTol)
{
  const double absoluteTol = relativeTol * (gold.length() + result.length());
  expect_eq_absolute(gold, result, absoluteTol);
}

void expect_eq_absolute(const Vector3d & gold, const Vector3d & result, const double absoluteTol)
{
  for (int i=0; i<3; ++i)
    EXPECT_NEAR(gold[i], result[i], absoluteTol) <<"gold: " << gold << " actual:" << result;
}

}


