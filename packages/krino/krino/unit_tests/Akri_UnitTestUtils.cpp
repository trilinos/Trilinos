// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include <Akri_UnitTestUtils.hpp>
#include <stk_math/StkVector.hpp>
#include <Akri_Facet.hpp>
#include <gtest/gtest.h>

namespace krino {

static bool is_near_absolute(const stk::math::Vector3d & gold, const stk::math::Vector3d & result, const double absoluteTol)
{
  for (int i=0; i<3; ++i)
    if (std::abs(gold[i]-result[i]) > absoluteTol)
      return false;
  return true;
}

static bool is_near_relative(const stk::math::Vector3d & gold, const stk::math::Vector3d & result, const double relativeTol)
{
  const double absoluteTol = relativeTol * (gold.length() + result.length());
  return is_near_absolute(gold, result, absoluteTol);
}

void expect_eq(const stk::math::Vector3d & gold, const stk::math::Vector3d & result, const double relativeTol)
{
  EXPECT_TRUE(is_near_relative(gold, result, relativeTol)) << "Failed vector comparison: gold: " << gold << " actual:" << result << " relative tol:" << relativeTol;
}

void expect_eq_absolute(const stk::math::Vector3d & gold, const stk::math::Vector3d & result, const double absoluteTol)
{
  EXPECT_TRUE(is_near_absolute(gold, result, absoluteTol)) << "Failed vector comparison: gold: " << gold << " actual:" << result << " absolute tol:" << absoluteTol;;
}

bool is_near_relative(const Facet2d & gold, const Facet2d & result, const double relativeTol)
{
  return is_near_relative(gold.facet_vertex(0), result.facet_vertex(0), relativeTol) &&
      is_near_relative(gold.facet_vertex(1), result.facet_vertex(1), relativeTol);
}

void expect_eq(const Facet2d & gold, const Facet2d & result, const double relativeTol)
{
  EXPECT_TRUE(is_near_relative(gold, result, relativeTol)) << "Failed Facet2d comparison: gold: " << gold << " actual:" << result << " relative tol:" << relativeTol;;
}

bool is_near_relative(const Facet3d & gold, const Facet3d & result, const double relativeTol)
{
  const double err0 = (gold.facet_vertex(0) - result.facet_vertex(0)).length_squared();
  const double err1 = (gold.facet_vertex(1) - result.facet_vertex(0)).length_squared();
  const double err2 = (gold.facet_vertex(2) - result.facet_vertex(0)).length_squared();
  if (err1 < err0 && err1 < err2)
  {
    return is_near_relative(gold.facet_vertex(1), result.facet_vertex(0), relativeTol) &&
      is_near_relative(gold.facet_vertex(2), result.facet_vertex(1), relativeTol) &&
      is_near_relative(gold.facet_vertex(0), result.facet_vertex(2), relativeTol);
  }
  else if (err2 < err0 && err2 < err1)
  {
    return is_near_relative(gold.facet_vertex(2), result.facet_vertex(0), relativeTol) &&
      is_near_relative(gold.facet_vertex(0), result.facet_vertex(1), relativeTol) &&
      is_near_relative(gold.facet_vertex(1), result.facet_vertex(2), relativeTol);
  }

  return is_near_relative(gold.facet_vertex(0), result.facet_vertex(0), relativeTol) &&
    is_near_relative(gold.facet_vertex(1), result.facet_vertex(1), relativeTol) &&
    is_near_relative(gold.facet_vertex(2), result.facet_vertex(2), relativeTol);
}

void expect_eq(const Facet3d & gold, const Facet3d & result, const double relativeTol)
{
  EXPECT_TRUE(is_near_relative(gold, result, relativeTol)) << "Failed Facet3d comparison: gold: " << gold << " actual:" << result << " relative tol:" << relativeTol;;
}

}


