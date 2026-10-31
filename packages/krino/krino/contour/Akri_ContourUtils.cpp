// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include <Akri_ContourUtils.hpp>

#include <Akri_MathUtil.hpp>
#include <stk_math/StkVector.hpp>
#include <stk_util/util/ReportHandler.hpp>

namespace krino {

template <size_t NCOORDS>
stk::math::Vector3d compute_linear_edge_crossing(const std::array<stk::math::Vector3d,NCOORDS> & coords,
  const std::array<double,NCOORDS> & distance,
  const unsigned i0, const unsigned i1)
{
  STK_ThrowAssert(distance[i0]*distance[i1] < 0.0); // Insist on a crossing
  const double loc = distance[i0]/(distance[i0]-distance[i1]);
  return (1.-loc) * coords[i0] + loc * coords[i1];
}

template <size_t NCOORDS, size_t NDIST>
stk::math::Vector3d compute_quadratic_edge_crossing(const std::array<stk::math::Vector3d,NCOORDS> & coords,
  const std::array<double,NDIST> & distance,
  const unsigned i0, const unsigned i1, const unsigned i2)
{
  const double loc = find_quadratic_crossing(distance[i0], distance[i1], distance[i2]);
  return (1.-loc) * coords[i0] + loc * coords[i1];
}

// Explicit template instantiation
template stk::math::Vector3d compute_linear_edge_crossing(const std::array<stk::math::Vector3d,3> & coords, const std::array<double,3> & distance, const unsigned i0, const unsigned i1);
template stk::math::Vector3d compute_linear_edge_crossing(const std::array<stk::math::Vector3d,4> & coords, const std::array<double,4> & distance, const unsigned i0, const unsigned i1);
template stk::math::Vector3d compute_quadratic_edge_crossing(const std::array<stk::math::Vector3d,3> & coords, const std::array<double,6> & distance, const unsigned i0, const unsigned i1, const unsigned i2);
template stk::math::Vector3d compute_quadratic_edge_crossing(const std::array<stk::math::Vector3d,4> & coords, const std::array<double,10> & distance, const unsigned i0, const unsigned i1, const unsigned i2);

}


