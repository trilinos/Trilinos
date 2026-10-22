// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.
#ifndef KRINO_KRINO_KRINO_LIB_AKRI_CONTOURUTILS_HPP_
#define KRINO_KRINO_KRINO_LIB_AKRI_CONTOURUTILS_HPP_
#include <array>

#include <stk_math/StkVector.hpp>

namespace krino {

template <size_t NCOORDS>
stk::math::Vector3d compute_linear_edge_crossing(const std::array<stk::math::Vector3d,NCOORDS> & coords,
  const std::array<double,NCOORDS> & distance,
  const unsigned i0, const unsigned i1);

template <size_t NCOORDS, size_t NDIST>
stk::math::Vector3d compute_quadratic_edge_crossing(const std::array<stk::math::Vector3d,NCOORDS> & coords,
  const std::array<double,NDIST> & distance,
  const unsigned i0, const unsigned i1, const unsigned i2);

}

#endif /* KRINO_KRINO_KRINO_LIB_AKRI_CONTOURUTILS_HPP_ */
