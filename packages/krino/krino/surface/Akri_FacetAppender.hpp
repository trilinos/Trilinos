// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#ifndef KRINO_KRINO_SURFACE_AKRI_FACETAPPENDER_HPP_
#define KRINO_KRINO_SURFACE_AKRI_FACETAPPENDER_HPP_

#include <stk_math/StkVector.hpp>

#include <cstddef>

namespace krino
{

// Where facets are written one at a time by the coordinates of their vertices (segments in 2D,
// triangles in 3D), as by contouring, without the writer knowing how they are stored.
class FacetAppender
{
public:
  virtual ~FacetAppender() = default;
  virtual void emplace_back_2d(const stk::math::Vector3d & x0, const stk::math::Vector3d & x1) = 0;
  virtual void emplace_back_3d(const stk::math::Vector3d & x0,
      const stk::math::Vector3d & x1,
      const stk::math::Vector3d & x2) = 0;
  virtual void clear() = 0;
  // The number of facets it holds.
  virtual size_t num_facets() const = 0;
};

} // namespace krino

#endif /* KRINO_KRINO_SURFACE_AKRI_FACETAPPENDER_HPP_ */
