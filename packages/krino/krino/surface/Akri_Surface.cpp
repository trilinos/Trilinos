// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include <Akri_Surface.hpp>
#include <Akri_SurfaceIntersectionFromSignedDistance.hpp>
#include <Akri_DiagWriter.hpp>
#include <stk_util/util/ReportHandler.hpp>

namespace krino{

void
Surface::prepare_to_compute(const double /*time*/, const BoundingBoxType & /*point_bbox*/, const double /*truncation_length*/)
{
  STK_ThrowErrorMsgIf(NULL != my_transformation,
      "This surface with type (" << type()
      << ") has motion specified, but the prepare_to_compute() method has not been implemented yet to support motion.");
}

void Surface::insert_into(BoundingBoxType & /*bbox*/) const
{
  ThrowRuntimeError("This surface with type (" << type() << ") has not implemented insert_into().");
}

bool Surface::does_intersect(const BoundingBoxType & /*bbox*/) const
{
  ThrowRuntimeError("This surface with type (" << type() << ") has not implemented does_intersect().");
}

std::pair<int, double> Surface::compute_intersection_with_segment(const stk::math::Vector3d &pt0, const stk::math::Vector3d &pt1, const double edgeCrossingTol) const
{
  return compute_surface_intersection_with_segment_from_signed_distance(*this, pt0, pt1, edgeCrossingTol);
}

} // namespace krino
