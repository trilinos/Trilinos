// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include <Akri_Surface.hpp>
#include <Akri_CubitJournal.hpp>
#include <Akri_SurfaceIntersectionFromSignedDistance.hpp>
#include <Akri_DiagWriter.hpp>
#include <stk_util/util/ReportHandler.hpp>
#include <cmath>
#include <ostream>

namespace krino{

void
Surface::prepare_to_compute(const double /*time*/, const BoundingBoxType & /*point_bbox*/, const double /*truncation_length*/)
{
}

void Surface::insert_into(BoundingBoxType & /*bbox*/) const
{
  STK_ThrowErrorMsg("This surface with type (" << type() << ") has not implemented insert_into().");
}

bool Surface::does_intersect(const BoundingBoxType & /*bbox*/) const
{
  STK_ThrowErrorMsg("This surface with type (" << type() << ") has not implemented does_intersect().");
  return false; // unreachable; STK_ThrowErrorMsg does not return
}

std::pair<int, double> Surface::compute_intersection_with_segment(const stk::math::Vector3d &pt0, const stk::math::Vector3d &pt1, const double edgeCrossingTol) const
{
  return compute_surface_intersection_with_segment_from_signed_distance(*this, pt0, pt1, edgeCrossingTol);
}

void Surface::fill_segment_crossings(const stk::math::Vector3d & pt0,
    const stk::math::Vector3d & pt1,
    const double edgeCrossingTol,
    std::vector<double> & crossingPositions) const
{
  crossingPositions.clear();
  // The distance at an end on the surface, such as a node snapped onto it, is roundoff of either
  // sign.  That sign would decide whether a crossing beyond the end is found: the segment can run
  // to the other side from the end and cross back, which the signs at its ends show only if the
  // roundoff happens to have the sign of the other side.  So the sign of such an end is taken from
  // the point the tolerance into the segment instead.
  const double onSurfaceTol = edgeCrossingTol * (pt1 - pt0).length();
  const double begin = (std::abs(point_signed_distance(pt0)) <= onSurfaceTol) ? edgeCrossingTol : 0.;
  const double end = (std::abs(point_signed_distance(pt1)) <= onSurfaceTol) ? 1. - edgeCrossingTol : 1.;
  if (begin >= end) return;
  const auto [crossingSign, positionInPart] = compute_intersection_with_segment(
      (1. - begin) * pt0 + begin * pt1, (1. - end) * pt0 + end * pt1, edgeCrossingTol);
  if (crossingSign == 0) return;
  const double position = begin + positionInPart * (end - begin);
  if (position > edgeCrossingTol && position < 1. - edgeCrossingTol)
    crossingPositions.push_back(position);
}

void Surface::add_to_journal(CubitJournal & journal, const BoundingBoxType & region) const
{
  constexpr int numCellsPerSide = 16;
  const stk::math::Vector3d boxMin(region.get_min()[0], region.get_min()[1], region.get_min()[2]);
  const stk::math::Vector3d boxMax(region.get_max()[0], region.get_max()[1], region.get_max()[2]);
  journal.add_comment("Approximate: the zero contour of the signed distance of a surface of type " +
      std::to_string(type()) + " on a grid of " + std::to_string(numCellsPerSide) + " cells a side.");
  add_zero_contour_in_box(journal,
      [this](const stk::math::Vector3d & x) { return point_signed_distance(x); },
      boxMin, boxMax, numCellsPerSide);
}

stk::math::Vector3d Surface::closest_point_normal(const stk::math::Vector3d &x) const
{
  STK_ThrowErrorMsg("This surface with type (" << type() << ") has not implemented closest_point_normal().");
  return stk::math::Vector3d::ZERO; // unreachable; STK_ThrowErrorMsg does not return
}

std::string Surface::output() const
{
  std::ostringstream os;
  os << "surface with type " << type();
  return os.str();
}

} // namespace krino
