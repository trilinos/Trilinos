// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include <Akri_TransformedSurface.hpp>

#include <Akri_AnalyticSurf.hpp>

#include <stk_util/util/ReportHandler.hpp>

#include <cmath>
#include <limits>
#include <sstream>

namespace krino {

// The box around the given box's corners mapped by fn, or the entire domain if the given box is
// unbounded.
template <class FN>
static BoundingBox box_around_mapped_corners(const BoundingBox & bbox, const FN & fn)
{
  if (!bbox.valid())
    return bbox;
  const auto & lo = bbox.get_min();
  const auto & hi = bbox.get_max();
  for (unsigned i = 0; i < 3; ++i)
    if (!std::isfinite(lo[i]) || !std::isfinite(hi[i]) || std::abs(lo[i]) == std::numeric_limits<BoundingBox::Real>::max() || std::abs(hi[i]) == std::numeric_limits<BoundingBox::Real>::max())
      return BoundingBox::ENTIRE_DOMAIN;

  BoundingBox mapped;
  for (int i = 0; i < 2; ++i)
    for (int j = 0; j < 2; ++j)
      for (int k = 0; k < 2; ++k)
        mapped.accommodate(fn(stk::math::Vector3d(i ? hi[0] : lo[0], j ? hi[1] : lo[1], k ? hi[2] : lo[2])));
  return mapped;
}

TransformedSurface::TransformedSurface(std::shared_ptr<Surface> surface, const Transformation & motion)
: mySurface(std::move(surface)),
  myMotion(motion),
  myPlacement(motion.rigid_transformation_at(0.))
{
  STK_ThrowRequire(mySurface);
  const Plane * plane = dynamic_cast<const Plane *>(mySurface.get());
  STK_ThrowRequireMsg(nullptr == plane || !plane->needs_offset_from_extent(),
      "Motion is not supported for a plane placed at a fraction of the extent of the mesh.");
}

std::string TransformedSurface::output() const
{
  std::ostringstream os;
  os << "moving " << mySurface->output();
  return os.str();
}

void TransformedSurface::prepare_to_compute(const double time, const BoundingBox & point_bbox, const double truncation_length)
{
  myPlacement = myMotion.rigid_transformation_at(time);
  const BoundingBox surfacePointBbox = box_around_mapped_corners(point_bbox,
      [&](const stk::math::Vector3d & x) { return myPlacement.inverse_transform_point(x); });
  mySurface->prepare_to_compute(time, surfacePointBbox, truncation_length);
}

void TransformedSurface::insert_into(BoundingBox & bbox) const
{
  BoundingBox surfaceBbox;
  mySurface->insert_into(surfaceBbox);
  bbox.accommodate(box_around_mapped_corners(surfaceBbox,
      [&](const stk::math::Vector3d & x) { return myPlacement.transform_point(x); }));
}

bool TransformedSurface::does_intersect(const BoundingBox & bbox) const
{
  return mySurface->does_intersect(box_around_mapped_corners(bbox,
      [&](const stk::math::Vector3d & x) { return myPlacement.inverse_transform_point(x); }));
}

double TransformedSurface::point_signed_distance(const stk::math::Vector3d &x) const
{
  return mySurface->point_signed_distance(myPlacement.inverse_transform_point(x));
}

double TransformedSurface::possibly_truncated_point_signed_distance(const stk::math::Vector3d &x, const double truncation_length, const double far_field_value) const
{
  return mySurface->possibly_truncated_point_signed_distance(myPlacement.inverse_transform_point(x), truncation_length, far_field_value);
}

stk::math::Vector3d TransformedSurface::closest_point_normal(const stk::math::Vector3d &x) const
{
  return myPlacement.transform_vector(mySurface->closest_point_normal(myPlacement.inverse_transform_point(x)));
}

std::pair<int, double> TransformedSurface::compute_intersection_with_segment(const stk::math::Vector3d &pt0, const stk::math::Vector3d &pt1, const double edgeCrossingTol) const
{
  // The crossing is the same fraction of the way along the segment either way.
  return mySurface->compute_intersection_with_segment(myPlacement.inverse_transform_point(pt0), myPlacement.inverse_transform_point(pt1), edgeCrossingTol);
}

void TransformedSurface::fill_segment_crossings(const stk::math::Vector3d & pt0,
    const stk::math::Vector3d & pt1,
    const double edgeCrossingTol,
    std::vector<double> & crossingPositions) const
{
  mySurface->fill_segment_crossings(myPlacement.inverse_transform_point(pt0),
      myPlacement.inverse_transform_point(pt1),
      edgeCrossingTol,
      crossingPositions);
}

void TransformedSurface::fill_triangle_intersection_parametric_coordinates(const std::array<stk::math::Vector3d,3> & faceNodes, std::vector<stk::math::Vector3d> & intersectionParamCoords) const
{
  // Parametric coordinates are the same either way.
  const std::array<stk::math::Vector3d,3> surfaceFaceNodes{{myPlacement.inverse_transform_point(faceNodes[0]),
      myPlacement.inverse_transform_point(faceNodes[1]), myPlacement.inverse_transform_point(faceNodes[2])}};
  mySurface->fill_triangle_intersection_parametric_coordinates(surfaceFaceNodes, intersectionParamCoords);
}

void TransformedSurface::fill_tetrahedon_intersection_parametric_coordinates(const std::array<stk::math::Vector3d,4> & tetNodes, std::vector<stk::math::Vector3d> & intersectionParamCoords) const
{
  const std::array<stk::math::Vector3d,4> surfaceTetNodes{{myPlacement.inverse_transform_point(tetNodes[0]),
      myPlacement.inverse_transform_point(tetNodes[1]), myPlacement.inverse_transform_point(tetNodes[2]), myPlacement.inverse_transform_point(tetNodes[3])}};
  mySurface->fill_tetrahedon_intersection_parametric_coordinates(surfaceTetNodes, intersectionParamCoords);
}

} // namespace krino
