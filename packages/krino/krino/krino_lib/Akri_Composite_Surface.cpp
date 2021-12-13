// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include <Akri_Composite_Surface.hpp>
#include <Akri_DiagWriter.hpp>
#include <stk_util/environment/EnvData.hpp>
#include <stk_util/parallel/ParallelReduce.hpp>
#include <limits>

namespace krino{

Composite_Surface::Composite_Surface(const std::string & sn)
   : SurfaceThatTakesAdvantageOfNarrowBandAndThereforeMightHaveWrongSign(),
     my_name(sn),
     my_composition_method(MINIMUM_SIGNED_DISTANCE) {}

Composite_Surface::~Composite_Surface() {}

void 
Composite_Surface::prepare_to_compute(const double time, const BoundingBox & point_bbox, const double truncation_length)
{ /* %TRACE[ON]% */ Trace trace__("krino::Composite_Surface::prepare_to_compute(const double time, const BoundingBox & point_bbox, const double truncation_length)"); /* %TRACE% */
  
  for ( auto&& surface : my_subsurfaces )
  {
    surface->prepare_to_compute(time, point_bbox, truncation_length);
  }

  if (truncation_length > 0.)
  {
    BoundingBox padded_point_bbox = point_bbox;
    padded_point_bbox.pad(truncation_length);

    SurfaceAutoVec nearby_surfaces;
    for ( auto&& surface : my_subsurfaces )
    {
      // keep nearby or moving surfaces
      if (FACETED_SURFACE == surface->type() ||
          nullptr != surface->get_transformation() ||
          point_bbox.intersects(surface->get_bounding_box()))
      {
        nearby_surfaces.emplace_back(std::move(surface));
      }
    }
    nearby_surfaces.swap(my_subsurfaces);
  }

  int subsurfaces_might_produce_wrong_sign = false;
  for ( auto&& surface : my_subsurfaces )
  {
    if (surface->truncated_distance_may_have_wrong_sign())
    {
      subsurfaces_might_produce_wrong_sign = true;
      break;
    }
  }
  int local_subsurfaces_might_produce_wrong_sign = subsurfaces_might_produce_wrong_sign;
  stk::all_reduce_max(stk::EnvData::parallel_comm(), &local_subsurfaces_might_produce_wrong_sign, &subsurfaces_might_produce_wrong_sign, 1);
  my_subsurfaces_might_produce_wrong_sign = subsurfaces_might_produce_wrong_sign;
}

double
Composite_Surface::truncated_point_signed_distance(const Vector3d &x, const double narrow_band, const double far_field_value) const
{ /* %TRACE% */  /* %TRACE% */

  if (my_composition_method == MINIMUM_SIGNED_DISTANCE)
  {
    ThrowRequireMsg(far_field_value >= narrow_band, "Composite surfaces have a specific requirement for far_field_value due to min/max operations.");
    double dist = (narrow_band == 0.) ? std::numeric_limits<double>::max() : far_field_value;
    for ( auto&& surface : my_subsurfaces )
    {
      dist = std::min(dist, surface->truncated_point_signed_distance(x, narrow_band, far_field_value));
    }
    return dist;
  }
  else
  {
    ThrowRequire(my_composition_method == MAXIMUM_SIGNED_DISTANCE);
    ThrowRequireMsg(far_field_value <= -narrow_band, "Composite surfaces have a specific requirement for far_field_value due to min/max operations.");
    double dist = (narrow_band == 0.) ? -std::numeric_limits<double>::max() : far_field_value;
    for ( auto&& surface : my_subsurfaces )
    {
      dist = std::max(dist, surface->truncated_point_signed_distance(x, narrow_band, far_field_value));
    }
    return dist;
  }
}

BoundingBox
Composite_Surface::get_bounding_box()
{
  BoundingBox bbox;
  for (auto && subsurf : my_subsurfaces)
  {
    bbox.accommodate(subsurf->get_bounding_box());
  }
  return bbox;
}


} // namespace krino



