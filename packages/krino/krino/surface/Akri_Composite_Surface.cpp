// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include <Akri_Composite_Surface.hpp>
#include <Akri_BoundingBox.hpp>
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
{
  
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
      if (FACETED_SURFACE_2D == surface->type() ||
          FACETED_SURFACE_3D == surface->type() ||
          nullptr != surface->get_transformation() ||
          surface->does_intersect(point_bbox))
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
Composite_Surface::truncated_point_signed_distance(const stk::math::Vector3d &x, const double narrow_band, const double far_field_value) const
{

  if (my_composition_method == MINIMUM_SIGNED_DISTANCE)
  {
    STK_ThrowRequireMsg(far_field_value >= narrow_band, "Composite surfaces have a specific requirement for far_field_value due to min/max operations.");
    double dist = (narrow_band == 0.) ? std::numeric_limits<double>::max() : far_field_value;
    for ( auto&& surface : my_subsurfaces )
    {
      dist = std::min(dist, surface->truncated_point_signed_distance(x, narrow_band, far_field_value));
    }
    return dist;
  }
  else
  {
    STK_ThrowRequire(my_composition_method == MAXIMUM_SIGNED_DISTANCE);
    STK_ThrowRequireMsg(far_field_value <= -narrow_band, "Composite surfaces have a specific requirement for far_field_value due to min/max operations.");
    double dist = (narrow_band == 0.) ? -std::numeric_limits<double>::max() : far_field_value;
    for ( auto&& surface : my_subsurfaces )
    {
      dist = std::max(dist, surface->truncated_point_signed_distance(x, narrow_band, far_field_value));
    }
    return dist;
  }
}

void Composite_Surface::insert_into(BoundingBox & bbox) const
{
  for (auto && subsurf : my_subsurfaces)
    subsurf->insert_into(bbox);
}

bool Composite_Surface::does_intersect(const BoundingBox & bbox) const
{
  for (auto && subsurf : my_subsurfaces)
    if (subsurf->does_intersect(bbox))
      return true;
  return false;
}


} // namespace krino



