// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#ifndef Akri_Surface_h
#define Akri_Surface_h

#include <Akri_BoundingBox.hpp>
#include <stk_math/StkVector.hpp>
#include <vector>
#include <memory>

namespace krino {

template<class T>
inline size_t
storage_size(const std::vector<T> & vec) {return (sizeof(std::vector<T>) + vec.capacity()*sizeof(T)); }


enum Surface_Type
{
  POINT=0,
  SPHERE,
  ELLIPSOID,
  CSG_SURFACE,
  CYLINDER,
  CUBOID,
  COMPOSITE_SURFACE,
  PLANE,
  RANDOM,
  STRING_FUNCTION,
  FACETED_SURFACE_2D,
  FACETED_SURFACE_3D,
  INR,
  TRANSFORMED_SURFACE,
  // Never, ever, ever add an entry after MAX_SURFACE_TYPE.  Never.
  MAX_SURFACE_TYPE
};

class CubitJournal;
class Surface;
typedef std::vector<Surface *> SurfaceVec;
typedef std::vector< std::unique_ptr<Surface> > SurfaceAutoVec;

// abstract class used to define a surface.

class Surface {

public:
  Surface() = default;
  virtual ~Surface() {}

  typedef BoundingBox BoundingBoxType;
  virtual void insert_into(BoundingBoxType & bbox) const;
  virtual bool does_intersect(const BoundingBoxType & bbox) const;

  // pre-calculations needed to compute distance
  virtual void prepare_to_compute(const double time, const BoundingBoxType & point_bbox, const double truncation_length);
  
  // compute signed distance from specific point to surface
  virtual double point_signed_distance(const stk::math::Vector3d &x) const = 0;
  
  // compute signed distance from specific point to surface
  // For distances larger than truncation_length (if truncation_length > 0.), the surface is allowed to return
  // far_field_value instead of the actual signed distance.
  virtual double possibly_truncated_point_signed_distance(const stk::math::Vector3d &x, const double truncation_length, const double far_field_value) const = 0;
  // If surface does return far_field_value instead of actual signed distance, the sign may be wrong.
  virtual bool truncated_distance_may_have_wrong_sign() const = 0;

  virtual bool can_approximate_closest_point_normal() const { return false; }
  virtual stk::math::Vector3d closest_point_normal(const stk::math::Vector3d &x) const;

  // For debugging: adds the surface in (or near) the region to the journal's current group.  By
  // default, an approximation: the zero contour of point_signed_distance on a grid of the region.
  virtual void add_to_journal(CubitJournal & journal, const BoundingBoxType & region) const;

  virtual std::pair<int, double> compute_intersection_with_segment(const stk::math::Vector3d &pt0, const stk::math::Vector3d &pt1, const double edgeCrossingTol) const;
  // Every position, as the fraction of the way from pt0 to pt1, where the segment crosses the
  // surface, in order, leaving out those within edgeCrossingTol of its ends.  Unlike
  // compute_intersection_with_segment these include crossings that the signs at the ends do not
  // show, as where a segment crosses the surface twice beside a sharp edge.  By default, the
  // crossing given by compute_intersection_with_segment, with the sign of an end within
  // edgeCrossingTol of the surface taken from the point that far into the segment.
  virtual void fill_segment_crossings(const stk::math::Vector3d & pt0,
      const stk::math::Vector3d & pt1,
      const double edgeCrossingTol,
      std::vector<double> & crossingPositions) const;

  // for debugging memory usage
  virtual size_t storage_size() const = 0;

  // Whether the surface is different at different times given to prepare_to_compute.
  virtual bool changes_with_time() const { return false; }

  virtual bool has_corners() const { return false; }
  virtual bool has_edges() const { return false; }
  bool has_sharp_features() const { return has_edges() || has_corners(); }
  virtual void fill_triangle_intersection_parametric_coordinates(const std::array<stk::math::Vector3d,3> & /*faceNodes*/, std::vector<stk::math::Vector3d> & /*intersectionParamCoords*/) const {}
  virtual void fill_tetrahedon_intersection_parametric_coordinates(const std::array<stk::math::Vector3d,4> & /*tetNodes*/, std::vector<stk::math::Vector3d> & /*intersectionParamCoords*/) const {}

  // queries
  virtual Surface_Type type() const = 0;
  virtual std::string output() const;
};

class SurfaceThatTakesAdvantageOfNarrowBandAndThereforeMightHaveWrongSign : public Surface {
public:
  virtual bool truncated_distance_may_have_wrong_sign() const override { return true; }
  virtual double point_signed_distance(const stk::math::Vector3d &x) const override
  {
    return possibly_truncated_point_signed_distance(x, 0., 0.);
  }
  virtual double possibly_truncated_point_signed_distance(const stk::math::Vector3d &x,
      const double truncationLength,
      const double farFieldValue) const override = 0;
};

class SurfaceThatDoesntTakeAdvantageOfNarrowBandAndThereforeHasCorrectSign : public Surface {
public:
  virtual bool truncated_distance_may_have_wrong_sign() const override { return false; }
  virtual double point_signed_distance(const stk::math::Vector3d &x) const override = 0;
  virtual double possibly_truncated_point_signed_distance(const stk::math::Vector3d &x, const double /*truncationLength*/, const double /*farFieldValue*/) const override
  {
    return point_signed_distance(x);
  }
};
} // namespace krino

#endif // Akri_Surface_h
