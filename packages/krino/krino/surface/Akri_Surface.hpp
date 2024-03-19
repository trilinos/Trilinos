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

class Transformation;

enum Surface_Type
{
  POINT=0,
  SPHERE,
  ELLIPSOID,
  CYLINDER,
  COMPOSITE_SURFACE,
  PLANE,
  RANDOM,
  STRING_FUNCTION,
  FACETED_SURFACE_2D,
  FACETED_SURFACE_3D,
  // Never, ever, ever add an entry after MAX_SURFACE_TYPE.  Never.
  MAX_SURFACE_TYPE
};

class Surface;
typedef std::vector<Surface *> SurfaceVec;
typedef std::vector< std::unique_ptr<Surface> > SurfaceAutoVec;

// abstract class used to define a surface.

class Surface {

public:
  Surface() { Surface::set_transformation(nullptr); }
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
  virtual double truncated_point_signed_distance(const stk::math::Vector3d &x, const double truncation_length, const double far_field_value) const = 0;
  // If surface does return far_field_value instead of actual signed distance, the sign may be wrong.
  virtual bool truncated_distance_may_have_wrong_sign() const = 0;

  virtual std::pair<int, double> compute_intersection_with_segment(const stk::math::Vector3d &pt0, const stk::math::Vector3d &pt1, const double edgeCrossingTol) const;

  // for debugging memory usage
  virtual size_t storage_size() const = 0;

  // methods related to moving surfaces (transformations)
  virtual void set_transformation(Transformation * trans) { my_transformation = trans; }
  Transformation * get_transformation() { return my_transformation; }
  const Transformation * get_transformation() const { return my_transformation; }

  // queries
  virtual Surface_Type type() const = 0;

protected:
  Transformation * my_transformation;
};

class SurfaceThatTakesAdvantageOfNarrowBandAndThereforeMightHaveWrongSign : public Surface {
public:
  virtual bool truncated_distance_may_have_wrong_sign() const override { return true; }
  virtual double point_signed_distance(const stk::math::Vector3d &x) const override
  {
    return truncated_point_signed_distance(x, 0., 0.);
  }
  virtual double truncated_point_signed_distance(const stk::math::Vector3d &x, const double truncation_length, const double far_field_value) const override = 0;
};

class SurfaceThatDoesntTakeAdvantageOfNarrowBandAndThereforeHasCorrectSign : public Surface {
public:
  virtual bool truncated_distance_may_have_wrong_sign() const override { return false; }
  virtual double point_signed_distance(const stk::math::Vector3d &x) const override = 0;
  virtual double truncated_point_signed_distance(const stk::math::Vector3d &x, const double truncation_length, const double far_field_value) const override
  {
    return point_signed_distance(x);
  }
};
} // namespace krino

#endif // Akri_Surface_h
