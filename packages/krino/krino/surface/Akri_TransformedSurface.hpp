// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#ifndef Akri_TransformedSurface_h
#define Akri_TransformedSurface_h

#include <Akri_Surface.hpp>
#include <Akri_Transformation.hpp>

#include <memory>

namespace krino {

// A surface moved by a motion: at a time it is its surface, as given at time 0, where the motion puts
// it then (see Transformation).  Its queries take the query points back to the surface, which does
// not change, so where it is at a time does not depend on the times it was at before, and any surface
// can move.  The time is the one given to prepare_to_compute, or 0 before that.
class TransformedSurface : public Surface {
public:
  TransformedSurface(std::shared_ptr<Surface> surface, const Transformation & motion);
  virtual ~TransformedSurface() {}

  virtual Surface_Type type() const override { return TRANSFORMED_SURFACE; }
  virtual size_t storage_size() const override { return sizeof(TransformedSurface) + mySurface->storage_size(); }
  virtual std::string output() const override;

  virtual void prepare_to_compute(const double time, const BoundingBox & point_bbox, const double truncation_length) override;
  virtual void insert_into(BoundingBox & bbox) const override;
  virtual bool does_intersect(const BoundingBox & bbox) const override;

  virtual double point_signed_distance(const stk::math::Vector3d &x) const override;
  virtual double possibly_truncated_point_signed_distance(const stk::math::Vector3d &x, const double truncation_length, const double far_field_value) const override;
  virtual bool truncated_distance_may_have_wrong_sign() const override { return mySurface->truncated_distance_may_have_wrong_sign(); }

  virtual bool can_approximate_closest_point_normal() const override { return mySurface->can_approximate_closest_point_normal(); }
  virtual stk::math::Vector3d closest_point_normal(const stk::math::Vector3d &x) const override;
  virtual std::pair<int, double> compute_intersection_with_segment(const stk::math::Vector3d &pt0, const stk::math::Vector3d &pt1, const double edgeCrossingTol) const override;
  virtual void fill_segment_crossings(const stk::math::Vector3d & pt0,
      const stk::math::Vector3d & pt1,
      const double edgeCrossingTol,
      std::vector<double> & crossingPositions) const override;

  virtual bool changes_with_time() const override { return true; }
  virtual bool has_corners() const override { return mySurface->has_corners(); }
  virtual bool has_edges() const override { return mySurface->has_edges(); }
  virtual void fill_triangle_intersection_parametric_coordinates(const std::array<stk::math::Vector3d,3> & faceNodes, std::vector<stk::math::Vector3d> & intersectionParamCoords) const override;
  virtual void fill_tetrahedon_intersection_parametric_coordinates(const std::array<stk::math::Vector3d,4> & tetNodes, std::vector<stk::math::Vector3d> & intersectionParamCoords) const override;

  const Surface & get_surface() const { return *mySurface; }

private:
  std::shared_ptr<Surface> mySurface;
  Transformation myMotion;
  RigidTransformation myPlacement;
};

} // namespace krino

#endif // Akri_TransformedSurface_h
