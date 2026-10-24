// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.
#ifndef KRINO_KRINO_SURFACE_AKRI_CSG_SURFACE_HPP_
#define KRINO_KRINO_SURFACE_AKRI_CSG_SURFACE_HPP_

#include <Akri_Surface.hpp>
#include <Akri_CSG_Tree.hpp>

namespace krino {

class CSG_Surface : public Surface
{
public:
  CSG_Surface() = default;
  CSG_Surface(const CSGTree & tree);
  virtual ~CSG_Surface();

  virtual Surface_Type type() const override { return CSG_SURFACE; }
  virtual std::string output() const override;
  virtual void prepare_to_compute(const double time, const BoundingBox & point_bbox, const double truncation_length) override;
  virtual size_t storage_size() const override;

  void set_tree(const CSGTree & tree);
  const std::vector<std::shared_ptr<Surface>> & get_subsurfaces() const { return myAllSurfaces; }
  void set_inside_to_be_positive() { myIsInsideNegative = false; }
  bool is_inside_negative() const { return myIsInsideNegative; }
  virtual void set_transformation(Transformation * trans) override;
  void insert_into(BoundingBox & bbox) const override;
  bool does_intersect(const BoundingBox & bbox) const override;
  virtual double possibly_truncated_point_signed_distance(const stk::math::Vector3d &x, const double truncationLength, const double farFieldValue) const override;
  virtual bool truncated_distance_may_have_wrong_sign() const override { return my_subsurfaces_might_produce_wrong_sign; }
  void clear();
  bool empty() const;
  double get_signed_narrow_band_size(const double positiveNarrowBand) const { return myIsInsideNegative ? positiveNarrowBand : -positiveNarrowBand; }
  virtual double point_signed_distance(const stk::math::Vector3d &x) const override
  {
    return possibly_truncated_point_signed_distance(x, 0., 0.);
  }
  double point_signed_distance_with_narrow_band(const stk::math::Vector3d &x, const double positiveNarrowBand) const
  {
    // Special to this class, this version of point_signed_distance automatically determines the far_field_value from myIsInsideNegative
    return possibly_truncated_point_signed_distance(x, positiveNarrowBand, get_signed_narrow_band_size(positiveNarrowBand));
  }

private:
  bool myIsInsideNegative{true};
  CSGTree myAllTree;
  CSGTree myNearbyTree;
  std::vector<std::shared_ptr<Surface>> myAllSurfaces;
  std::vector<std::shared_ptr<Surface>> myNearbySurfaces;
  bool my_subsurfaces_might_produce_wrong_sign = true;
};

} // namespace krino

#endif /* KRINO_KRINO_SURFACE_AKRI_CSG_SURFACE_HPP_ */
