// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#ifndef Akri_Composite_Surface_h
#define Akri_Composite_Surface_h

#include <Akri_Surface.hpp>

namespace krino {

class SurfaceTree;
class Transformation;

class Composite_Surface : public SurfaceThatTakesAdvantageOfNarrowBandAndThereforeMightHaveWrongSign {

public:
  Composite_Surface(const std::string & sn);
  virtual ~Composite_Surface();
  
  enum CompositionMethod{ MINIMUM_SIGNED_DISTANCE, MAXIMUM_SIGNED_DISTANCE };

  virtual Surface_Type type() const override { return COMPOSITE_SURFACE; }
  virtual void prepare_to_compute(const double time, const BoundingBox & point_bbox, const double truncation_length) override;
  virtual size_t storage_size() const override { size_t tot_size = 0; for (auto && surf : my_subsurfaces) tot_size += surf->storage_size(); return tot_size; }
  
  // query/modify subsurfaces
  void add( Surface * surf ) { my_subsurfaces.emplace_back(surf); }
  void reserve(unsigned size_) { my_subsurfaces.reserve(size_); }
  unsigned size() const { return my_subsurfaces.size(); }
  Surface * operator()( const unsigned index ) const { return my_subsurfaces[index].get(); }
  void clear() { my_subsurfaces.clear(); }
  void swap(Composite_Surface & other) { my_subsurfaces.swap(other.my_subsurfaces); }
  SurfaceAutoVec & get_surfaces() { return my_subsurfaces; }
  const SurfaceAutoVec & get_surfaces() const { return my_subsurfaces; }
  
  virtual void set_transformation(Transformation * trans) override { Surface::set_transformation(trans); for (auto && subsurf : my_subsurfaces) subsurf->set_transformation(trans);}
  void set_composition_method(const CompositionMethod composition_method) { my_composition_method = composition_method; }
  void insert_into(BoundingBox & bbox) const override;
  bool does_intersect(const BoundingBox & bbox) const override;
  virtual double truncated_point_signed_distance(const stk::math::Vector3d &x, const double truncation_length, const double far_field_value) const override;
  double point_signed_distance_with_narrow_band(const stk::math::Vector3d &x, const double narrow_band) const
  {
    // Special to this class, this version of point_signed_distance automatically determines the far_field_value from the composition_method
    return truncated_point_signed_distance(x, narrow_band, get_signed_narrow_band_size(narrow_band));
  }

  virtual bool truncated_distance_may_have_wrong_sign() const override { return my_subsurfaces_might_produce_wrong_sign; }
  double get_signed_narrow_band_size(const double pos_narrow_band) const { return (my_composition_method == MINIMUM_SIGNED_DISTANCE) ? pos_narrow_band : -pos_narrow_band; }
  
protected:
  std::string my_name;
  CompositionMethod my_composition_method;
  SurfaceAutoVec my_subsurfaces;
  bool my_subsurfaces_might_produce_wrong_sign = true;
};

} // namespace krino


#endif // Akri_Composite_Surface_h
