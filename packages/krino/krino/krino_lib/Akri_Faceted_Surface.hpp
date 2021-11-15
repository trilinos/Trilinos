// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#ifndef Akri_Faceted_Surface_h
#define Akri_Faceted_Surface_h

#include <Akri_Facet.hpp>
#include <Akri_SearchTree.hpp>
#include <Akri_Surface.hpp>

namespace krino {

class Faceted_Surface : public SurfaceThatTakesAdvantageOfNarrowBandAndThereforeMightHaveWrongSign {

public:
  Faceted_Surface(const std::string & sn);

  virtual Surface_Type type() const override { return FACETED_SURFACE; }
  virtual size_t storage_size() const override;
  virtual void pack_into_buffer(stk::CommBuffer & b) const; // pack into buffer for off-processor communication
  virtual void prepare_to_compute(const double time, const BoundingBox & point_bbox, const double truncation_length) override;
  virtual double truncated_point_signed_distance(const Vector3d &x, const double narrow_band_size, const double far_field_value) const override
  {
    return point_distance(x, narrow_band_size, far_field_value, true);
  }
  
  // query/modify facets
  void add( std::unique_ptr<Facet> facet ) { my_local_facets.emplace_back(std::move(facet)); }
  void reserve(unsigned size_) { my_local_facets.reserve(size_); }
  unsigned size() const { return my_local_facets.size(); }
  unsigned nonlocal_size() const { return my_nonlocal_facets.size(); }
  Facet * operator()( const unsigned index ) const { return my_local_facets[index].get(); }
  void clear() { my_local_facets.clear(); }
  void swap(Faceted_Surface & other) { my_local_facets.swap(other.my_local_facets); }
  const FacetOwningVec & get_facets() const { return my_local_facets; }
  void parallel_distribute_facets(const size_t batch_size, const std::vector<BoundingBox> & proc_bboxes);
  
public:
  void prepare_to_compute(const BoundingBox & point_bbox, const double truncation_length)
    { ThrowAssert(nullptr == my_transformation); prepare_to_compute(0.0, point_bbox, truncation_length); }
  double point_unsigned_distance(const Vector3d &x, const double narrow_band_size, const double far_field_value) const
  {
    return point_distance(x, narrow_band_size, far_field_value, false);
  }

private:
  virtual void build_local_facets(const BoundingBox & proc_bbox) {}
  double point_distance(const Vector3d &x, const double narrow_band_size, const double far_field_value, const bool compute_signed_distance) const;
  double compute_point_to_facets_distance_by_average_normal(const Vector3d &x, const FacetVec & facets) const;
  Vector3d compute_pseudo_normal(const unsigned dim, const std::vector<FacetDistanceQuery> & facet_queries, const unsigned nearest) const;
  void gather_nonlocal_facets(const BoundingBox & local_bbox, const double truncation_length);
  
  std::string my_name;
  FacetOwningVec my_local_facets;

  mutable std::unique_ptr<SearchTree<Facet*>> my_facet_tree;
  mutable FacetOwningVec my_nonlocal_facets;
  mutable FacetVec my_all_facets;
  BoundingBox my_bounding_box;
};

} // namespace krino

#endif // Akri_Faceted_Surface_h
