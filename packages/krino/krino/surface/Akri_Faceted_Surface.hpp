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
#include <type_traits>
#include <stk_util/parallel/ParallelComm.hpp>

namespace krino {

template<class FACET>
class Faceted_Surface;

class FacetedSurfaceBase : public SurfaceThatTakesAdvantageOfNarrowBandAndThereforeMightHaveWrongSign
{
public:
  FacetedSurfaceBase(const int dim)
  : SurfaceThatTakesAdvantageOfNarrowBandAndThereforeMightHaveWrongSign(), myDim(dim)
  { STK_ThrowRequire(myDim==2 || myDim==3); }

  static std::unique_ptr<FacetedSurfaceBase> build(const int dim);
  static std::unique_ptr<FacetedSurfaceBase> build_with_velocity(const int dim);

  virtual Surface_Type type() const override { return (myDim == 3) ? FACETED_SURFACE_3D : FACETED_SURFACE_2D; }

  virtual void swap(FacetedSurfaceBase & other) = 0;
  virtual void clear() = 0;
  virtual size_t size() const = 0;
  virtual size_t nonlocal_size() const = 0;
  virtual double point_distance(const stk::math::Vector3d &x, const double narrow_band_size, const double far_field_value, const bool compute_signed_distance) const = 0;
  void prepare_to_compute(const double time,
      const BoundingBox & point_bbox,
      const double truncation_length) override = 0;
  virtual std::string print_sizes() const = 0;
  virtual stk::math::Vector3d closest_point(const stk::math::Vector3d &x) const = 0;

  void prepare_to_compute(const BoundingBox & point_bbox, const double truncation_length)
  {
    STK_ThrowAssert(nullptr == my_transformation);
    prepare_to_compute(0.0, point_bbox, truncation_length);
  }

  double point_unsigned_distance(const stk::math::Vector3d &x, const double narrow_band_size, const double far_field_value) const
  {
    return point_distance(x, narrow_band_size, far_field_value, false);
  }

  template<class FACET>
  const Faceted_Surface<FACET> & as_derived_type() const
  {
    const auto * derived = dynamic_cast<const Faceted_Surface<FACET> *>(this);
    STK_ThrowRequireMsg(derived, "Can't access FacetedSurfaceBase as Faceted_Surface<" << typeid(this).name() << ">");
    return *derived;
  }

  template<class FACET>
  Faceted_Surface<FACET> & as_derived_type()
  {
    auto * derived = dynamic_cast<Faceted_Surface<FACET> *>(this);
    STK_ThrowRequireMsg(derived, "Can't access FacetedSurfaceBase as Faceted_Surface<" << typeid(this).name() << ">");
    return *derived;
  }

  const std::vector<Facet2d> & get_facets_2d() const;
  const std::vector<Facet3d> & get_facets_3d() const;
  std::vector<Facet2d> & get_facets_2d();
  std::vector<Facet3d> & get_facets_3d();

  void emplace_back_2d(const stk::math::Vector3d & x0, const stk::math::Vector3d & x1) { get_facets_2d().emplace_back(x0,x1); }
  void emplace_back_3d(const stk::math::Vector3d & x0, const stk::math::Vector3d & x1, const stk::math::Vector3d & x2)  { get_facets_3d().emplace_back(x0,x1,x2); }

private:
  int myDim;
};

template<class FACET>
class Faceted_Surface : public FacetedSurfaceBase
{
public:
  Faceted_Surface() : FacetedSurfaceBase(FACET::DIM) {}

  virtual size_t storage_size() const override;
  virtual void prepare_to_compute(const double time, const BoundingBox & point_bbox, const double truncation_length) override;
  virtual double truncated_point_signed_distance(const stk::math::Vector3d &x, const double narrow_band_size, const double far_field_value) const override
  {
    return point_distance(x, narrow_band_size, far_field_value, true);
  }
  virtual std::pair<int, double> compute_intersection_with_segment(const stk::math::Vector3d &pt0, const stk::math::Vector3d &pt1, const double edgeCrossingTol) const override;
  
  // query/modify facets
  virtual void clear() override { myLocalFacets.clear(); }
  virtual size_t size() const override { return myLocalFacets.size(); }
  virtual void swap(FacetedSurfaceBase & other) override;

  virtual size_t nonlocal_size() const override { return myNonLocalFacets.size(); }
  std::string print_sizes() const override;

  void parallel_distribute_facets(const size_t batch_size, const std::vector<BoundingBox> & proc_bboxes);
  double point_distance(const stk::math::Vector3d &x, const double narrow_band_size, const double far_field_value, const bool compute_signed_distance) const override;
  
  const std::vector<FACET> & get_facets() const { return myLocalFacets; }
  std::vector<FACET> & get_facets() { return myLocalFacets; }

  virtual stk::math::Vector3d closest_point(const stk::math::Vector3d &x) const override;

public:
  stk::math::Vector3d pseudo_normal_at_closest_point(const stk::math::Vector3d &x) const;
  const FACET * get_closest_facet(const stk::math::Vector3d &x) const;

private:
  virtual void build_local_facets(const BoundingBox & proc_bbox) {}
  
  std::vector<FACET> myLocalFacets;

  mutable std::unique_ptr<SearchTree<const FACET*>> my_facet_tree;
  mutable std::vector<FACET> myNonLocalFacets;
  mutable std::vector<const FACET*> myAllFacetPtrs;
  BoundingBox my_bounding_box;
};

inline std::unique_ptr<FacetedSurfaceBase> FacetedSurfaceBase::build(const int dim)
{
  if (dim == 2)
    return std::make_unique<Faceted_Surface<Facet2d>>();
  return std::make_unique<Faceted_Surface<Facet3d>>();
}

inline std::unique_ptr<FacetedSurfaceBase> FacetedSurfaceBase::build_with_velocity(const int dim)
{
  if (dim == 2)
    return std::make_unique<Faceted_Surface<FacetWithVelocity2d>>();
  return std::make_unique<Faceted_Surface<FacetWithVelocity3d>>();
}

} // namespace krino

#endif // Akri_Faceted_Surface_h
