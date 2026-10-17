// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#ifndef KRINO_KRINO_GEOMETRY_AKRI_WINDINGNUMBERFACETOPS_HPP_
#define KRINO_KRINO_GEOMETRY_AKRI_WINDINGNUMBERFACETOPS_HPP_

#include <array>

#include <stk_math/StkVector.hpp>
#include <Akri_BoundingBox.hpp>
#include <Akri_WindingNumberCalcs.hpp>

namespace krino {

// Facet-ops policy for the Barnes-Hut winding-number trees.
//
// The trees (WindingNumberTreeBase and its concrete subclasses) are templated
// on a policy type that names the facet representation and supplies every
// facet-touching operation, so a tree can be built over any facet type without
// the tree ever reaching into facet internals. This mirrors the existing krino
// convention in Akri_SearchTree.hpp (templated on get_centroid / insert_into_bbox
// callables, fed by Facet3d::get_centroid / Facet3d::insert_into_bounding_box).
//
// The tree holds a pointer to a policy instance (see WindingNumberTreeBase::
// myOps), so the operations are const member functions rather than static: a
// policy may carry state shared by all facets in the tree (e.g. a pointer to the
// facet source). Because the tree only points at the policy, the policy instance
// must outlive the tree.
//
// A policy must provide:
//
//   using FacetType = ...;   // the stored facet type
//
//   stk::math::Vector3d
//     get_centroid(const FacetType & facet) const;
//       // facet centroid, used for the object-median split
//
//   void
//     insert_into_bbox(const FacetType & facet, BoundingBox_T<double,3> & bbox) const;
//       // expand bbox to contain the facet (note: double-precision box)
//
//   double
//     compute_exact_winding_number(const FacetType & facet,
//                                  const stk::math::Vector3d & queryPt) const;
//       // exact per-facet winding-number contribution
//
//   void
//     accumulate_approximation(const FacetType & facet,
//                              const stk::math::Vector3d & center,
//                              FacetClusterApproximation & approx) const;
//       // accumulate this facet's moments into a top-down cluster about center
//
//   void
//     compute_facet_approximation(const FacetType & facet,
//                                 ClusterApproximation & approx) const;
//       // build a single-facet cluster (the bottom-up leaf primitive)
//
// The child->parent merge accumulate_cluster_approximation(child, parent)
// operates on ClusterApproximation objects, not facets, so it is NOT part of
// this policy; the bottom-up build calls it directly.

// CRTP mixin that implements the full WindingNumberFacetOps contract for any
// facet type that can be reduced to three vertices. A derived policy supplies
// only:
//
//   using FacetType = <your facet>;                    // (== the Facet param)
//   std::array<stk::math::Vector3d,3>
//     get_vertices(const FacetType & facet) const;     // the one required method
//
// Every op converts the facet to a vertex triple via derived().get_vertices and
// delegates to the free functions in Akri_WindingNumberCalcs. Facet is a
// template parameter (not Derived::FacetType) so the method signatures don't
// require Derived to be complete when it inherits.
//
// get_vertices is a const member (not static), so a derived policy may resolve
// the vertices through state it holds (e.g. a pointer to the facet source),
// letting the facet itself be a bare handle such as an index.
//
// e.g. to build a tree over krino::Facet3d:
//
//   struct Facet3dOps
//     : WindingNumberVertexFacetOps<Facet3dOps, Facet3d>
//   {
//     std::array<stk::math::Vector3d,3> get_vertices(const Facet3d & f) const
//     { return f.facet_vertices(); }
//   };
template <typename Derived, typename Facet>
struct WindingNumberVertexFacetOps
{
  using FacetType = Facet;

  stk::math::Vector3d get_centroid(const FacetType & facet) const
  {
    const auto & v = derived().get_vertices(facet);
    return (v[0] + v[1] + v[2]) / 3.0;
  }

  void insert_into_bbox(const FacetType & facet, BoundingBox_T<double, 3> & bbox) const
  {
    const auto & v = derived().get_vertices(facet);
    bbox.accommodate(v[0]);
    bbox.accommodate(v[1]);
    bbox.accommodate(v[2]);
  }

  double compute_exact_winding_number(const FacetType & facet, const stk::math::Vector3d & queryPt) const
  {
    const auto & v = derived().get_vertices(facet);
    return compute_exact_facet_winding_number(v[0], v[1], v[2], queryPt);
  }

  void accumulate_approximation(const FacetType & facet, const stk::math::Vector3d & center, FacetClusterApproximation & approx) const
  {
    const auto & v = derived().get_vertices(facet);
    accumulate_cluster_approximation(v[0], v[1], v[2], center, approx);
  }

  void compute_facet_approximation(const FacetType & facet, ClusterApproximation & approx) const
  {
    krino::compute_facet_approximation(derived().get_vertices(facet), approx);
  }

private:
  const Derived & derived() const
  {
    return static_cast<const Derived &>(*this);
  }
};

// Default policy: facets are already plain vertex triples, so get_vertices is
// the identity. Behavior on this path is identical to the pre-generalization
// trees.
struct WindingNumberArrayFacetOps
    : WindingNumberVertexFacetOps<WindingNumberArrayFacetOps, std::array<stk::math::Vector3d, 3>>
{
  const std::array<stk::math::Vector3d, 3> & get_vertices(const FacetType & facet) const
  {
    return facet;
  }
};

}

#endif /* KRINO_KRINO_GEOMETRY_AKRI_WINDINGNUMBERFACETOPS_HPP_ */
