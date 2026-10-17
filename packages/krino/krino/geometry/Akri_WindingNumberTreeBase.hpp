// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#ifndef KRINO_KRINO_GEOMETRY_AKRI_WINDINGNUMBERTREEBASE_HPP_
#define KRINO_KRINO_GEOMETRY_AKRI_WINDINGNUMBERTREEBASE_HPP_

#include <array>
#include <cstddef>
#include <limits>
#include <vector>

#include <stk_math/StkVector.hpp>
#include <stk_util/util/ReportHandler.hpp>
#include <Akri_BoundingBox.hpp>
#include <Akri_WindingNumberCalcs.hpp>
#include <Akri_WindingNumberFacetOps.hpp>

namespace krino {

namespace detail {

// Shared tree-building helpers. These are pure geometry and are used by the
// build routines of every concrete winding-number tree.

// Distance from the AABB center to its farthest corner. For an axis-aligned
// box this is exactly half the diagonal length, so it bounds every facet in
// the box. It is a conservative (upper) bound on the multipole radius.
inline double aabb_radius(
    const BoundingBox_T<double, 3> & bbox)
{
  return 0.5 * std::sqrt(bbox.size_squared());
}

template <typename FacetOps, typename FacetVector>
BoundingBox_T<double, 3> compute_range_aabb(
    const FacetVector & facets,
    const std::vector<std::size_t> & facetIndices,
    const std::size_t begin,
    const std::size_t end,
    const FacetOps & ops)
{
  BoundingBox_T<double, 3> bbox;

  for (std::size_t i = begin; i < end; ++i)
  {
    ops.insert_into_bbox(facets[facetIndices[i]], bbox);
  }

  return bbox;
}

// Static-storage default for a stateless facet-ops policy. The trees keep only a
// pointer to their ops (see WindingNumberTreeBase::myOps), so a defaulted
// temporary would dangle; this returns a reference to an object with static
// lifetime that outlives any tree. Only used as the build()/ctor default (i.e.
// the stateless array policy); a stateful policy is always passed its own lvalue.
template <typename FacetOps>
const FacetOps & default_facet_ops()
{
  static const FacetOps ops{};
  return ops;
}

} // namespace detail

// Common base for the Barnes-Hut winding-number trees. It owns the shared value
// types, the tree scaffolding (facet permutation + node array), and the entire
// query path. It does NOT own the approximation storage or the build strategy:
// those differ between the concrete trees and stay in the derived classes.
//
// This is a CRTP base (no virtual dispatch): the derived class is passed as
// Derived and supplies exactly one hook,
//
//     const Approximation * approximation_for_node(std::size_t nodeIndex) const
//
// returning the node's multipole approximation, or nullptr when the node has
// none and must be evaluated exactly. That single hook is what lets one query
// implementation serve both the sparse (top-down) and dense (bottom-up)
// approximation-storage schemes.
template <typename Derived, typename Approximation, typename FacetOps>
class WindingNumberTreeBase
{
public:
  using FacetOpsType = FacetOps;
  using Facet = typename FacetOps::FacetType;
  using FacetVector = std::vector<Facet>;
  using ApproximationType = Approximation;

  // Double-precision bounding box; the box center is used as the multipole
  // expansion point, so we keep double precision here rather than using the
  // float BoundingBox typedef.
  using Aabb = BoundingBox_T<double, 3>;

  struct QueryStats
  {
    std::size_t nodesVisited = 0;
    std::size_t approximateNodesUsed = 0;
    std::size_t leavesVisited = 0;
    std::size_t exactFacetsEvaluated = 0;

    void clear()
    {
      nodesVisited = 0;
      approximateNodesUsed = 0;
      leavesVisited = 0;
      exactFacetsEvaluated = 0;
    }
  };

  struct Node
  {
    Aabb bbox;

    double radius = 0.0;

    // Range into myFacetIndices.
    std::size_t begin = 0;
    std::size_t end = 0;

    // Preorder binary tree layout:
    //   left child  = this node index + 1
    //   right child = rightIndex
    //
    // rightIndex == 0 means leaf. The root is index 0 and is never a right child,
    // so zero is an unambiguous sentinel.
    std::size_t rightIndex = 0;

    bool is_leaf() const
    {
      return rightIndex == 0;
    }

    bool have_children() const
    {
      return rightIndex != 0;
    }

    std::size_t num_facets() const
    {
      return end - begin;
    }
  };

  bool empty() const
  {
    return myNodes.empty();
  }

  std::size_t num_nodes() const
  {
    return myNodes.size();
  }

  std::size_t num_facets() const
  {
    return myFacets ? myFacets->size() : 0;
  }

  std::size_t leaf_size() const
  {
    return myLeafSize;
  }

  const std::vector<Node> & nodes() const
  {
    return myNodes;
  }

  // Approximate winding number.
  //
  // theta is the Barnes-Hut acceptance ratio, NOT an angle: a node's cached
  // multipole approximation is used when radius / distance < theta, where
  // radius is the node's bounding radius and distance is from the query point
  // to the node center. It is dimensionless and typically in (0, 1].
  //
  // Smaller theta:
  //   stricter acceptance, more nodes expanded -> more accurate, slower
  //
  // Larger theta:
  //   looser acceptance, more nodes approximated -> faster, less accurate
  //
  // theta <= 0 disables approximation entirely (exact evaluation everywhere).
  double compute_winding_number(
      const stk::math::Vector3d & queryPt,
      double theta = 0.5,
      QueryStats * stats = nullptr) const
  {
    if (stats != nullptr)
    {
      stats->clear();
    }

    if (myNodes.empty())
    {
      return 0.0;
    }

    return compute_winding_number_from_node(0, queryPt, theta, stats);
  }

  // Always computes the exact winding number by summing all facets.
  // Useful for testing the tree approximation.
  double compute_exact_winding_number(
      const stk::math::Vector3d & queryPt) const
  {
    if (myFacets == nullptr)
    {
      return 0.0;
    }

    double windingNumber = 0.0;

    for (const Facet & facet : *myFacets)
    {
      windingNumber += myOps->compute_exact_winding_number(facet, queryPt);
    }

    return windingNumber;
  }

protected:
  // Sentinel for "this node has no approximation" in a sparse node ->
  // approximation index mapping (used by the top-down tree).
  static constexpr std::size_t NO_APPROXIMATION =
      std::numeric_limits<std::size_t>::max();

  const Derived & derived() const
  {
    return static_cast<const Derived &>(*this);
  }

  double compute_winding_number_from_node(
      const std::size_t nodeIndex,
      const stk::math::Vector3d & queryPt,
      const double theta,
      QueryStats * stats) const
  {
    STK_ThrowAssert(myFacets != nullptr);
    STK_ThrowAssert(nodeIndex < myNodes.size());

    if (stats != nullptr)
    {
      ++stats->nodesVisited;
    }

    const Node & node = myNodes[nodeIndex];

    // The derived class decides whether this node carries an approximation; a
    // node without one (nullptr) is either summed exactly (leaf) or descended.
    const Approximation * approximation = derived().approximation_for_node(nodeIndex);

    if (approximation != nullptr && should_use_approximation(node, queryPt, theta))
    {
      if (stats != nullptr)
      {
        ++stats->approximateNodesUsed;
      }

      return compute_approximate_winding_number(*approximation, queryPt);
    }

    // Not approximated: a leaf is summed exactly, an internal node descends.
    if (node.is_leaf())
    {
      if (stats != nullptr)
      {
        ++stats->leavesVisited;
        stats->exactFacetsEvaluated += node.num_facets();
      }

      const FacetVector & facets = *myFacets;

      double windingNumber = 0.0;

      for (std::size_t i = node.begin; i < node.end; ++i)
      {
        const Facet & facet = facets[myFacetIndices[i]];

        windingNumber += myOps->compute_exact_winding_number(facet, queryPt);
      }

      return windingNumber;
    }

    const std::size_t leftIndex = nodeIndex + 1;
    const std::size_t rightIndex = node.rightIndex;

    return compute_winding_number_from_node(leftIndex, queryPt, theta, stats)
         + compute_winding_number_from_node(rightIndex, queryPt, theta, stats);
  }

  bool should_use_approximation(
      const Node & node,
      const stk::math::Vector3d & queryPt,
      const double theta) const
  {
    if (theta <= 0.0)
    {
      return false;
    }

    // The approximation's expansion point is the bbox center by construction, so
    // the Barnes-Hut acceptance test needs only the (small) node, not the
    // approximation object itself.
    const stk::math::Vector3d center = node.bbox.center();

    const double d2 = (queryPt-center).length_squared();

    if (d2 <= 0.0)
    {
      return false;
    }

    // Use squared form of:
    //
    //     radius / distance < theta
    //
    // to avoid a sqrt during traversal.
    const double theta2 = theta * theta;
    const double r2 = node.radius * node.radius;

    return r2 < theta2 * d2;
  }

protected:
  const FacetVector * myFacets = nullptr;

  // The facet-ops policy, held by pointer: the tree does not own it, so the
  // pointed-to instance must outlive the tree (mirroring myFacets). Set by
  // build() from its ops argument. A stateful policy carries data shared by all
  // facets in the tree, such as a pointer to the facet source, held once here
  // instead of per facet.
  const FacetOps * myOps = nullptr;

  // Permutation of facet indices. Each node owns a contiguous subrange of this.
  std::vector<std::size_t> myFacetIndices;

  std::vector<Node> myNodes;

  std::size_t myLeafSize = 16;
};

}

#endif /* KRINO_KRINO_GEOMETRY_AKRI_WINDINGNUMBERTREEBASE_HPP_ */
