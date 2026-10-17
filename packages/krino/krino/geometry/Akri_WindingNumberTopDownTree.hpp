// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#ifndef KRINO_KRINO_GEOMETRY_AKRI_WINDINGNUMBERTOPDOWNTREE_HPP_
#define KRINO_KRINO_GEOMETRY_AKRI_WINDINGNUMBERTOPDOWNTREE_HPP_

#include <algorithm>
#include <cstddef>
#include <numeric>
#include <vector>

#include <stk_math/StkVector.hpp>
#include <stk_util/util/ReportHandler.hpp>
#include <Akri_WindingNumberCalcs.hpp>
#include <Akri_WindingNumberFacetOps.hpp>
#include <Akri_WindingNumberTreeBase.hpp>

namespace krino {

// Barnes-Hut winding-number tree built top-down: each internal node's multipole
// approximation is accumulated by re-summing ALL facets in its range, so a facet
// is visited once per tree level (O(N log N) build).
//
// Approximations are stored sparsely: internal nodes always carry one, leaves
// only under the approximateLeaves policy. A node with no approximation is
// summed exactly. See WindingNumberBottomUpTree for the dense/O(N)-build variant.
//
// FacetOps is the facet-ops policy (see Akri_WindingNumberFacetOps.hpp): it
// names the facet type and supplies the centroid / bbox / winding / approximation
// operations. See WindingNumberFacetArrayTopDownTree (below) for the plain
// vertex-triple instantiation.
template <typename FacetOps>
class WindingNumberTopDownTree
    : public WindingNumberTreeBase<WindingNumberTopDownTree<FacetOps>, FacetClusterApproximation, FacetOps>
{
  using Base = WindingNumberTreeBase<WindingNumberTopDownTree<FacetOps>, FacetClusterApproximation, FacetOps>;

public:
  using typename Base::Facet;
  using typename Base::FacetVector;
  using typename Base::ApproximationType;

  WindingNumberTopDownTree() = default;

  explicit WindingNumberTopDownTree(
      const FacetVector & facets,
      std::size_t leafSize = 16,
      bool approximateLeaves = false,
      const FacetOps & ops = detail::default_facet_ops<FacetOps>())
  {
    build(facets, leafSize, approximateLeaves, ops);
  }

  // The tree stores a pointer to facets, not a copy.
  // The caller must ensure facets outlive this tree.
  //
  // approximateLeaves controls the leaf policy for far-away leaves that pass
  // the Barnes-Hut acceptance test:
  //   false (default): always evaluate leaf facets exactly (more accurate).
  //   true:            build a multipole approximation for leaves too, so an
  //                    accepted leaf uses it (fewer atan2/sqrt for large leaves,
  //                    at the cost of storing an approximation per leaf).
  //
  // ops is the facet-ops policy; the tree keeps only a POINTER to it (like
  // facets), so the passed instance must outlive the tree. The default is a
  // static-storage stateless instance, which suffices for stateless policies. A
  // stateful policy (e.g. one holding a pointer to the facet source) must be
  // passed its own lvalue here so the tree resolves facets through it.
  void build(
      const FacetVector & facets,
      std::size_t leafSize = 16,
      bool approximateLeaves = false,
      const FacetOps & ops = detail::default_facet_ops<FacetOps>());

private:
  // CRTP hook used by WindingNumberTreeBase's query path. Returns this node's
  // approximation, or nullptr when the node has none (sparse storage: exact-only
  // leaves keep the NO_APPROXIMATION sentinel).
  friend Base;

  const ApproximationType * approximation_for_node(std::size_t nodeIndex) const
  {
    const std::size_t approximationIndex = myApproximationIndices[nodeIndex];
    if (approximationIndex == Base::NO_APPROXIMATION)
    {
      return nullptr;
    }
    return &myApproximations[approximationIndex];
  }

  std::size_t build_node(
      std::size_t begin,
      std::size_t end);

private:
  // Pull the base-class members into scope so the templated method bodies can
  // name them unqualified (they are dependent names in a template base).
  using Base::myFacets;
  using Base::myOps;
  using Base::myFacetIndices;
  using Base::myNodes;
  using Base::myLeafSize;
  using Node = typename Base::Node;

  // Sparse node -> approximation mapping, indexed 1:1 with myNodes. Entries are
  // NO_APPROXIMATION for exact-only leaves; otherwise an index into
  // myApproximations. Kept out of Node to keep traversal cache-friendly and to
  // avoid storing an unused approximation on exact-only leaves.
  std::vector<std::size_t> myApproximationIndices;
  std::vector<ApproximationType> myApproximations;

  // See build(): whether leaves carry a multipole approximation.
  bool myApproximateLeaves = false;
};

template <typename FacetOps>
void WindingNumberTopDownTree<FacetOps>::build(
    const FacetVector & facets,
    const std::size_t leafSize,
    const bool approximateLeaves,
    const FacetOps & ops)
{
  myFacets = &facets;
  myOps = &ops;
  myLeafSize = std::max<std::size_t>(std::size_t{1}, leafSize);
  myApproximateLeaves = approximateLeaves;

  myFacetIndices.resize(facets.size());
  std::iota(myFacetIndices.begin(), myFacetIndices.end(), std::size_t{0});

  myNodes.clear();
  myApproximationIndices.clear();
  myApproximations.clear();

  if (facets.empty())
  {
    return;
  }

  const std::size_t estimatedLeaves =
      (facets.size() + myLeafSize - 1) / myLeafSize;

  const std::size_t estimatedNodes =
      std::max<std::size_t>(std::size_t{1}, 2 * estimatedLeaves - 1);

  myNodes.reserve(estimatedNodes);
  myApproximationIndices.reserve(estimatedNodes);

  const std::size_t rootIndex = build_node(0, facets.size());
  STK_ThrowRequire(rootIndex == 0);
  STK_ThrowAssert(myNodes.size() == myApproximationIndices.size());
}

template <typename FacetOps>
std::size_t WindingNumberTopDownTree<FacetOps>::build_node(
    const std::size_t begin,
    const std::size_t end)
{
  STK_ThrowAssert(myFacets != nullptr);
  STK_ThrowAssert(begin < end);
  STK_ThrowAssert(end <= myFacetIndices.size());

  const FacetVector & facets = *myFacets;

  // Reserve this node's slot (preorder index) and its 1:1 approximation-index
  // slot. Both vectors may reallocate on deeper emplace_back calls during
  // recursion, so we re-index into them after any recursion below.
  const std::size_t nodeIndex = myNodes.size();
  myNodes.emplace_back();
  myApproximationIndices.emplace_back(Base::NO_APPROXIMATION);

  {
    Node & node = myNodes[nodeIndex];
    node.begin = begin;
    node.end = end;
    node.bbox = detail::compute_range_aabb(facets, myFacetIndices, begin, end, *myOps);
    node.radius = detail::aabb_radius(node.bbox);
  }

  const stk::math::Vector3d center = myNodes[nodeIndex].bbox.center();

  const std::size_t count = end - begin;
  const bool isLeaf = count <= myLeafSize;

  // Build a multipole approximation when the node is internal (always needed),
  // or when the leaf policy asks leaves to carry one. Exact-only leaves store
  // nothing, keeping approximation storage and traversal lean.
  if (!isLeaf || myApproximateLeaves)
  {
    const std::size_t approximationIndex = myApproximations.size();
    myApproximations.emplace_back();
    myApproximationIndices[nodeIndex] = approximationIndex;

    ApproximationType & approximation = myApproximations[approximationIndex];
    approximation.clear();
    approximation.center = center;

    for (std::size_t i = begin; i < end; ++i)
    {
      const Facet & facet = facets[myFacetIndices[i]];

      myOps->accumulate_approximation(facet, center, approximation);
    }
  }

  if (isLeaf)
  {
    return nodeIndex;
  }

  const unsigned axis = myNodes[nodeIndex].bbox.max_span_direction();
  const std::size_t mid = begin + count / 2;

  std::nth_element(
      myFacetIndices.begin() + begin,
      myFacetIndices.begin() + mid,
      myFacetIndices.begin() + end,
      [&](const std::size_t lhs, const std::size_t rhs)
      {
        const stk::math::Vector3d lhsCentroid = myOps->get_centroid(facets[lhs]);
        const stk::math::Vector3d rhsCentroid = myOps->get_centroid(facets[rhs]);

        return lhsCentroid[axis] < rhsCentroid[axis];
      });

  const std::size_t leftIndex = build_node(begin, mid);
  STK_ThrowRequire(leftIndex == nodeIndex + 1);

  const std::size_t rightIndex = build_node(mid, end);

  // Do this after recursion. References into myNodes may have been invalidated
  // by emplace_back during recursive construction.
  myNodes[nodeIndex].rightIndex = rightIndex;

  return nodeIndex;
}

// Top-down tree over plain vertex-triple facets (std::array<Vector3d,3>).
using WindingNumberFacetArrayTopDownTree = WindingNumberTopDownTree<WindingNumberArrayFacetOps>;

}

#endif /* KRINO_KRINO_GEOMETRY_AKRI_WINDINGNUMBERTOPDOWNTREE_HPP_ */
