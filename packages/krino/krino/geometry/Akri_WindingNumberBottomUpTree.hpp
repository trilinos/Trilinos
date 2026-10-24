// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#ifndef KRINO_KRINO_GEOMETRY_AKRI_WINDINGNUMBERBOTTOMUPTREE_HPP_
#define KRINO_KRINO_GEOMETRY_AKRI_WINDINGNUMBERBOTTOMUPTREE_HPP_

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

// Barnes-Hut winding-number tree that stores a ClusterApproximation for EVERY
// node (leaves included) and builds it bottom-up.
//
// This differs from WindingNumberTopDownTree in two related ways:
//
//   * Storage: every node carries a multipole approximation, held 1:1 in a
//     parallel vector (node i <-> myApproximations[i]). There is no
//     approximation-index or sentinel, and no "exact leaves" policy: a leaf can
//     always be approximated when it passes the acceptance test.
//
//   * Build cost: WindingNumberTopDownTree assembles each internal node by
//     re-summing ALL facets in its range (O(N log N) build). This tree instead
//     computes each leaf's moments once from its facets, then forms every
//     internal node's moments as the O(1) parallel-axis merge of its two
//     children, for an O(N) build overall. Because the evaluated expansion is
//     identical, query cost and accuracy match WindingNumberTopDownTree; only
//     build time and per-node storage differ (ClusterApproximation is ~24%
//     smaller than FacetClusterApproximation).
//
// FacetOps is the facet-ops policy (see Akri_WindingNumberFacetOps.hpp): it
// names the facet type and supplies the centroid / bbox / winding / approximation
// operations. See WindingNumberFacetArrayBottomUpTree (below) for the plain
// vertex-triple instantiation.
template <typename FacetOps>
class WindingNumberBottomUpTree
    : public WindingNumberTreeBase<WindingNumberBottomUpTree<FacetOps>, ClusterApproximation, FacetOps>
{
  using Base = WindingNumberTreeBase<WindingNumberBottomUpTree<FacetOps>, ClusterApproximation, FacetOps>;

public:
  using typename Base::Facet;
  using typename Base::FacetVector;
  using typename Base::ApproximationType;

  WindingNumberBottomUpTree() = default;

  explicit WindingNumberBottomUpTree(
      const FacetVector & facets,
      std::size_t leafSize = 16,
      const FacetOps & ops = detail::default_facet_ops<FacetOps>())
  {
    build(facets, leafSize, ops);
  }

  // The tree stores a pointer to facets, not a copy.
  // The caller must ensure facets outlive this tree.
  //
  // ops is the facet-ops policy; the tree keeps only a POINTER to it (like
  // facets), so the passed instance must outlive the tree. The default is a
  // static-storage stateless instance, which suffices for stateless policies. A
  // stateful policy (e.g. one holding a pointer to the facet source) must be
  // passed its own lvalue here so the tree resolves facets through it.
  void build(
      const FacetVector & facets,
      std::size_t leafSize = 16,
      const FacetOps & ops = detail::default_facet_ops<FacetOps>());

private:
  // CRTP hook used by WindingNumberTreeBase's query path. Every node carries an
  // approximation, indexed 1:1 with myNodes, so this never returns nullptr.
  friend Base;

  const ApproximationType * approximation_for_node(std::size_t nodeIndex) const
  {
    return &myApproximations[nodeIndex];
  }

  // Recursively partitions [begin, end) and, on the way back up, fills in this
  // node's bbox/radius/approximation. Returns the node's index. The
  // approximation is assembled from children (internal) or facets (leaf), so
  // this must run post-order.
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

  // Multipole approximations, one per node and indexed 1:1 with myNodes (node i
  // uses myApproximations[i]). Kept in a parallel array rather than inside Node
  // to keep traversal (bbox/radius/children) cache-friendly.
  std::vector<ApproximationType> myApproximations;
};

template <typename FacetOps>
void WindingNumberBottomUpTree<FacetOps>::build(
    const FacetVector & facets,
    const std::size_t leafSize,
    const FacetOps & ops)
{
  myFacets = &facets;
  myOps = &ops;
  myLeafSize = std::max<std::size_t>(std::size_t{1}, leafSize);

  myFacetIndices.resize(facets.size());
  std::iota(myFacetIndices.begin(), myFacetIndices.end(), std::size_t{0});

  myNodes.clear();
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
  myApproximations.reserve(estimatedNodes);

  const std::size_t rootIndex = build_node(0, facets.size());
  STK_ThrowRequire(rootIndex == 0);
  STK_ThrowAssert(myNodes.size() == myApproximations.size());
}

template <typename FacetOps>
std::size_t WindingNumberBottomUpTree<FacetOps>::build_node(
    const std::size_t begin,
    const std::size_t end)
{
  STK_ThrowAssert(myFacets != nullptr);
  STK_ThrowAssert(begin < end);
  STK_ThrowAssert(end <= myFacetIndices.size());

  const FacetVector & facets = *myFacets;

  // Reserve this node's slot (preorder index) along with its 1:1 approximation
  // slot. Both vectors may reallocate on deeper emplace_back calls during
  // recursion, so we re-index into them after any recursion below rather than
  // holding references across it.
  const std::size_t nodeIndex = myNodes.size();
  myNodes.emplace_back();
  myApproximations.emplace_back();

  {
    Node & node = myNodes[nodeIndex];
    node.begin = begin;
    node.end = end;
    node.bbox = detail::compute_range_aabb(facets, myFacetIndices, begin, end, *myOps);
    node.radius = detail::aabb_radius(node.bbox);
  }

  // The expansion point is the bbox center, matching should_use_approximation.
  const stk::math::Vector3d center = myNodes[nodeIndex].bbox.center();

  const std::size_t count = end - begin;
  const bool isLeaf = count <= myLeafSize;

  if (isLeaf)
  {
    // Leaf: assemble the approximation directly from its facets. Each facet's
    // moments are formed about the facet centroid, then shifted to the leaf
    // center by the parallel-axis merge.
    ApproximationType & approximation = myApproximations[nodeIndex];
    approximation.clear();
    approximation.center = center;

    for (std::size_t i = begin; i < end; ++i)
    {
      const Facet & facet = facets[myFacetIndices[i]];

      ApproximationType facetApproximation;
      myOps->compute_facet_approximation(facet, facetApproximation);
      accumulate_cluster_approximation(facetApproximation, approximation);
    }

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

  // Do this after recursion. References into myNodes/myApproximations may have
  // been invalidated by emplace_back during recursive construction.
  myNodes[nodeIndex].rightIndex = rightIndex;

  // Bottom-up assembly: this internal node's moments are the parallel-axis
  // merge of its two children's, an O(1) combine that avoids re-summing facets.
  ApproximationType & approximation = myApproximations[nodeIndex];
  approximation.clear();
  approximation.center = center;
  accumulate_cluster_approximation(myApproximations[leftIndex], approximation);
  accumulate_cluster_approximation(myApproximations[rightIndex], approximation);

  return nodeIndex;
}

// Bottom-up tree over plain vertex-triple facets (std::array<Vector3d,3>).
using WindingNumberFacetArrayBottomUpTree = WindingNumberBottomUpTree<WindingNumberArrayFacetOps>;

}

#endif /* KRINO_KRINO_GEOMETRY_AKRI_WINDINGNUMBERBOTTOMUPTREE_HPP_ */
