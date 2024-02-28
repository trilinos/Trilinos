// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#ifndef KRINO_INCLUDE_AKRI_INTERSECTION_POINTS_H_
#define KRINO_INCLUDE_AKRI_INTERSECTION_POINTS_H_

#include <Akri_InterfaceID.hpp>
#include <Akri_MasterElement.hpp>
#include <Akri_NodeToCapturedDomains.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <vector>

namespace krino {

class CDFEM_Support;
class Phase_Support;
class CDFEM_Parent_Edge;
class InterfaceGeometry;

class IntersectionPoint
{
public:
  IntersectionPoint(const bool owned, const std::vector<stk::mesh::Entity> & nodes, const std::vector<double> & weights, const std::vector<int> & sortedDomains)
  : mOwned(owned), mNodes(nodes), mWeights(weights), mSortedDomains(sortedDomains) {STK_ThrowAssert(mNodes.size() == mWeights.size());}
  bool is_owned() const { return mOwned; }
  const std::vector<stk::mesh::Entity> & get_nodes() const { return mNodes; }
  const std::vector<double> & get_weights() const { return mWeights; }
  const std::vector<int> & get_sorted_domains() const { return mSortedDomains; }
private:
  bool mOwned;
  std::vector<stk::mesh::Entity> mNodes;
  std::vector<double> mWeights;
  std::vector<int> mSortedDomains;
};

class InterpolationPoint
{
public:
  InterpolationPoint(const std::vector<stk::mesh::Entity> & nodes, const std::vector<double> & weights)
  : mNodes(nodes), mWeights(weights) {STK_ThrowAssert(mNodes.size() == mWeights.size());}
  const std::vector<stk::mesh::Entity> & get_nodes() const { return mNodes; }
  const std::vector<double> & get_weights() const { return mWeights; }
private:
  std::vector<stk::mesh::Entity> mNodes;
  std::vector<double> mWeights;
};

struct EdgeIntersection
{
  EdgeIntersection(const IntersectionPoint & intersectionPt);
  std::array<stk::mesh::Entity,2> nodes;
  double crossingLocation;
  InterfaceID interface;
};

template <typename T>
bool any_entity_in_first_vector_is_contained_in_second_sorted_vector(const std::vector<T> & firstVec, const std::vector<T> & secondVec)
{
  for (auto && first : firstVec)
    if (std::binary_search(secondVec.begin(), secondVec.end(), first))
      return true;
  return false;
}

typedef std::function<bool(const std::vector<stk::mesh::Entity> &, const std::vector<int> &)> IntersectionPointFilter;

bool first_sorted_vector_of_domains_contains_all_domains_in_second_vector(const std::vector<int> & firstVec, const std::vector<int> & secondVec);
bool domains_already_snapped_to_node_are_also_at_intersection_point(const NodeToCapturedDomainsMap & nodesToCapturedDomains, stk::mesh::Entity node, const std::vector<int> & intersectionPointDomains);
std::string debug_output(const stk::mesh::BulkData & mesh, const IntersectionPoint & intersectionPoint);

IntersectionPointFilter keep_all_intersection_points_filter();

std::vector<IntersectionPoint> build_all_intersection_points(const stk::mesh::BulkData & mesh,
    const stk::mesh::Selector & elementSelector,
    const InterfaceGeometry & geometry,
    const NodeToCapturedDomainsMap & nodesToCapturedDomains);

std::vector<IntersectionPoint> build_uncaptured_intersection_points(const stk::mesh::BulkData & mesh,
    const stk::mesh::Selector & elementSelector,
    const InterfaceGeometry & geometry,
    const NodeToCapturedDomainsMap & nodesToCapturedDomains);

std::vector<size_t> update_intersection_points_after_snap_iteration(const stk::mesh::BulkData & mesh,
  const stk::mesh::Selector & elementSelector,
  const InterfaceGeometry & geometry,
  const std::vector<stk::mesh::Entity> & iterationSortedSnapNodes,
  const NodeToCapturedDomainsMap & nodesToCapturedDomains,
  std::vector<IntersectionPoint> & intersectionPoints);

}

#endif /* KRINO_INCLUDE_AKRI_INTERSECTION_POINTS_H_ */
