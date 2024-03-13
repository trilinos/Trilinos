// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include <Akri_DiagWriter.hpp>
#include <Akri_Intersection_Points.hpp>
#include <Akri_MeshHelpers.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/GetEntities.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_util/parallel/CommSparse.hpp>
#include <Akri_InterfaceGeometry.hpp>

namespace krino {

EdgeIntersection::EdgeIntersection(const IntersectionPoint & intersectionPt)
{
  const std::vector<stk::mesh::Entity> & intPtNodes = intersectionPt.get_nodes();
  STK_ThrowAssert(intPtNodes.size() == 2);
  nodes = {intPtNodes[0], intPtNodes[1]};
  crossingLocation = intersectionPt.get_weights()[1];
  const auto & domains = intersectionPt.get_sorted_domains();
  STK_ThrowAssert(domains.size() == 1 || domains.size() == 2);
  interface = (domains.size() == 1) ? InterfaceID(domains[0],domains[0]) : InterfaceID(domains[0], domains[1]);
}

std::string debug_output(const stk::mesh::BulkData & mesh, const IntersectionPoint & intersectionPoint)
{
  const std::vector<stk::mesh::Entity> & nodes = intersectionPoint.get_nodes();
  const std::vector<double> & weights = intersectionPoint.get_weights();
  const std::vector<int> & domains = intersectionPoint.get_sorted_domains();
  std::ostringstream os;
  os << "intersection point domains={ ";
  for (int domain : domains)
    os << domain << " ";
  os << "} stencil={ ";
  for (size_t i=0; i<nodes.size(); ++i)
    os << mesh.identifier(nodes[i]) << "@" << weights[i] << " ";
  os << "} ";
  return os.str();
}

bool first_sorted_vector_of_domains_contains_all_domains_in_second_vector(const std::vector<int> & firstVec, const std::vector<int> & secondVec)
{
  for (int domain : secondVec)
     if (!std::binary_search(firstVec.begin(), firstVec.end(), domain))
       return false;
  return true;
}

bool any_node_already_captures_intersection_point(const std::vector<stk::mesh::Entity> & nodes,
    const std::vector<int> & intersectionPointSortedDomains,
    const NodeToCapturedDomainsMap & nodesToCapturedDomains)
{
  for (auto node : nodes)
  {
    auto iter = nodesToCapturedDomains.find(node);
    if (iter != nodesToCapturedDomains.end())
    {
      const std::vector<int> & nodeSortedDomains = iter->second;
      if (first_sorted_vector_of_domains_contains_all_domains_in_second_vector(nodeSortedDomains, intersectionPointSortedDomains))
        return true;
    }
  }
  return false;
}

bool domains_already_snapped_to_node_are_also_at_intersection_point(const NodeToCapturedDomainsMap & nodesToCapturedDomains, stk::mesh::Entity node, const std::vector<int> & intersectionPointDomains)
{
  const auto iter = nodesToCapturedDomains.find(node);
  if (iter == nodesToCapturedDomains.end())
    return true;
  return first_sorted_vector_of_domains_contains_all_domains_in_second_vector(intersectionPointDomains, iter->second);
}

IntersectionPointFilter
keep_all_intersection_points_filter()
{
  auto filter =
  [](const std::vector<stk::mesh::Entity> & intersectionPointNodes, const std::vector<int> & intersectionPointSortedDomains)
  {
    return true;
  };
  return filter;
}

static
void pack_intersection_points_for_owners_of_nodes(const stk::mesh::BulkData & mesh,
    const size_t indexOfFirstToCommunicate,
    const std::vector<IntersectionPoint> & intersectionPoints,
    stk::CommSparse &commSparse)
{
  std::vector<int> intersectionPointNodeOwners;
  stk::pack_and_communicate(commSparse,[&]()
  {
    for (size_t intPtIndex=indexOfFirstToCommunicate; intPtIndex<intersectionPoints.size(); ++intPtIndex)
    {
      const auto & intersectionPoint = intersectionPoints[intPtIndex];
      intersectionPointNodeOwners.clear();
      for (auto && node : intersectionPoint.get_nodes())
      {
        const int procId = mesh.parallel_owner_rank(node);
        if (procId != commSparse.parallel_rank())
          intersectionPointNodeOwners.push_back(procId);
      }
      stk::util::sort_and_unique(intersectionPointNodeOwners);
      for (int procId : intersectionPointNodeOwners)
      {
        commSparse.send_buffer(procId).pack(intersectionPoint.get_nodes().size());
        for (auto && node : intersectionPoint.get_nodes())
          commSparse.send_buffer(procId).pack(mesh.identifier(node));
        for (auto && weight : intersectionPoint.get_weights())
          commSparse.send_buffer(procId).pack(weight);
        commSparse.send_buffer(procId).pack(intersectionPoint.get_sorted_domains().size());
        for (auto && domain : intersectionPoint.get_sorted_domains())
          commSparse.send_buffer(procId).pack(domain);
      }
    }
  });
}

static
void unpack_intersection_points(const stk::mesh::BulkData & mesh,
    std::vector<IntersectionPoint> & intersectionPoints,
    stk::CommSparse &commSparse)
{
  std::vector<stk::mesh::Entity> intersectionPointNodes;
  std::vector<double> intersectionPointWeights;
  std::vector<int> intersectionPointDomains;
  stk::mesh::EntityId nodeId;
  const bool intersectionPointIsOwned = false;

  stk::unpack_communications(commSparse, [&](int procId)
  {
    stk::CommBuffer & buffer = commSparse.recv_buffer(procId);

    while ( buffer.remaining() )
    {
      size_t numNodes;
      commSparse.recv_buffer(procId).unpack(numNodes);
      intersectionPointNodes.resize(numNodes);
      for (auto && node : intersectionPointNodes)
      {
        commSparse.recv_buffer(procId).unpack(nodeId);
        node = mesh.get_entity(stk::topology::NODE_RANK, nodeId);
      }
      intersectionPointWeights.resize(numNodes);
      for (auto && weight : intersectionPointWeights)
        commSparse.recv_buffer(procId).unpack(weight);
      size_t numDomains;
      commSparse.recv_buffer(procId).unpack(numDomains);
      intersectionPointDomains.resize(numDomains);
      for (auto && domain : intersectionPointDomains)
      {
        commSparse.recv_buffer(procId).unpack(domain);
      }
      intersectionPoints.emplace_back(intersectionPointIsOwned, intersectionPointNodes, intersectionPointWeights, intersectionPointDomains);
    }
  });
}

static void communicate_intersection_points(const stk::mesh::BulkData & mesh, const size_t indexOfFirstToCommunicate, std::vector<IntersectionPoint> & intersectionPoints)
{
  stk::CommSparse commSparse(mesh.parallel());
  pack_intersection_points_for_owners_of_nodes(mesh, indexOfFirstToCommunicate, intersectionPoints, commSparse);
  unpack_intersection_points(mesh, intersectionPoints, commSparse);
}

static void communicate_all_intersection_points(const stk::mesh::BulkData & mesh, std::vector<IntersectionPoint> & intersectionPoints)
{
  communicate_intersection_points(mesh, 0, intersectionPoints);
}

static std::vector<IntersectionPoint> build_intersection_points(const stk::mesh::BulkData & mesh,
    const stk::mesh::Selector & elementSelector,
    const InterfaceGeometry & geometry,
    const NodeToCapturedDomainsMap & nodesToCapturedDomains,
    const IntersectionPointFilter & intersectionPointFilter)
{
  std::vector<stk::mesh::Entity> elementsToIntersect;
  stk::mesh::get_entities( mesh, stk::topology::ELEMENT_RANK, mesh.mesh_meta_data().locally_owned_part() & elementSelector, elementsToIntersect, false);

  std::vector<IntersectionPoint> intersectionPoints;
  geometry.append_element_intersection_points(mesh, nodesToCapturedDomains, elementsToIntersect, intersectionPointFilter, intersectionPoints);
  communicate_all_intersection_points(mesh, intersectionPoints);

  return intersectionPoints;
}

std::vector<IntersectionPoint> build_all_intersection_points(const stk::mesh::BulkData & mesh,
    const stk::mesh::Selector & elementSelector,
    const InterfaceGeometry & geometry,
    const NodeToCapturedDomainsMap & nodesToCapturedDomains)
{
  const IntersectionPointFilter intersectionPointFilter = keep_all_intersection_points_filter();
  return build_intersection_points(mesh, elementSelector, geometry, nodesToCapturedDomains, intersectionPointFilter);
}

static IntersectionPointFilter
filter_intersection_points_to_those_not_already_handled(const NodeToCapturedDomainsMap & nodesToCapturedDomains)
{
  auto filter =
  [&nodesToCapturedDomains](const std::vector<stk::mesh::Entity> & intersectionPointNodes, const std::vector<int> & intersectionPointSortedDomains)
  {
    return !any_node_already_captures_intersection_point(intersectionPointNodes, intersectionPointSortedDomains, nodesToCapturedDomains);
  };
  return filter;
}

std::vector<IntersectionPoint> build_uncaptured_intersection_points(const stk::mesh::BulkData & mesh,
    const stk::mesh::Selector & elementSelector,
    const InterfaceGeometry & geometry,
    const NodeToCapturedDomainsMap & nodesToCapturedDomains)
{
  IntersectionPointFilter intersectionPointFilter = filter_intersection_points_to_those_not_already_handled(nodesToCapturedDomains);

  return build_intersection_points(mesh, elementSelector, geometry, nodesToCapturedDomains, intersectionPointFilter);
}

static IntersectionPointFilter
filter_intersection_points_to_those_using_previous_iteration_snap_nodes_but_not_already_handled(const NodeToCapturedDomainsMap & nodesToCapturedDomains, const std::vector<stk::mesh::Entity> & iterationSortedSnapNodes)
{
  auto filter =
  [&nodesToCapturedDomains, &iterationSortedSnapNodes](const std::vector<stk::mesh::Entity> & intersectionPointNodes, const std::vector<int> & intersectionPointSortedDomains)
  {
    return any_entity_in_first_vector_is_contained_in_second_sorted_vector(intersectionPointNodes, iterationSortedSnapNodes) &&
        !any_node_already_captures_intersection_point(intersectionPointNodes, intersectionPointSortedDomains, nodesToCapturedDomains);
  };
  return filter;
}

std::vector<stk::mesh::Entity> get_owned_elements_using_nodes_knowing_that_nodes_dont_have_common_elements(const stk::mesh::BulkData & mesh,
    const stk::mesh::Selector & elementSelector,
    const std::vector<stk::mesh::Entity> & nodes)
{
  std::vector<stk::mesh::Entity> nodeElements;
  for (auto node : nodes)
  {
    for (auto element : StkMeshEntities{mesh.begin_elements(node), mesh.end_elements(node)})
    {
      const auto & bucket = mesh.bucket(element);
      if (bucket.owned() && elementSelector(bucket))
      {
        nodeElements.push_back(element);
      }
    }
  }
  return nodeElements;
}

std::vector<size_t> update_intersection_points_after_snap_iteration(const stk::mesh::BulkData & mesh,
    const stk::mesh::Selector & elementSelector,
    const InterfaceGeometry & geometry,
    const std::vector<stk::mesh::Entity> & iterationSortedSnapNodes,
    const NodeToCapturedDomainsMap & nodesToCapturedDomains,
    std::vector<IntersectionPoint> & intersectionPoints)
{
  const auto intersectionPointFilter = filter_intersection_points_to_those_using_previous_iteration_snap_nodes_but_not_already_handled(nodesToCapturedDomains, iterationSortedSnapNodes);
  const size_t badIndex = std::numeric_limits<size_t>::max();

  std::vector<size_t> oldToNewIntPts;
  oldToNewIntPts.reserve(intersectionPoints.size());
  size_t newSize=0;
  for (auto && intersectionPoint : intersectionPoints)
  {
    if (!any_entity_in_first_vector_is_contained_in_second_sorted_vector(intersectionPoint.get_nodes(), iterationSortedSnapNodes))
    {
      oldToNewIntPts.push_back(newSize);
      std::swap(intersectionPoint, intersectionPoints[newSize++]);
    }
    else
    {
      oldToNewIntPts.push_back(badIndex);
    }
  }
  intersectionPoints.erase(intersectionPoints.begin()+newSize, intersectionPoints.end());

  if (geometry.snapped_elements_may_have_new_intersections())
  {
    const std::vector<stk::mesh::Entity> updateElements = get_owned_elements_using_nodes_knowing_that_nodes_dont_have_common_elements(mesh, elementSelector, iterationSortedSnapNodes);
    geometry.append_element_intersection_points(mesh, nodesToCapturedDomains, updateElements, intersectionPointFilter, intersectionPoints);
    communicate_intersection_points(mesh, newSize, intersectionPoints);
  }

  return oldToNewIntPts;
}

}
