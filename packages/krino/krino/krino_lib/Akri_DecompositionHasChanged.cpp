// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include <Akri_CDFEM_Parent_Edges.hpp>
#include <Akri_CDFEM_Support.hpp>
#include <Akri_CDMesh_Utils.hpp>
#include <Akri_DecompositionHasChanged.hpp>
#include <Akri_DiagWriter.hpp>
#include <Akri_InterfaceGeometry.hpp>
#include <Akri_Intersection_Points.hpp>
#include <Akri_MeshHelpers.hpp>
#include <Akri_NodeToCapturedDomains.hpp>
#include <Akri_ParentsToChildMapper.hpp>
#include <Akri_Phase_Support.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_util/parallel/CommSparse.hpp>
#include <stk_util/parallel/ParallelReduceBool.hpp>
#include <array>

namespace krino {

static bool owned_nodes_are_handled_by_other_procs(const stk::mesh::BulkData & mesh,
    const std::vector<stk::mesh::Entity> & interfaceNodesWithoutMatchingCrossing,
    const std::vector<stk::mesh::Entity> & interfaceNodesWithMatchingCrossing)
{
  stk::CommSparse commSparse(mesh.parallel());
  pack_entities_for_owning_proc(mesh, interfaceNodesWithMatchingCrossing, commSparse);
  std::set<stk::mesh::Entity> ownedNodesWithRemoteMatchingCrossing;
  unpack_entities_from_other_procs(mesh, ownedNodesWithRemoteMatchingCrossing, commSparse);

  for (auto node : interfaceNodesWithoutMatchingCrossing)
    if (mesh.bucket(node).owned() && ownedNodesWithRemoteMatchingCrossing.count(node) == 0)
      return false;

  return true;
}

static bool any_node_of_edge_including_children_is_on_interface(const stk::mesh::BulkData & mesh,
    const Phase_Support & phaseSupport,
    const std::vector<Surface_Identifier> & surfaceIDs,
    const ParentsToChildMapper & parentToChildMapper,
    const EdgeIntersection & edgeCrossing)
{
  std::vector<stk::mesh::Entity> edgeNodesIncludingChildren;
  fill_edge_nodes(mesh, edgeCrossing.nodes[0], edgeCrossing.nodes[1], parentToChildMapper, edgeNodesIncludingChildren);
  for (auto node : edgeNodesIncludingChildren)
    if (node_is_on_interface(mesh, phaseSupport, surfaceIDs, node, edgeCrossing.interface))
      return true;
  return false;
}

static bool edge_crossings_have_matching_interface_node(const stk::mesh::BulkData & mesh,
    const Phase_Support & phaseSupport,
    const std::vector<Surface_Identifier> & surfaceIDs,
    const ParentsToChildMapper & parentToChildMapper,
    const EdgeIntersection & edgeCrossing,
    const double snapTol)
{
  const auto & edgeNodes = edgeCrossing.nodes;
  const auto & crossingInterface = edgeCrossing.interface;
  const double crossingPosition = edgeCrossing.crossingLocation;
  if (crossingPosition < snapTol)
  {
    if (!node_is_on_interface(mesh, phaseSupport, surfaceIDs, edgeNodes[0], crossingInterface))
      return false;
  }
  else if (crossingPosition > 1.-snapTol)
  {
    if (!node_is_on_interface(mesh, phaseSupport, surfaceIDs, edgeNodes[1], crossingInterface))
      return false;
  }
  else if (!any_node_of_edge_including_children_is_on_interface(mesh, phaseSupport, surfaceIDs, parentToChildMapper, edgeCrossing))
  {
    return false;
  }
  return true;
}

static bool edges_with_crossings_have_matching_interface_nodes(const stk::mesh::BulkData & mesh,
    const CDFEM_Support & cdfemSupport,
    const Phase_Support & phaseSupport,
    const std::vector<Surface_Identifier> & surfaceIDs,
    const ParentsToChildMapper & parentToChildMapper,
    const std::vector<IntersectionPoint> & edgeIntersections)
{
  const double snapTol = cdfemSupport.get_snapper().get_edge_tolerance();
  bool edgesWithCrossingHaveMatchingInterfaceNodes = true;
  for (auto && edgeIntersection : edgeIntersections)
  {
    const EdgeIntersection edge(edgeIntersection);
    if (!edge_crossings_have_matching_interface_node(mesh, phaseSupport, surfaceIDs, parentToChildMapper, edge, snapTol))
    {
      edgesWithCrossingHaveMatchingInterfaceNodes = false;
      break;
    }
  }
  return stk::is_true_on_all_procs(mesh.parallel(), edgesWithCrossingHaveMatchingInterfaceNodes);
}

static std::map<stk::mesh::Entity, std::set<InterfaceID>> build_nodes_to_interfaces_within_tolerance(const stk::mesh::BulkData & mesh,
    const CDFEM_Support & cdfemSupport,
    const std::vector<IntersectionPoint> & edgeIntersections,
    const ParentsToChildMapper & parentToChildMapper)
{

  const double snapTol = cdfemSupport.get_snapper().get_edge_tolerance();
  std::vector<stk::mesh::Entity> parentEdgeNodes;
  std::vector<double> parentEdgeNodePositions;

  std::map<stk::mesh::Entity, std::set<InterfaceID>> nodesToInterfacesWithinTolerance;
  for (auto && edgeIntersection : edgeIntersections)
  {
    const EdgeIntersection edgeCrossing(edgeIntersection);
    const auto & crossingInterface = edgeCrossing.interface;
    const double crossingPosition = edgeCrossing.crossingLocation;

    parentEdgeNodes.clear();
    parentEdgeNodePositions.clear();
    fill_edge_nodes_and_positions(mesh, edgeCrossing.nodes[0], edgeCrossing.nodes[1], parentToChildMapper, parentEdgeNodes, parentEdgeNodePositions);

    for (size_t iNode=0; iNode<parentEdgeNodes.size(); ++iNode)
      if (std::abs(parentEdgeNodePositions[iNode]-crossingPosition) < snapTol)
        nodesToInterfacesWithinTolerance[parentEdgeNodes[iNode]].insert(crossingInterface);
  }
  return nodesToInterfacesWithinTolerance;
}

static bool node_has_matching_interface_within_tolerance(const stk::mesh::BulkData & mesh,
    const Phase_Support & phaseSupport,
    const std::vector<Surface_Identifier> & surfaceIDs,
    const stk::mesh::Entity node,
    const std::map<stk::mesh::Entity, std::set<InterfaceID>> & nodesToInterfacesWithinTolerance)
{
  const auto mapIter = nodesToInterfacesWithinTolerance.find(node);
  if (mapIter == nodesToInterfacesWithinTolerance.end())
    return false;

  const PhaseTag nodePhase = determine_phase_for_entity(mesh, node, phaseSupport);
  const auto & nodeInterfacesWithinTolerance = mapIter->second;
  for (auto && interfaceWithinTolerance : nodeInterfacesWithinTolerance)
    if (phase_matches_interface(phaseSupport.has_one_levelset_per_phase(), surfaceIDs, nodePhase, interfaceWithinTolerance))
      return true;

  return false;
}

static void fill_interface_nodes_with_and_without_matching_crossing(const stk::mesh::BulkData & mesh,
    const CDFEM_Support & cdfemSupport,
    const Phase_Support & phaseSupport,
    const std::vector<Surface_Identifier> & surfaceIDs,
    const ParentsToChildMapper & parentToChildMapper,
    const std::vector<IntersectionPoint> & edgeIntersections,
    std::vector<stk::mesh::Entity> & interfaceNodesWithMatchingCrossing,
    std::vector<stk::mesh::Entity> & interfaceNodesWithoutMatchingCrossing)
{
  interfaceNodesWithMatchingCrossing.clear();
  interfaceNodesWithoutMatchingCrossing.clear();

  const std::map<stk::mesh::Entity, std::set<InterfaceID>> nodesToInterfacesWithinTolerance = build_nodes_to_interfaces_within_tolerance(mesh, cdfemSupport, edgeIntersections, parentToChildMapper);

  for ( auto && bucket : mesh.buckets(stk::topology::NODE_RANK) )
  {
    if (nodes_are_on_any_interface(mesh, phaseSupport, *bucket))
    {
      for ( auto && node : *bucket )
      {
        if (node_has_matching_interface_within_tolerance(mesh, phaseSupport, surfaceIDs, node, nodesToInterfacesWithinTolerance))
          interfaceNodesWithMatchingCrossing.push_back(node);
        else
          interfaceNodesWithoutMatchingCrossing.push_back(node);
      }
    }
  }
}

static bool interface_nodes_have_matching_edge_crossing(const stk::mesh::BulkData & mesh,
    const CDFEM_Support & cdfemSupport,
    const Phase_Support & phaseSupport,
    const std::vector<Surface_Identifier> & surfaceIDs,
    const ParentsToChildMapper & parentToChildMapper,
    const std::vector<IntersectionPoint> & edgeIntersections)
{
  std::vector<stk::mesh::Entity> interfaceNodesWithMatchingCrossing;
  std::vector<stk::mesh::Entity> interfaceNodesWithoutMatchingCrossing;
  fill_interface_nodes_with_and_without_matching_crossing(mesh, cdfemSupport, phaseSupport, surfaceIDs, parentToChildMapper, edgeIntersections, interfaceNodesWithMatchingCrossing, interfaceNodesWithoutMatchingCrossing);

  if (stk::is_true_on_all_procs(mesh.parallel(), interfaceNodesWithoutMatchingCrossing.empty()))
    return true;

  bool allOwnedNotSharedInterfaceNodesHaveMatchingEdgeCrossing = true;
  for (auto node : interfaceNodesWithoutMatchingCrossing)
  {
    if (mesh.bucket(node).owned() && !mesh.bucket(node).shared())
    {
      allOwnedNotSharedInterfaceNodesHaveMatchingEdgeCrossing = false;
      break;
    }
  }

  if (!stk::is_true_on_all_procs(mesh.parallel(), allOwnedNotSharedInterfaceNodesHaveMatchingEdgeCrossing))
    return false;

  const bool allInterfaceNodesHaveMatchingEdgeCrossing =
      owned_nodes_are_handled_by_other_procs(mesh, interfaceNodesWithoutMatchingCrossing, interfaceNodesWithMatchingCrossing);

  return stk::is_true_on_all_procs(mesh.parallel(), allInterfaceNodesHaveMatchingEdgeCrossing);
}

void fill_neighbor_nodes(const stk::mesh::BulkData & mesh,
    stk::mesh::Entity node,
    const stk::mesh::Selector & elementSelector,
    std::vector<stk::mesh::Entity> & neighborNodes)
{
  neighborNodes.clear();
  for (auto element : StkMeshEntities{mesh.begin_elements(node), mesh.end_elements(node)})
    if (elementSelector(mesh.bucket(element)))
      for (auto nbr : StkMeshEntities{mesh.begin_nodes(element), mesh.end_nodes(element)})
        if (nbr != node)
          neighborNodes.push_back(nbr);
  stk::util::sort_and_unique(neighborNodes);
}

static bool snap_displacements_at_node_are_small_compared_to_parent_edges(const stk::mesh::BulkData & mesh,
    stk::mesh::Entity node,
    FieldRef cdfemSnapField,
    FieldRef coordsField,
    const stk::mesh::Selector & parentElementSelector,
    const double snapTol)
{
  const stk::math::Vector3d snapDisplacements(field_data<double>(cdfemSnapField, node), mesh.mesh_meta_data().spatial_dimension());
  const double snapDisplacmentsSqrLength = snapDisplacements.length_squared();
  const stk::math::Vector3d nodeCoords(field_data<double>(coordsField, node), mesh.mesh_meta_data().spatial_dimension());
  std::vector<stk::mesh::Entity> neighborNodes;
  fill_neighbor_nodes(mesh, node, parentElementSelector, neighborNodes);
  for (auto && nbr : neighborNodes)
  {
    const stk::math::Vector3d nbrCoords(field_data<double>(coordsField, nbr), mesh.mesh_meta_data().spatial_dimension());
    const double edgeSqrLength = (nbrCoords-nodeCoords).length_squared();
    if (snapDisplacmentsSqrLength > snapTol*snapTol * edgeSqrLength)
      return false;
  }
  return true;
}

static bool locally_snap_displacements_are_small_on_nodes_with_interfaces(const stk::mesh::BulkData & mesh,
    const stk::mesh::Part & activePart,
    const CDFEM_Support & cdfemSupport,
    const Phase_Support & phaseSupport)
{
  FieldRef cdfemSnapField = cdfemSupport.get_cdfem_snap_displacements_field();
  if (cdfemSnapField.valid())
  {
    const double snapTol = cdfemSupport.get_snapper().get_edge_tolerance();
    FieldRef oldCdfemSnapField = cdfemSupport.get_cdfem_snap_displacements_field().field_state(stk::mesh::StateOld);
    const stk::mesh::Selector ownedParentElementSelector = get_potential_cdfem_parent_element_selector(activePart, cdfemSupport) & mesh.mesh_meta_data().locally_owned_part();
    const FieldRef coordsField(mesh.mesh_meta_data().coordinate_field());
    for ( auto && bucket : mesh.buckets(stk::topology::NODE_RANK) )
      if (bucket->field_data_is_allocated(oldCdfemSnapField) && nodes_are_on_any_interface(mesh, phaseSupport, *bucket))
        for ( auto && node : *bucket )
          if (!snap_displacements_at_node_are_small_compared_to_parent_edges(mesh, node, oldCdfemSnapField, coordsField, ownedParentElementSelector, snapTol))
            return false;
  }
  return true;
}

static bool snap_displacements_are_small_on_nodes_with_interfaces(const stk::mesh::BulkData & mesh,
    const stk::mesh::Part & activePart,
    const CDFEM_Support & cdfemSupport,
    const Phase_Support & phaseSupport)
{
  return stk::is_true_on_all_procs(mesh.parallel(), locally_snap_displacements_are_small_on_nodes_with_interfaces(mesh, activePart, cdfemSupport, phaseSupport));
}

bool decomposition_has_changed(const stk::mesh::BulkData & mesh,
    const InterfaceGeometry & interfaceGeometry,
    const stk::mesh::Part & activePart,
    const CDFEM_Support & cdfemSupport,
    const Phase_Support & phaseSupport)
{
  if (!snap_displacements_are_small_on_nodes_with_interfaces(mesh, activePart, cdfemSupport, phaseSupport))
    return true;

  const NodeToCapturedDomainsMap nodesToCapturedDomains;
  const std::vector<IntersectionPoint> edgeIntersections = interfaceGeometry.get_edge_intersection_points(mesh, nodesToCapturedDomains);

  const bool addHigherOrderMidSideNodes = false;
  ParentsToChildMapper parentToChildMapper;
  parentToChildMapper.build_map(mesh, activePart, cdfemSupport, addHigherOrderMidSideNodes);

  const std::vector<Surface_Identifier> & surfaceIDs = interfaceGeometry.get_surface_identifiers();

  if (edges_with_crossings_have_matching_interface_nodes(mesh, cdfemSupport, phaseSupport, surfaceIDs, parentToChildMapper, edgeIntersections) &&
      interface_nodes_have_matching_edge_crossing(mesh, cdfemSupport, phaseSupport, surfaceIDs, parentToChildMapper, edgeIntersections))
  {
    return false;
  }

  return true;
}


} // namespace
