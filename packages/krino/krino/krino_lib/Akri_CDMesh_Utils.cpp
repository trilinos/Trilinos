// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include <Akri_CDMesh_Utils.hpp>
#include <Akri_CDFEM_Support.hpp>
#include <Akri_InterfaceID.hpp>
#include <Akri_PhaseTag.hpp>
#include <Akri_Phase_Support.hpp>
#include <stk_io/IossBridge.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/MetaData.hpp>

namespace krino {

bool
parts_are_compatible_for_snapping(const stk::mesh::BulkData & mesh, stk::mesh::Entity possible_snap_node, stk::mesh::Entity fixed_node)
{
  const stk::mesh::PartVector & possible_snap_node_parts = mesh.bucket(possible_snap_node).supersets();
  const stk::mesh::PartVector & fixed_node_parts = mesh.bucket(fixed_node).supersets();
  for (auto && possible_snap_node_part : possible_snap_node_parts)
  {
    if ((possible_snap_node_part->primary_entity_rank() == stk::topology::ELEMENT_RANK ||
         possible_snap_node_part->primary_entity_rank() == mesh.mesh_meta_data().side_rank()) &&
        !stk::mesh::is_auto_declared_part(*possible_snap_node_part) &&
        !stk::mesh::contain(fixed_node_parts, *possible_snap_node_part))
    {
      return false;
    }
  }
  return true;
}

static stk::mesh::Part & get_part_to_check(const Phase_Support & phaseSupport, stk::mesh::Part & part)
{
  if (phaseSupport.is_decomposed(part))
    return phaseSupport.find_nonconformal_part(part);
  return part;
}

bool is_part_to_check_for_snapping_compatibility(const Phase_Support & phaseSupport, const AuxMetaData & auxMeta, const stk::mesh::EntityRank targetRank, const stk::mesh::Part & part)
{
  const stk::mesh::Part & exposedBoundaryPart = auxMeta.exposed_boundary_part();
  return part.primary_entity_rank() == targetRank &&
    !phaseSupport.is_interface(part) &&
    (&part == &exposedBoundaryPart || stk::io::is_part_io_part(part) || phaseSupport.is_nonconformal(part));
}

static stk::mesh::PartVector get_nonconformal_parts_to_check(const stk::mesh::BulkData & mesh, const AuxMetaData & auxMeta, const Phase_Support & phaseSupport, const stk::mesh::EntityRank targetRank, const std::vector<stk::mesh::Entity> & targetEntities)
{
  stk::mesh::PartVector partsToCheck;
  for (auto && targetEntity : targetEntities)
    for (auto * part : mesh.bucket(targetEntity).supersets())
      if (is_part_to_check_for_snapping_compatibility(phaseSupport, auxMeta, targetRank, *part))
        partsToCheck.push_back(&get_part_to_check(phaseSupport, *part));
  stk::util::sort_and_unique(partsToCheck, stk::mesh::PartLess());
  return partsToCheck;
}

bool
parts_are_compatible_for_snapping_when_ignoring_phase(const stk::mesh::BulkData & mesh,
    const AuxMetaData & auxMeta,
    const Phase_Support & phaseSupport,
    const stk::mesh::Entity possibleSnapNode,
    const stk::mesh::EntityRank targetRank,
    const stk::mesh::PartVector & nonconformalPartsToCheck)
{
  for (auto * possibleSnapNodePart : mesh.bucket(possibleSnapNode).supersets())
    if (is_part_to_check_for_snapping_compatibility(phaseSupport, auxMeta, targetRank, *possibleSnapNodePart))
      if (!stk::mesh::contain(nonconformalPartsToCheck, get_part_to_check(phaseSupport, *possibleSnapNodePart)))
        return false;
  return true;
}

static stk::topology get_simplex_element_topology(const stk::mesh::BulkData & mesh)
{
  return ((mesh.mesh_meta_data().spatial_dimension() == 2) ? stk::topology::TRIANGLE_3_2D : stk::topology::TETRAHEDRON_4);
}

static void fill_topology_entities(const stk::mesh::BulkData & mesh, const stk::topology & topology, const std::vector<stk::mesh::Entity> & nodes, std::vector<stk::mesh::Entity> & topologyEntities)
{
  topologyEntities.clear();
  if (nodes.size() <= topology.num_nodes())
  {
    stk::mesh::get_entities_through_relations(mesh, nodes, topology.rank(), topologyEntities);
  }
}

std::vector<bool> which_intersection_point_nodes_are_compatible_for_snapping(const stk::mesh::BulkData & mesh, const AuxMetaData & auxMeta, const Phase_Support & phaseSupport, const std::vector<stk::mesh::Entity> & intersectionPointNodes)
{
  std::vector<bool> areIntersectionPointsCompatibleForSnapping(intersectionPointNodes.size(), true);
  filter_which_intersection_point_nodes_are_compatible_for_snapping(mesh, auxMeta, phaseSupport, intersectionPointNodes, areIntersectionPointsCompatibleForSnapping);
  return areIntersectionPointsCompatibleForSnapping;
}

void filter_which_intersection_point_nodes_are_compatible_for_snapping(const stk::mesh::BulkData & mesh,
    const AuxMetaData & auxMeta,
    const Phase_Support & phaseSupport,
    const std::vector<stk::mesh::Entity> & intersectionPointNodes,
    std::vector<bool> & areIntersectionPointsCompatibleForSnapping)
{
  STK_ThrowAssert(intersectionPointNodes.size() == areIntersectionPointsCompatibleForSnapping.size());
  std::vector<stk::mesh::Entity> topologyEntities;
  stk::topology elemTopology = get_simplex_element_topology(mesh);
  std::array<stk::topology,2> sideAndElementTopology{{elemTopology.side_topology(), elemTopology}};
  for (stk::topology topo : sideAndElementTopology)
  {
    fill_topology_entities(mesh, topo, intersectionPointNodes, topologyEntities);
    const stk::mesh::PartVector nonconformalPartsToCheck = get_nonconformal_parts_to_check(mesh, auxMeta, phaseSupport, topo.rank(), topologyEntities);
    for(size_t iNode=0; iNode<intersectionPointNodes.size(); ++iNode)
    {
      areIntersectionPointsCompatibleForSnapping[iNode] =
          areIntersectionPointsCompatibleForSnapping[iNode] &&
          parts_are_compatible_for_snapping_when_ignoring_phase(mesh, auxMeta, phaseSupport, intersectionPointNodes[iNode], topo.rank(), nonconformalPartsToCheck);
    }
  }
}

bool phase_matches_interface(const bool oneLSPerPhase, const std::vector<Surface_Identifier> & surfaceIDs, const PhaseTag & phase, const InterfaceID interface)
{
  if(surfaceIDs.size() > 1 && oneLSPerPhase)
  {
    return (phase.contain(surfaceIDs[interface.first_ls()], -1) &&
            phase.contain(surfaceIDs[interface.second_ls()], -1));
  }
  return (phase.contain(surfaceIDs[interface.first_ls()], -1) &&
         phase.contain(surfaceIDs[interface.first_ls()], +1));
}

void determine_phase_from_parts(PhaseTag & phase, const stk::mesh::PartVector & parts, const Phase_Support & phaseSupport)
{
  STK_ThrowAssert(phase.empty());

  for (auto && part : parts)
    if (part->primary_entity_rank() == stk::topology::ELEMENT_RANK && phaseSupport.is_conformal(*part))
      phase.add(phaseSupport.get_iopart_phase(*part));
}

PhaseTag determine_phase_for_entity(const stk::mesh::BulkData & mesh, stk::mesh::Entity entity, const Phase_Support & phaseSupport)
{
  PhaseTag phase;
  const stk::mesh::PartVector & parts = mesh.bucket(entity).supersets();
  determine_phase_from_parts(phase, parts, phaseSupport);
  return phase;
}

bool node_is_on_interface(const stk::mesh::BulkData & mesh, const Phase_Support & phaseSupport, const std::vector<Surface_Identifier> & surfaceIDs, stk::mesh::Entity node, const InterfaceID & interface)
{
  const PhaseTag nodePhase = determine_phase_for_entity(mesh, node, phaseSupport);
  return phase_matches_interface(phaseSupport.has_one_levelset_per_phase(), surfaceIDs, nodePhase, interface);
}

bool nodes_are_on_any_interface(const stk::mesh::BulkData & /*mesh*/, const Phase_Support & phaseSupport, const stk::mesh::Bucket & nodeBucket)
{
  for(auto * part : nodeBucket.supersets())
    if(phaseSupport.is_interface(*part))
      return true;
  return false;
}

bool node_is_on_any_interface(const stk::mesh::BulkData & mesh, const Phase_Support & phaseSupport, const stk::mesh::Entity node)
{
  return nodes_are_on_any_interface(mesh, phaseSupport, mesh.bucket(node));
}

}


