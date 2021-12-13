// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include <Akri_CDFEM_Support.hpp>
#include <Akri_InterfaceID.hpp>
#include <Akri_PhaseTag.hpp>
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
        possible_snap_node_part->name().compare(0,7,"refine_") != 0 &&
        !stk::mesh::contain(fixed_node_parts, *possible_snap_node_part))
    {
      return false;
    }
  }
  return true;
}

static stk::mesh::Part * get_nonconformal_part(const Phase_Support & phaseSupport, stk::mesh::Part * part)
{
  return const_cast<stk::mesh::Part *>(phaseSupport.find_nonconformal_part(*part));
}

static bool is_part_to_check(const Phase_Support & phaseSupport, const AuxMetaData & auxMeta, const stk::mesh::Part & part)
{
  const stk::mesh::Part & exposedBoundaryPart = auxMeta.exposed_boundary_part();
  return part.primary_entity_rank() != stk::topology::INVALID_RANK &&
    (&part == &exposedBoundaryPart || stk::io::is_part_io_part(part)) &&
    part.name().compare(0,7,"refine_") != 0 &&
    !phaseSupport.is_interface(&part);
}

static stk::mesh::PartVector get_nonconformal_parts_to_check(const AuxMetaData & auxMeta, const Phase_Support & phaseSupport, const stk::mesh::PartVector & inputParts)
{
  stk::mesh::PartVector partsToCheck;
  partsToCheck.reserve(inputParts.size());
  for (auto && part : inputParts)
    if (is_part_to_check(phaseSupport, auxMeta, *part))
      partsToCheck.push_back(get_nonconformal_part(phaseSupport, part));
  stk::util::sort_and_unique(partsToCheck, stk::mesh::PartLess());
  return partsToCheck;
}

bool
parts_are_compatible_for_snapping_when_ignoring_phase(const stk::mesh::BulkData & mesh, const AuxMetaData & auxMeta, const Phase_Support & phaseSupport, stk::mesh::Entity possibleSnapNode, stk::mesh::Entity fixedNode)
{
  const stk::mesh::PartVector & possibleSnapNodeParts = mesh.bucket(possibleSnapNode).supersets();
  const stk::mesh::PartVector nonconformalPartsToCheck = get_nonconformal_parts_to_check(auxMeta, phaseSupport, mesh.bucket(fixedNode).supersets());
  for (auto && possibleSnapNodePart : possibleSnapNodeParts)
  {
    if (is_part_to_check(phaseSupport, auxMeta, *possibleSnapNodePart))
    {
      stk::mesh::Part * nonconformalPart = get_nonconformal_part(phaseSupport, possibleSnapNodePart);
      if (!stk::mesh::contain(nonconformalPartsToCheck, *nonconformalPart))
        return false;
    }
  }
  return true;
}

bool phase_matches_interface(const CDFEM_Support & cdfemSupport, const PhaseTag & phase, const InterfaceID interface)
{
  if(cdfemSupport.num_ls_fields() > 1 && Phase_Support::has_one_levelset_per_phase())
  {
    return (phase.contain(cdfemSupport.ls_field(interface.first_ls()).identifier, -1) &&
            phase.contain(cdfemSupport.ls_field(interface.second_ls()).identifier, -1));
  }
  return (phase.contain(cdfemSupport.ls_field(interface.first_ls()).identifier, -1) &&
         phase.contain(cdfemSupport.ls_field(interface.first_ls()).identifier, +1));
}

bool determine_phase_from_parts(PhaseTag & phase, const stk::mesh::PartVector & parts, const Phase_Support & phaseSupport)
{
  ThrowAssert(phase.empty());
  bool has_conformal_ioparts = false;

  for (auto && part : parts)
  {
    if (part->primary_entity_rank() != stk::topology::ELEMENT_RANK || // limit ourselves to phase-specific volumes
        !(stk::io::is_part_io_part(*part) ||
        phaseSupport.is_nonconformal(part)))
      continue;

    const PhaseTag & iopart_phase = phaseSupport.get_iopart_phase(*part);

    if (!iopart_phase.empty())
      has_conformal_ioparts = true;

    phase.add(iopart_phase);
  }

  return (has_conformal_ioparts);
}

PhaseTag determine_phase_for_entity(const stk::mesh::BulkData & mesh, stk::mesh::Entity entity, const Phase_Support & phaseSupport)
{
  PhaseTag phase;
  const stk::mesh::PartVector & parts = mesh.bucket(entity).supersets();
  determine_phase_from_parts(phase, parts, phaseSupport);
  return phase;
}

bool node_is_on_interface(const stk::mesh::BulkData & mesh, const CDFEM_Support & cdfemSupport, const Phase_Support & phaseSupport, stk::mesh::Entity node, const InterfaceID & interface)
{
  const PhaseTag nodePhase = determine_phase_for_entity(mesh, node, phaseSupport);
  return phase_matches_interface(cdfemSupport, nodePhase, interface);
}

bool nodes_are_on_any_interface(const stk::mesh::BulkData & mesh, const Phase_Support & phaseSupport, const stk::mesh::Bucket & nodeBucket)
{
  auto side_rank = mesh.mesh_meta_data().side_rank();
  for(auto && part : nodeBucket.supersets())
    if(part->primary_entity_rank() == side_rank && stk::io::is_part_io_part(*part) && phaseSupport.is_interface(part))
      return true;
  return false;
}

bool node_is_on_any_interface(const stk::mesh::BulkData & mesh, const Phase_Support & phaseSupport, const stk::mesh::Entity node)
{
  return nodes_are_on_any_interface(mesh, phaseSupport, mesh.bucket(node));
}

}


