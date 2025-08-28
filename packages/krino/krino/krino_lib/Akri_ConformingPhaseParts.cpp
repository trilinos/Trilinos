#include <Akri_ConformingPhaseParts.hpp>

#include <Akri_Phase_Support.hpp>
#include <stk_mesh/base/Types.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_util/util/ReportHandler.hpp>
#include <stk_io/IossBridge.hpp>

namespace krino {

void determine_original_undecomposed_part_changes_for_current_parts(
    const Phase_Support & phaseSupport,
    const stk::mesh::PartVector & currentParts,
    const stk::mesh::EntityRank entityRank,
    stk::mesh::PartVector & addParts,
    stk::mesh::PartVector & removeParts)
{
  for(auto * partPtr : currentParts)
  {
    if( partPtr->primary_entity_rank() == entityRank && phaseSupport.is_decomposed(*partPtr) )
    {
      stk::mesh::Part & originalPart = phaseSupport.find_original_part(*partPtr);
      if (&originalPart != partPtr)
      {
        addParts.push_back(&originalPart);
        removeParts.push_back(partPtr);

        for(auto * supersetPartPtr : partPtr->supersets())
          if (!stk::mesh::is_auto_declared_part(*supersetPartPtr))
            removeParts.push_back(supersetPartPtr);
      }
    }
  }
}

void determine_original_undecomposed_part_changes_for_entities(
    const stk::mesh::BulkData & /*mesh*/,
    const stk::mesh::Bucket & bucket,
    const Phase_Support & phaseSupport,
    stk::mesh::Part & childPart,
    stk::mesh::Part & parentPart,
    stk::mesh::Part & activePart,
    stk::mesh::PartVector & addParts,
    stk::mesh::PartVector & removeParts)
{
  addParts.clear();
  removeParts.clear();
  const stk::mesh::PartVector & currentParts = bucket.supersets();
  determine_original_undecomposed_part_changes_for_current_parts(phaseSupport, currentParts, bucket.entity_rank(), addParts, removeParts);

  if (!bucket.member(activePart)) addParts.push_back(&activePart);
  if (bucket.member(childPart)) removeParts.push_back(&childPart);
  if (bucket.member(parentPart)) removeParts.push_back(&parentPart);
}

static void append_part_and_user_supersets(stk::mesh::PartVector & parts, stk::mesh::Part * part)
{
  parts.push_back(part);

  for(auto * superset : part->supersets())
    if (!stk::mesh::is_auto_declared_part(*superset))
      parts.push_back(superset);
}

void append_conforming_part_changes_for_current_parts(const Phase_Support & phaseSupport,
    const stk::mesh::PartVector & currentParts,
    const stk::mesh::EntityRank entityRank,
    const PhaseTag & phase,
    stk::mesh::PartVector & addParts,
    stk::mesh::PartVector & removeParts)
{
  for(auto * part : currentParts)
  {
    if( part->primary_entity_rank() == entityRank && phaseSupport.is_decomposed(*part) )
    {
      stk::mesh::Part & conformal_elem_io_part = phaseSupport.find_conformal_io_part(*part, phase);
      if (&conformal_elem_io_part != part)
      {
        addParts.push_back(&conformal_elem_io_part);
        append_part_and_user_supersets(removeParts, part);
      }
    }
  }
}

void append_nonconforming_part_changes_for_current_parts(
    const Phase_Support & phaseSupport,
    const stk::mesh::PartVector & currentParts,
    const stk::mesh::EntityRank entityRank,
    stk::mesh::Part & childPart,
    stk::mesh::Part & parentPart,
    stk::mesh::Part & activePart,
    stk::mesh::PartVector & addParts,
    stk::mesh::PartVector & removeParts)
{
  for(auto * part : currentParts)
  {
    if( part->primary_entity_rank() == entityRank && phaseSupport.is_decomposed(*part))
    {
      stk::mesh::Part & nonconformal_io_part = phaseSupport.find_nonconformal_part(*part);
      if (&nonconformal_io_part != part)
      {
        addParts.push_back(&nonconformal_io_part);
        append_part_and_user_supersets(removeParts, part);
      }
    }
  }

  removeParts.push_back(&activePart);

  if (entityRank == stk::topology::ELEMENT_RANK)
  {
    addParts.push_back(&parentPart);
    removeParts.push_back(&childPart);
  }
}

void append_active_part_changes(
    const bool isEntityActive,
    stk::mesh::Part & activePart,
    stk::mesh::PartVector & addParts,
    stk::mesh::PartVector & removeParts)
{
  if(isEntityActive)
    addParts.push_back(&activePart);
  else
    removeParts.push_back(&activePart);
}

void probe_volume_parts_of_side(const Phase_Support & phaseSupport,
    const stk::mesh::PartVector & currentSideParts,
    int &  numVolumeParts,
    int &  numNonconformingVolumeParts,
    stk::mesh::PartVector & conformingVolumeParts)
{
  numVolumeParts = 0;
  numNonconformingVolumeParts = 0;
  conformingVolumeParts.clear();

  for(auto * sidePart : currentSideParts)
  {
    if (sidePart->primary_entity_rank() == stk::topology::ELEMENT_RANK)
    {
      if (phaseSupport.is_conformal(*sidePart))
        conformingVolumeParts.push_back(sidePart);

      if (phaseSupport.is_nonconformal(*sidePart))
        ++numNonconformingVolumeParts;
      else if (stk::io::is_part_io_part(*sidePart) && !stk::io::is_part_assembly_io_part(*sidePart))
        ++numVolumeParts;
    }
  }
}

static void add_correct_interface_parts_and_remove_other_interface_parts(const Phase_Support & phaseSupport, const stk::mesh::PartVector & currentSideParts, const stk::mesh::PartVector & conformingVolumeParts, stk::mesh::PartVector & addParts, stk::mesh::PartVector & removeParts)
{
  stk::mesh::Part * conformingSide0 = const_cast<stk::mesh::Part *>(phaseSupport.find_interface_part(*conformingVolumeParts[0], *conformingVolumeParts[1]));
  if (nullptr != conformingSide0) addParts.push_back(conformingSide0);
  stk::mesh::Part * conformingSide1 = const_cast<stk::mesh::Part *>(phaseSupport.find_interface_part(*conformingVolumeParts[1], *conformingVolumeParts[0]));
  if (nullptr != conformingSide1) addParts.push_back(conformingSide1);

  for(auto * part : currentSideParts)
    if (phaseSupport.is_interface(*part) && part != conformingSide0 && part != conformingSide1)
      append_part_and_user_supersets(removeParts, part);
}

static void remove_all_interface_parts(const Phase_Support & phaseSupport, const stk::mesh::PartVector & existingParts, stk::mesh::PartVector & removeParts)
{
  for(auto * part : existingParts)
    if (phaseSupport.is_interface(*part))
      append_part_and_user_supersets(removeParts, part);
}

void determine_part_changes_for_side(
    const Phase_Support & phaseSupport,
    const stk::mesh::EntityRank sideRank,
    const stk::mesh::PartVector & currentSideParts,
    const bool doesSideHaveActiveElements,
    stk::mesh::Part & blockBoundarySide,
    stk::mesh::Part & childPart,
    stk::mesh::Part & parentPart,
    stk::mesh::Part & activePart,
    stk::mesh::PartVector & addParts,
    stk::mesh::PartVector & removeParts)
{
  addParts.clear();
  removeParts.clear();

  int numVolumeParts = 0;
  int numNonconformingVolumeParts = 0;
  stk::mesh::PartVector conformingVolumeParts;
  probe_volume_parts_of_side(phaseSupport, currentSideParts, numVolumeParts, numNonconformingVolumeParts, conformingVolumeParts);

  if (numVolumeParts == 2)
    addParts.push_back(&blockBoundarySide);
  else
    removeParts.push_back(&blockBoundarySide);

  if (conformingVolumeParts.empty())
  {
    /* There are two possible cases where no conformal volume parts are found:
     *   1) This side is part of a surface that does not touch any blocks that are being decomposed.
     *      Only the active parts for these sides should be updated.
     *   2) This side is a parent side that should be deactivated and moved to the nonconformal part.
     *      These sides will have at least 1 nonconformal volume part from the parent volume element.
     */
    if(0 == numNonconformingVolumeParts)
      append_active_part_changes(doesSideHaveActiveElements, activePart, addParts, removeParts);
    else
      append_nonconforming_part_changes_for_current_parts(phaseSupport, currentSideParts, sideRank, childPart, parentPart, activePart, addParts, removeParts);

    return;
  }

  STK_ThrowRequire(conformingVolumeParts.size() == 1 || conformingVolumeParts.size() == 2);

  std::vector<PhaseTag> side_phases(conformingVolumeParts.size());
  for (unsigned iphase = 0; iphase<side_phases.size(); ++iphase)
  {
    side_phases[iphase] = phaseSupport.get_iopart_phase(*conformingVolumeParts[iphase]);
    STK_ThrowRequire(!side_phases[iphase].empty());
  }

  if (conformingVolumeParts.size() == 2 && side_phases[0] != side_phases[1])
    add_correct_interface_parts_and_remove_other_interface_parts(phaseSupport, currentSideParts, conformingVolumeParts, addParts, removeParts);
  else
    remove_all_interface_parts(phaseSupport, currentSideParts, removeParts);

  for (auto && side_phase : side_phases)
    append_conforming_part_changes_for_current_parts(phaseSupport, currentSideParts, sideRank, side_phase, addParts, removeParts);

  append_active_part_changes(doesSideHaveActiveElements, activePart, addParts, removeParts);
}

void determine_child_conforming_parts(const stk::mesh::MetaData & meta,
    const Phase_Support & phaseSupport,
    const stk::topology topology,
    const stk::mesh::PartVector & parentParts,
    const stk::mesh::PartVector & attributeParts,
    stk::mesh::Part & childPart,
    stk::mesh::Part & activePart,
    const PhaseTag & phase,
    stk::mesh::PartVector & childParts)
{
  childParts.clear();

  stk::mesh::EntityRank entityRank = topology.rank();
  for(auto * part : parentParts)
  {
    if( part->primary_entity_rank() == entityRank &&
        (stk::io::is_part_io_part(part) || phaseSupport.is_decomposed(*part)) &&
        !phaseSupport.is_interface(*part) )
      childParts.push_back(&phaseSupport.find_conformal_io_part(*part, phase));
    else if (stk::mesh::contain(attributeParts, *part))
      childParts.push_back(part);
  }

  childParts.push_back(&meta.get_topology_root_part(topology));

  if (entityRank == stk::topology::ELEMENT_RANK)
    childParts.push_back(&childPart);

  childParts.push_back(&activePart);
}

}
