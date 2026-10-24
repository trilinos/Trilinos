// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include <Akri_Phase_Support.hpp>

#include <algorithm>
#include <cctype>
#include <map>
#include <set>
#include <sstream>

#include <stk_mesh/base/MetaData.hpp>

#include <Akri_AuxMetaData.hpp>
#include <Akri_DiagWriter.hpp>
#include <Akri_TopologyData.hpp>
#include <stk_util/util/ReportHandler.hpp>

#include <stk_util/environment/RuntimeDoomed.hpp>
#include <Akri_PhaseUtils.hpp>
#include <Akri_Utility.hpp>

namespace krino{

namespace {

std::string build_part_name(const Phase_Support & ps,
  const stk::mesh::Part & part, const PhaseTag & phaseTag,
  const std::string & suffix)
{
  const stk::mesh::Part * parentPart = get_parent_sideset(part);
  if(parentPart == nullptr)
    return part.name() + suffix;

  const std::set<stk::mesh::PartOrdinal> touchingBlockOrdinals = ps.get_input_block_surface_connectivity().get_blocks_touching_surface(part.mesh_meta_data_ordinal());
  STK_ThrowRequireMsg(touchingBlockOrdinals.size() > 0,
      "krino::Akri_Phase_Support: Side block must be touching at least 1 block");

  stk::mesh::Part & touchingBlock = part.mesh_meta_data().get_part(*touchingBlockOrdinals.begin());

  std::string blockOrTopologyName;
  if (is_sideblock_split_by_topology(part.name(), touchingBlock.topology()))
  {
    blockOrTopologyName = get_ioss_topology_name(touchingBlock.topology());
  }
  else
  {
    STK_ThrowRequireMsg(touchingBlockOrdinals.size() == 1,
      "krino::Akri_Phase_Support: Side blocks split by element block should touch exactly 1 block");
    auto * conformingBlockPtr = ps.find_conforming_part_if_decomposed_by_phase(touchingBlock, phaseTag);
    if(!conformingBlockPtr)
    {
      conformingBlockPtr = part.mesh_meta_data().get_part(touchingBlock.name() + suffix);
      if(!conformingBlockPtr)
        conformingBlockPtr = &touchingBlock;
    }
    blockOrTopologyName = conformingBlockPtr->name();
  }

  return generate_ioss_sideblock_name(parentPart->name() + suffix, blockOrTopologyName, get_ioss_topology_name(part.topology()));
}

stk::mesh::Selector select_union_from_part_ordinals(const stk::mesh::MetaData & meta, const std::vector<unsigned> & partOrdinals)
{
  const stk::mesh::PartVector & allParts = meta.get_parts();
  stk::mesh::PartVector parts;
  for (auto && partOrdinal : partOrdinals)
    parts.push_back(allParts[partOrdinal]);
  return stk::mesh::selectUnion(parts);
}

void insert_into_touching_blocks(const stk::mesh::Part & touchingBlock, std::vector<const stk::mesh::Part*> & touchingBlocks)
{
  if (std::find(touchingBlocks.begin(), touchingBlocks.end(), &touchingBlock) == touchingBlocks.end())
    touchingBlocks.push_back(&touchingBlock);
}

bool do_blocks_have_topology_that_can_be_decomposed(const stk::mesh::PartVector & blocksToDecompose)
{
  for(auto * block : blocksToDecompose)
  {
    stk::topology topo = block->topology();
    if(topo != stk::topology::TET_4 && topo != stk::topology::TET_10
      && topo != stk::topology::SHELL_TRI_3 && topo != stk::topology::SHELL_TRI_6
      && topo != stk::topology::SHELL_TRI_3_ALL_FACE_SIDES && topo != stk::topology::SHELL_TRI_6_ALL_FACE_SIDES
      && topo != stk::topology::TRI_3_2D && topo != stk::topology::TRI_6_2D)
    {
      stk::RuntimeDoomedAdHoc() << "Cannot use CDFEM with block " << block->name() << " because only Tetrahedron 4 or 10, Shell_Triangle 3 or 6 (in 3D) and Triangle 3 or 6 (in 2D) element types are supported.\n";
      return false;
    }
  }
  return true;
}

bool does_block_have_any_of_surface_phases(const std::set<unsigned> & blockPhases, const std::vector<unsigned> surfacePhases)
{
  for(auto surfPhase : surfacePhases)
    if (blockPhases.find(surfPhase) != blockPhases.end())
      return true;
  return false;
}

struct BlockDecomposition
{
  stk::mesh::PartOrdinal block;
  unsigned decompIndex;
};

using PhaseTagToBlockDecompositions = std::map<PhaseTag, std::vector<BlockDecomposition>>;

std::string describe_phase_name(const std::string & phaseName)
{
  if (phaseName.empty())
    return "\"\" (the unnamed alive phase of CDFEM death)";
  return "\"" + phaseName + "\"";
}

std::vector<std::pair<std::string, PhaseTagToBlockDecompositions>>
get_conflicting_phase_names(const DecompositionPackage & decomps, const std::set<stk::mesh::PartOrdinal> & touchingBlocks)
{
  std::map<std::string, PhaseTagToBlockDecompositions> phaseNameToTagToBlocks;
  for (unsigned i=0; i<decomps.size(); ++i)
    for (auto * block : decomps.get_blocks_to_decompose(i))
      if (touchingBlocks.count(block->mesh_meta_data_ordinal()) > 0)
        for (auto && phase : decomps.get_phases(i))
          phaseNameToTagToBlocks[phase.name()][phase.tag()].push_back({block->mesh_meta_data_ordinal(), i});

  std::vector<std::pair<std::string, PhaseTagToBlockDecompositions>> conflictingPhaseNames;
  for (auto && [phaseName, tagToBlocks] : phaseNameToTagToBlocks)
    if (tagToBlocks.size() > 1)
      conflictingPhaseNames.emplace_back(phaseName, tagToBlocks);
  return conflictingPhaseNames;
}

}

bool Phase_Support::is_cdfem_use_case() const
{
  return !(get_all_decomposed_blocks_selector() == stk::mesh::Selector());
}

void Phase_Support::force_cdfem_use_case_for_minimal_unit_tests()
{
  all_decomposed_blocks_selector = meta().universal_part();
}

std::vector<unsigned> Phase_Support::get_negative_levelset_interface_ordinals(const Surface_Identifier levelSetIdentifier) const
{
  std::vector<unsigned> negLevelsetInterfaceOrdinals;
  for (auto * part : meta().get_mesh_parts())
  {
    const unsigned partOrd = part->mesh_meta_data_ordinal();
    if (myPhasePartInfo.is_interface_part(partOrd) &&
        myPhasePartInfo.get_interface_part_touching_phase(partOrd).contain(levelSetIdentifier, -1) &&
        myPhasePartInfo.get_interface_part_opposite_phase(partOrd).contain(levelSetIdentifier, +1))
      negLevelsetInterfaceOrdinals.push_back(partOrd);
  }
  return negLevelsetInterfaceOrdinals;
}

std::vector<unsigned> Phase_Support::get_negative_levelset_block_ordinals(const Surface_Identifier levelSetIdentifier) const
{
  std::vector<unsigned> negLevelsetBlockOrdinals;
  for (auto * part : get_conformal_parts_of_rank(stk::topology::ELEMENT_RANK))
  {
    const unsigned partOrd = part->mesh_meta_data_ordinal();
    if (myPhasePartInfo.get_conforming_part_phase(partOrd).contain(levelSetIdentifier, -1))
      negLevelsetBlockOrdinals.push_back(partOrd);
  }
  return negLevelsetBlockOrdinals;
}

stk::mesh::Selector
Phase_Support::get_negative_levelset_interface_selector(const Surface_Identifier levelSetIdentifier) const
{
  return select_union_from_part_ordinals(meta(), get_negative_levelset_interface_ordinals(levelSetIdentifier));
}

stk::mesh::Selector
Phase_Support::get_negative_levelset_block_selector(const Surface_Identifier levelSetIdentifier) const
{
  return select_union_from_part_ordinals(meta(), get_negative_levelset_block_ordinals(levelSetIdentifier));
}

std::vector<unsigned> Phase_Support::get_levelset_decomposed_block_ordinals(const Surface_Identifier levelSetIdentifier) const
{
  std::vector<unsigned> levelsetBlockOrdinals;
  for (auto * partPtr : meta().get_mesh_parts())
  {
    if (partPtr->primary_entity_rank() == stk::topology::ELEMENT_RANK &&
        (is_decomposed(*partPtr) && level_set_is_used_by_nonconforming_part(levelSetIdentifier, get_nonconforming_part(*partPtr))))
        levelsetBlockOrdinals.push_back(partPtr->mesh_meta_data_ordinal());
  }
  return levelsetBlockOrdinals;
}

stk::mesh::Selector
Phase_Support::get_levelset_decomposed_blocks_selector(const Surface_Identifier levelSetIdentifier) const
{
  return select_union_from_part_ordinals(meta(), get_levelset_decomposed_block_ordinals(levelSetIdentifier));
}

Phase_Support::Phase_Support(stk::mesh::MetaData & meta)
: myMeta(&meta),
  myAuxMeta(&AuxMetaData::get(meta)),
  oneLevelSetPerPhase(false) {}

Phase_Support &
Phase_Support::get(const stk::mesh::MetaData & meta)
{
  Phase_Support * support = const_cast<Phase_Support *>(meta.get_attribute<Phase_Support>());
  STK_ThrowRequireMsg(nullptr != support, "No Phase_Support found for MetaData.");
  return *support;
}

Phase_Support &
Phase_Support::get(stk::mesh::MetaData & meta)
{
  Phase_Support * support = const_cast<Phase_Support *>(meta.get_attribute<Phase_Support>());
  if (nullptr == support)
  {
    support = new Phase_Support(meta);
    meta.declare_attribute_with_delete<Phase_Support>(support);
  }
  return *support;
}

void
Phase_Support::set_phase_input(const PhaseInput & phaseInput)
{
  STK_ThrowRequireMsg(myMeshPhases.empty() && myMeshBlockPhasesByName.empty(), "Phase_Support::set_phase_input called after phases were already set.");
  myMeshPhases = phaseInput.meshPhases;
  myMeshBlockPhasesByName = phaseInput.blockPhasesByName;
  oneLevelSetPerPhase = phaseInput.oneLevelSetPerPhase;
}

bool
Phase_Support::exists_and_has_phases_defined(const stk::mesh::MetaData & meta)
{
  Phase_Support * support = const_cast<Phase_Support *>(meta.get_attribute<Phase_Support>());
  return support != nullptr && support->phases_defined();
}

stk::mesh::Selector Phase_Support::get_all_interface_surfaces_selector() const
{
  std::vector<stk::mesh::Part *> interfaceParts;
  for (auto * part : meta().get_mesh_parts())
    if(is_interface(*part))
      interfaceParts.push_back(part);
  return selectUnion(interfaceParts);
}

std::string
Phase_Support::get_phase_part_name(const stk::mesh::Part & originalPart, const NamedPhase & phase) const
{
  // The alive phase for death problems (empty name) just uses the original part itself.
  const std::string suffix = phase.name().empty() ? "" : "_" + phase.name();
  return build_part_name(*this, originalPart, phase.tag(), suffix);
}

stk::mesh::Part &
Phase_Support::declare_phase_part(const stk::mesh::Part & originalPart, const NamedPhase & phase)
{
  stk::mesh::Part & phasePart = aux_meta().declare_io_part(get_phase_part_name(originalPart, phase), originalPart.primary_entity_rank());
  STK_ThrowAssert(!phase.name().empty() || phasePart.mesh_meta_data_ordinal() == originalPart.mesh_meta_data_ordinal());
  if (stk::topology::INVALID_TOPOLOGY != originalPart.topology())
    stk::mesh::set_topology(phasePart, originalPart.topology());
  return phasePart;
}

void
Phase_Support::add_conforming_phase_part(stk::mesh::Part & originalPart, const NamedPhase & phase)
{
  stk::mesh::Part & conformingPart = declare_phase_part(originalPart, phase);
  stk::mesh::Part & nonconformingPart = get_nonconforming_part_for_original_part(originalPart);

  if(krinolog.shouldPrint(LOG_PARTS))
    krinolog << "Adding conforming part " << conformingPart.name() << " with topology " << conformingPart.topology().name()
             << ", nonconforming part " << nonconformingPart.name() << ", original part " << originalPart.name()
             << " " << phase.tag() << stk::diag::dendl;

  myPhasePartInfo.setup_conforming_part(meta().get_parts().size(), conformingPart.mesh_meta_data_ordinal(), nonconformingPart.mesh_meta_data_ordinal(), originalPart.mesh_meta_data_ordinal(), phase.tag());
  // Conforming parts are looked up (build_part_name, add_interface_part) while more are still being added.
  myPhasePartInfo.sort_conforming_parts();
}

std::string
Phase_Support::get_nonconforming_part_name(const stk::mesh::Part & originalPart) const
{
  static const std::string nonconformingSuffix = "_nonconformal";
  return build_part_name(*this, originalPart, PhaseTag(), nonconformingSuffix);
}

stk::mesh::Part &
Phase_Support::get_nonconforming_part_for_original_part(const stk::mesh::Part & originalPart) const
{
  // By name, since this is used before the nonconforming part is registered in PhasePartInfo.
  return aux_meta().get_part(get_nonconforming_part_name(originalPart));
}

void
Phase_Support::create_nonconforming_parts(const PartSet & decomposedIoParts)
{
  for(auto && decomposedIoPartOrd : decomposedIoParts)
  {
    const stk::mesh::Part & ioPart = aux_meta().get_part(meta().get_part(decomposedIoPartOrd).name());

    const std::string nonconformingPartName = get_nonconforming_part_name(ioPart);

    if(!aux_meta().has_part(nonconformingPartName))
    {
      const bool restartOnlyIOPart = true;
      if (ioPart.topology() == stk::topology::INVALID_TOPOLOGY)
        aux_meta().declare_io_part(nonconformingPartName, ioPart.primary_entity_rank(), restartOnlyIOPart);
      else
        aux_meta().declare_io_part_with_topology(nonconformingPartName, ioPart.topology(), restartOnlyIOPart);
      if(krinolog.shouldPrint(LOG_PARTS)) krinolog << "Created nonconformal part " << nonconformingPartName << stk::diag::dendl;
    }
  }
}

stk::mesh::PartVector
Phase_Support::get_nonconformal_parts() const
{
  stk::mesh::PartVector result;
  for (auto * part : meta().get_mesh_parts())
    if (myPhasePartInfo.is_nonconforming_part(part->mesh_meta_data_ordinal()))
      result.push_back(&meta().get_part(part->mesh_meta_data_ordinal()));
  return result;
}

stk::mesh::PartVector
Phase_Support::get_nonconformal_parts_of_rank(const stk::mesh::EntityRank rank) const
{
  stk::mesh::PartVector result;
  for (auto & part : meta().get_mesh_parts())
    if (part->primary_entity_rank() == rank && myPhasePartInfo.is_nonconforming_part(part->mesh_meta_data_ordinal()))
      result.push_back(&meta().get_part(part->mesh_meta_data_ordinal()));
  return result;
}

stk::mesh::PartVector
Phase_Support::get_conformal_parts(const bool includeInterfaceParts) const
{
  stk::mesh::PartVector result;
  for (auto & part : meta().get_mesh_parts())
    if (is_conformal(*part) ||  (includeInterfaceParts && is_interface(*part)))
      result.push_back(part);
  return result;
}

stk::mesh::PartVector
Phase_Support::get_conformal_parts_of_rank(const stk::mesh::EntityRank rank, const bool includeInterfaceParts) const
{
  stk::mesh::PartVector result;
  for (auto & part : meta().get_mesh_parts())
    if (part->primary_entity_rank() == rank &&
        (is_conformal(*part) || (includeInterfaceParts && is_interface(*part))))
      result.push_back(part);
  return result;
}

stk::mesh::PartVector
Phase_Support::get_conformal_parts_of_rank(const stk::mesh::EntityRank rank) const
{
  STK_ThrowRequireMsg(rank != meta().side_rank(), "Do not call get_conformal_parts_of_rank without specifying whether to includeInterfaceParts for side rank.");
  return get_conformal_parts_of_rank(rank, false);
}

void
Phase_Support::update_touching_parts_for_phase_part(const stk::mesh::Part & origPart, const stk::mesh::Part & phasePart, const PhaseTag & phase)
{
  const bool isNonconforming = phase.empty();
  const std::set<stk::mesh::PartOrdinal> & origTouchingBlockOrdinals = my_input_block_surface_connectivity.get_blocks_touching_surface(origPart.mesh_meta_data_ordinal());

  std::vector<const stk::mesh::Part*> phaseTouchingBlocks = meta().get_blocks_touching_surface(&phasePart);

  for (auto && origTouchingBlockOrdinal : origTouchingBlockOrdinals)
  {
    stk::mesh::Part & origTouchingBlock = meta().get_part(origTouchingBlockOrdinal);
    if (isNonconforming)
    {
      if (is_decomposed(origTouchingBlock))
        insert_into_touching_blocks(get_nonconforming_part(origTouchingBlock), phaseTouchingBlocks);
      else
        insert_into_touching_blocks(origTouchingBlock, phaseTouchingBlocks);
    }
    else
    {
      stk::mesh::Part * phaseTouchingBlock = find_conforming_part_if_decomposed_by_phase(origTouchingBlock, phase);
      if (nullptr != phaseTouchingBlock)
        insert_into_touching_blocks(*phaseTouchingBlock, phaseTouchingBlocks);
      else
        insert_into_touching_blocks(origTouchingBlock, phaseTouchingBlocks);
    }
  }

  if(krinolog.shouldPrint(LOG_PARTS))
  {
    const std::string conformingType = isNonconforming ? "Nonconforming" : "Conforming";
    krinolog << conformingType << " surface " << phasePart.name() << " touches blocks ";
    for (auto && phaseTouchingBlock : phaseTouchingBlocks)
      krinolog << phaseTouchingBlock->name() << " ";
    krinolog << "\n";
  }

  meta().set_surface_to_block_mapping(&phasePart, phaseTouchingBlocks);
}

void
Phase_Support::build_decomposed_block_surface_connectivity()
{
  std::set<std::pair<unsigned, unsigned>> nonconformingAndOriginalPartOrdinalPairs;

  for (auto && part : meta().get_mesh_parts())
  {
    const unsigned partOrd = part->mesh_meta_data_ordinal();
    if (part->primary_entity_rank() == meta().side_rank())
    {
      if (myPhasePartInfo.is_conforming_part(partOrd))
      {
        stk::mesh::Part & origPart = find_original_part(*part);
        stk::mesh::Part & nonconformingPart = get_nonconforming_part(*part);
        update_touching_parts_for_phase_part(origPart, *part, myPhasePartInfo.get_conforming_part_phase(partOrd));

        // store off nonconforming and original parts for second pass below
        nonconformingAndOriginalPartOrdinalPairs.emplace(nonconformingPart.mesh_meta_data_ordinal(), origPart.mesh_meta_data_ordinal());
      }
      else if (myPhasePartInfo.is_interface_part(partOrd))
      {
        const stk::mesh::Part & touchingBlock = meta().get_part( myPhasePartInfo.get_interface_part_touching_part(partOrd));
        if(krinolog.shouldPrint(LOG_PARTS)) krinolog << "Interface surface " << part->name() << " touches block " << touchingBlock.name() << "\n";
        std::vector<const stk::mesh::Part*> touching_blocks = meta().get_blocks_touching_surface(part);
        if (std::find(touching_blocks.begin(), touching_blocks.end(), &touchingBlock) == touching_blocks.end())
          touching_blocks.push_back(&touchingBlock);

        meta().set_surface_to_block_mapping(part, touching_blocks);
      }
    }
  }

  const PhaseTag emptyPhaseToIndicateNonconformingPart;
  for (auto && nonconformingAndOriginalPartOrdinalPair : nonconformingAndOriginalPartOrdinalPairs)
  {
    const stk::mesh::Part & nonconformingPart = meta().get_part(nonconformingAndOriginalPartOrdinalPair.first);
    const stk::mesh::Part & origPart = meta().get_part(nonconformingAndOriginalPartOrdinalPair.second);
    update_touching_parts_for_phase_part(origPart, nonconformingPart, emptyPhaseToIndicateNonconformingPart);
  }

  if(krinolog.shouldPrint(LOG_PARTS))
  {
    Block_Surface_Connectivity::dump_surface_connectivity(meta());
  }
}

std::map<int,int> Phase_Support::build_phase_conversion_map(const PhaseTag & fromPhase, const PhaseTag & toPhase) const
{
  return build_phase_conversion_map_for_element_and_side_parts(meta(), myPhasePartInfo, fromPhase, toPhase);
}

stk::mesh::Selector Phase_Support::get_selector_of_blocks_with_phase(const PhaseTag & phase) const
{
  std::vector<unsigned> phaseBlockOrdinals;
  for (auto * part : meta().get_mesh_parts())
  {
    const unsigned partOrd = part->mesh_meta_data_ordinal();
      if (part->primary_entity_rank() == stk::topology::ELEMENT_RANK &&
          myPhasePartInfo.is_conforming_part(partOrd) &&
          myPhasePartInfo.get_conforming_part_phase(partOrd) == phase)
        phaseBlockOrdinals.push_back(partOrd);
  }
  return select_union_from_part_ordinals(meta(), phaseBlockOrdinals);
}

void
Phase_Support::create_interface_phase_parts(
    const PhaseVec& ls_phases,
    const PartSet& decomposed_ioparts,
    const Interface_Name_Generator& interface_name_gen)
{
  const int num_ls_phases = ls_phases.size();
  if (num_ls_phases > 1)
  {
    const int num_interfaces = num_ls_phases * (num_ls_phases - 1) / 2;
    std::vector<std::pair<unsigned, unsigned> > interface_phases;
    interface_phases.reserve(num_interfaces);
    for (int i = 0; i < num_ls_phases; ++i)
    {
      for (int j = i+1; j < num_ls_phases; ++j)
      {
        std::pair<unsigned, unsigned> phase_pair(i, j);
        interface_phases.push_back(phase_pair);
      }
    }
    STK_ThrowRequire((int) interface_phases.size() == num_interfaces);

    for (int i = 0; i < num_interfaces; ++i)
    {
      const int name_compare = ls_phases[interface_phases[i].first].name().compare(ls_phases[interface_phases[i].second].name());
      STK_ThrowRequire(name_compare != 0);
      // form phase intersection
      const auto & phase0 = (name_compare < 0) ? ls_phases[interface_phases[i].first] : ls_phases[interface_phases[i].second];
      const auto & phase1 = (name_compare < 0) ? ls_phases[interface_phases[i].second] : ls_phases[interface_phases[i].first];
      std::string superset_interface_phase_name = phase0.name() + "_" + phase1.name();
      const std::string superset_phase_part_name = interface_name_gen.interface_superset_name(superset_interface_phase_name);

      // assumes all sides have same topology
      stk::mesh::Part &superset_phase_part = aux_meta().declare_io_part(superset_phase_part_name, meta().side_rank());

      if(krinolog.shouldPrint(LOG_PARTS)) krinolog << "Adding interfacial phase part: conformal_part_name=" << superset_phase_part.name() << "\n";

      // Each phase pair defines an interface. When the generator gives every phase pairing
      // its own interface part we create both orientations; when it reuses a single interface
      // part for both orientations (e.g. cdfem death) we create only the first, since the
      // swapped orientation would re-register the same part.
      std::vector<std::pair<NamedPhase,NamedPhase>> orientation_pairs = {std::make_pair(phase0,phase1)};
      if (interface_name_gen.has_part_per_phase())
        orientation_pairs.push_back(std::make_pair(phase1,phase0));

      for (auto && phase_pair : orientation_pairs)
      {
        const auto & touching_phase = phase_pair.first.tag();
        const auto & opposite_phase = phase_pair.second.tag();

        add_interface_superset_part(superset_phase_part);

        for (auto && ioPartOrd : decomposed_ioparts)
        {
          const stk::mesh::Part & ioPart = meta().get_part(ioPartOrd);
          // only handle blocks for now -> creating interfacial surface "phases"
          if (ioPart.primary_entity_rank() == stk::topology::ELEMENT_RANK)
          {
            std::string subset_interface_phase_name = phase_pair.first.name() + "_" + phase_pair.second.name();
            const std::string subset_phase_part_name = interface_name_gen.interface_subset_name(ioPart.name(), subset_interface_phase_name);
            const bool restartOnlyIOPart = true; // Subset part is only output for restart
            stk::mesh::Part &subset_phase_part = aux_meta().declare_io_part_with_topology(subset_phase_part_name, interface_topology(ioPart.topology()), restartOnlyIOPart);
            meta().declare_part_subset(superset_phase_part, subset_phase_part);

            const stk::mesh::Part & nonconformingPart = get_nonconforming_part(ioPart);

            if(krinolog.shouldPrint(LOG_PARTS))
            {
              krinolog << "  Adding subset interfacial phase part: conformal_part_name=" << subset_phase_part.name()
                  << ", original_part_name=" << ioPart.name();
              krinolog << ", touching_phase=" << touching_phase
                       << ", opposite_phase=" << opposite_phase;
              krinolog << "\n";
            }

            add_interface_part(subset_phase_part, nonconformingPart, touching_phase, opposite_phase);
          }
        }
      }
    }
  }
  myPhasePartInfo.sort_interface_parts();
}

void Phase_Support::add_interface_part(const stk::mesh::Part & interfacePart,
    const stk::mesh::Part & nonconformingPart,
    const PhaseTag & touchingPhase,
    const PhaseTag & oppositePhase)
{
  const unsigned touchingPartOrdinal = get_conforming_part(nonconformingPart, touchingPhase).mesh_meta_data_ordinal();
  const unsigned oppositePartOrdinal = get_conforming_part(nonconformingPart, oppositePhase).mesh_meta_data_ordinal();

  myPhasePartInfo.setup_interface_part(meta().get_parts().size(), interfacePart.mesh_meta_data_ordinal(), touchingPartOrdinal, oppositePartOrdinal);
}

void Phase_Support::add_interface_superset_part(const stk::mesh::Part & interfacePart)
{
  myPhasePartInfo.setup_interface_superset_part(meta().get_parts().size(), interfacePart.mesh_meta_data_ordinal());
}

const NamedPhase & Phase_Support::get_phase_with_name(const std::string & phaseName) const
{
  const auto iter = std::find_if(myMeshPhases.begin(), myMeshPhases.end(),
      [&phaseName](const NamedPhase & x) { return phaseName == x.name(); });
  STK_ThrowRequireMsg(iter != myMeshPhases.end(), "Did not find phase " << phaseName);
  return *iter;
}

void Phase_Support::set_small_cluster_cleanup(const std::vector<std::tuple<std::string,std::string,unsigned>> & smallClusterCleanup)
{
  std::vector<std::tuple<NamedPhase,NamedPhase,unsigned>> fromPhaseToPhaseAndClusterSize;
  for (auto & [fromPhaseName,toPhaseName,clusterSize] : smallClusterCleanup)
  {
    mySmallClusterCleanup.emplace_back(get_phase_with_name(fromPhaseName), get_phase_with_name(toPhaseName), clusterSize);
  }
}

std::set<stk::mesh::PartOrdinal>
Phase_Support::get_blocks_touching_sideset(const stk::mesh::Part & sideset) const
{
  std::set<stk::mesh::PartOrdinal> touchingBlocks = my_input_block_surface_connectivity.get_blocks_touching_surface(sideset.mesh_meta_data_ordinal());
  for (auto * subset : sideset.subsets())
  {
    const std::set<stk::mesh::PartOrdinal> subsetTouchingBlocks = my_input_block_surface_connectivity.get_blocks_touching_surface(subset->mesh_meta_data_ordinal());
    touchingBlocks.insert(subsetTouchingBlocks.begin(), subsetTouchingBlocks.end());
  }
  return touchingBlocks;
}

std::vector<unsigned>
Phase_Support::get_sidesets_touching_decomposed_blocks(const DecompositionPackage & decomps) const
{
  std::set<unsigned> sidesetOrds;
  for (unsigned i=0; i<decomps.size(); ++i)
  {
    for (auto * block : decomps.get_blocks_to_decompose(i))
    {
      for (auto && surfOrd : my_input_block_surface_connectivity.get_surfaces_touching_block(block->mesh_meta_data_ordinal()))
      {
        stk::mesh::Part & surf = meta().get_part(surfOrd);
        stk::mesh::Part * parentSideset = get_parent_sideset(surf);
        sidesetOrds.insert((parentSideset == nullptr) ? surfOrd : parentSideset->mesh_meta_data_ordinal());
      }
    }
  }
  return std::vector<unsigned>(sidesetOrds.begin(), sidesetOrds.end());
}

std::vector<unsigned>
Phase_Support::find_sidesets_with_conflicting_phase_criteria(const DecompositionPackage & decomps) const
{
  std::vector<unsigned> conflictingSidesetOrds;
  for (auto && sidesetOrd : get_sidesets_touching_decomposed_blocks(decomps))
    if (!get_conflicting_phase_names(decomps, get_blocks_touching_sideset(meta().get_part(sidesetOrd))).empty())
      conflictingSidesetOrds.push_back(sidesetOrd);
  return conflictingSidesetOrds;
}

bool
Phase_Support::validate_surface_phase_connectivity(const DecompositionPackage & decomps) const
{
  // The phase parts of a sideset (e.g. surface_1_dead and its subsets) are shared by every block the
  // sideset touches, so all of those blocks must use the same criterion (phase tag) for a given phase name.
  const std::vector<unsigned> conflictingSidesetOrds = find_sidesets_with_conflicting_phase_criteria(decomps);
  for (auto && sidesetOrd : conflictingSidesetOrds)
  {
    const stk::mesh::Part & sideset = meta().get_part(sidesetOrd);
    const std::set<stk::mesh::PartOrdinal> touchingBlocks = get_blocks_touching_sideset(sideset);

    std::ostringstream msg;
    msg << "Sideset " << sideset.name() << " touches blocks that are decomposed with different criteria for a phase with the same name.\n";
    std::set<stk::mesh::PartOrdinal> decomposedTouchingBlocks;
    for (auto && [phaseName, tagToBlocks] : get_conflicting_phase_names(decomps, touchingBlocks))
    {
      msg << "  Phase " << describe_phase_name(phaseName) << " is defined differently on these blocks touched by " << sideset.name() << ":\n";
      unsigned criterion = 0;
      for (auto && [tag, blockDecomps] : tagToBlocks)
      {
        msg << "    criterion " << ++criterion << " (" << decomps.get_description(blockDecomps.front().decompIndex) << "): blocks";
        for (auto && blockDecomp : blockDecomps)
        {
          msg << " " << meta().get_part(blockDecomp.block).name();
          decomposedTouchingBlocks.insert(blockDecomp.block);
        }
        msg << "\n";
      }
    }

    std::vector<std::string> undecomposedTouchingBlockNames;
    for (auto && blockOrd : touchingBlocks)
      if (decomposedTouchingBlocks.count(blockOrd) == 0)
        undecomposedTouchingBlockNames.push_back(meta().get_part(blockOrd).name());
    if (!undecomposedTouchingBlockNames.empty())
    {
      msg << "  " << sideset.name() << " also touches blocks that are not decomposed:";
      for (auto && blockName : undecomposedTouchingBlockNames)
        msg << " " << blockName;
      msg << "\n";
    }

    msg << "Each phase of a sideset gets a single set of mesh parts (for example " << sideset.name()
        << " and " << sideset.name() << "_<phase>) shared by every block the sideset touches, so the phase of a side cannot depend on which block it is on.\n"
        << "To fix this, change the mesh so that " << sideset.name()
        << " is replaced by separate sidesets, each touching only blocks that use the same criterion"
        << " (for example, one sideset on the faces of each block listed above), and update any input commands that reference "
        << sideset.name() << " to use the new sidesets.\n";

    stk::RuntimeDoomedAdHoc() << msg.str();
  }
  return conflictingSidesetOrds.empty();
}

void
Phase_Support::create_phase_parts(const PhaseVec& phases, const PartSet& decomposedParts)
{
  for (auto && phase : phases)
  {
    for (auto && partOrd : decomposedParts)
    {
      stk::mesh::Part & part = meta().get_part(partOrd);
      if(get_parent_sideset(part) == nullptr)
        add_conforming_phase_part(part, phase);
    }
  }
}

void
Phase_Support::subset_and_alias_surface_phase_parts(const PhaseVec& phases,
    const PartSet& decomposedParts)
{
  for (auto && partOrd : decomposedParts)
  {
    stk::mesh::Part & origPart = meta().get_part(partOrd);
    if (!(origPart.subsets().empty()))
    {
      for (auto && phase : phases)
      {
        stk::mesh::Part & conformingSuperset = get_conforming_part(origPart, phase.tag());
        stk::mesh::Part & nonconformingSuperset = get_nonconforming_part(origPart);

        for (auto * partSubset : origPart.subsets())
        {
          STK_ThrowRequire(nullptr != partSubset);

          add_conforming_phase_part(*partSubset, phase);
          stk::mesh::Part & conformingSubset = get_conforming_part(*partSubset, phase.tag());

          if(krinolog.shouldPrint(LOG_PARTS)) krinolog << "Adding " << conformingSubset.name() << " as subset of " << conformingSuperset.name() << stk::diag::dendl;
          meta().declare_part_subset(conformingSuperset, conformingSubset);

          stk::mesh::Part & nonconformingSubset = get_nonconforming_part(*partSubset);

          if(krinolog.shouldPrint(LOG_PARTS)) krinolog << "Adding " << nonconformingSubset.name() << " as subset of " << nonconformingSuperset.name() << stk::diag::dendl;
          meta().declare_part_subset(nonconformingSuperset, nonconformingSubset);

          for (auto && touchingBlockOrdinal : my_input_block_surface_connectivity.get_blocks_touching_surface(partSubset->mesh_meta_data_ordinal()))
          {
            const std::string conformingPartAlias = conformingSuperset.name() + "_" + meta().get_part(touchingBlockOrdinal).name();
            if(krinolog.shouldPrint(LOG_PARTS)) krinolog << "Adding alias " << conformingPartAlias << " for conforming part " << conformingSubset.name() << stk::diag::dendl;
            aux_meta().define_part_alias(conformingSubset, conformingPartAlias);
          }
        }
      }
    }
  }
}

void Phase_Support::decompose_blocks(const DecompositionPackage & decomps)
{
  for(unsigned i=0; i<decomps.size(); ++i)
    if(!do_blocks_have_topology_that_can_be_decomposed(decomps.get_blocks_to_decompose(i)))
      return;

  const std::vector<unsigned> conflictingSidesetOrds = find_sidesets_with_conflicting_phase_criteria(decomps);
  STK_ThrowRequireMsg(conflictingSidesetOrds.empty(), "Sideset " << meta().get_part(conflictingSidesetOrds.front()).name()
      << " touches blocks decomposed with different criteria for a phase with the same name.  Call validate_surface_phase_connectivity at parse time for a detailed error.");

  std::vector<PartSet> blocksAndSurfaces(decomps.size());
  for(unsigned i=0; i<decomps.size(); ++i)
  {
    const PhaseVec & phases = decomps.get_phases(i);

    blocksAndSurfaces[i] = get_blocks_and_touching_surfaces(meta(), decomps.get_blocks_to_decompose(i), my_input_block_surface_connectivity);

    create_nonconforming_parts(blocksAndSurfaces[i]);

    create_phase_parts(phases, blocksAndSurfaces[i]);
  }

  // Second pass: subset_and_alias_surface_phase_parts names surface subset phase parts after the
  // conforming blocks they touch, so every decomposition's block phase parts must already exist.
  for(unsigned i=0; i<decomps.size(); ++i)
  {
    const PhaseVec & phases = decomps.get_phases(i);

    subset_and_alias_surface_phase_parts(phases, blocksAndSurfaces[i]);

    create_interface_phase_parts(phases, blocksAndSurfaces[i], decomps.get_interface_name_generator(i));
  }

  set_all_decomposed_blocks_selector();
}

//--------------------------------------------------------------------------------

void
Phase_Support::determine_block_phases(const std::set<std::string> & FEmodel_block_names)
{
  if (myMeshPhases.empty()) return;

  STK_ThrowAssertMsg(my_mesh_block_phases.empty(), "determine_block_phases should only be called once per mesh");

  for (auto && part : meta().get_parts())
  {
    if (part->primary_entity_rank() == stk::topology::ELEMENT_RANK)
    {
      auto & block_phases = my_mesh_block_phases[part->mesh_meta_data_ordinal()];
      for(unsigned phase_index=0; phase_index<myMeshPhases.size(); ++phase_index)
      {
        const std::string & phase_name = myMeshPhases[phase_index].name();
        std::string phase_part_name = part->name() + "_" + phase_name;
        std::transform(phase_part_name.begin(), phase_part_name.end(), phase_part_name.begin(), ::tolower);
        if (FEmodel_block_names.find(phase_part_name) != FEmodel_block_names.end() )
        {
          block_phases.insert(phase_index);
        }
      }
    }
  }
}

//--------------------------------------------------------------------------------

static std::vector<unsigned> get_matching_phase_indices(
    const std::string & phaseNameOrPattern,
    const std::string & blockName,
    const std::map<std::string, unsigned> & phaseNameToIndexMap)
{
  std::vector<unsigned> matchingPhaseIndices;

  if (!utility::has_wildcard(phaseNameOrPattern))
  {
    const auto iter = phaseNameToIndexMap.find(phaseNameOrPattern);

    if (iter != phaseNameToIndexMap.end())
    {
      matchingPhaseIndices.push_back(iter->second);
    }
    else
    {
      stk::RuntimeDoomedAdHoc()
          << "Failed to find phase " << phaseNameOrPattern
          << ", specified for block " << blockName
          << " in the list of all phases.\n";
    }

    return matchingPhaseIndices;
  }

  for (const auto & phaseEntry : phaseNameToIndexMap)
  {
    const std::string & phaseName = phaseEntry.first;
    const unsigned phaseIndex = phaseEntry.second;

    if (utility::is_wildcard_match(phaseNameOrPattern, phaseName))
      matchingPhaseIndices.push_back(phaseIndex);
  }

  if (matchingPhaseIndices.empty())
  {
    stk::RuntimeDoomedAdHoc()
        << "Failed to find any phases matching pattern " << phaseNameOrPattern
        << ", specified for block " << blockName
        << " in the list of all phases.\n";
  }

  return matchingPhaseIndices;
}

static bool is_supported_name_pattern(const std::string & nameOrPattern, const std::string & description)
{
  if (utility::is_supported_wildcard_pattern(nameOrPattern))
    return true;

  stk::RuntimeDoomedAdHoc()
      << "The " << description << " " << nameOrPattern
      << " in the finite element model specification is not supported.  Names may contain at most one '*' wildcard, and may not contain '?'.\n";
  return false;
}

static bool are_supported_name_patterns(const std::string & blockNameOrPattern, const std::vector<std::string> & phaseNamesOrPatterns)
{
  bool areSupported = is_supported_name_pattern(blockNameOrPattern, "part name");
  for (const std::string & phaseNameOrPattern : phaseNamesOrPatterns)
    areSupported &= is_supported_name_pattern(phaseNameOrPattern, "subdomain name");
  return areSupported;
}

void
Phase_Support::determine_block_phases()
{
  if (myMeshPhases.empty()) return;

  STK_ThrowAssertMsg(
      my_mesh_block_phases.empty(),
      "determine_block_phases should only be called once per mesh");

  if (myMeshBlockPhasesByName.empty())
  {
    stk::RuntimeDoomedAdHoc()
        << "Failure in determine_block_phases. Were the subdomains specified "
        << "for the decomposed blocks in the finite element model specification?\n";
    return;
  }

  std::map<std::string, unsigned> phaseNameToIndexMap;
  for (unsigned phaseIndex = 0; phaseIndex < myMeshPhases.size(); ++phaseIndex)
  {
    const std::string & phaseName = myMeshPhases[phaseIndex].name();
    phaseNameToIndexMap[phaseName] = phaseIndex;
  }

  for (const auto & mapEntry : myMeshBlockPhasesByName)
  {
    const std::string & blockNameOrPattern = mapEntry.first;
    const auto & phaseNamesOrPatterns = mapEntry.second;

    if (!are_supported_name_patterns(blockNameOrPattern, phaseNamesOrPatterns))
      continue;

    const std::vector<const stk::mesh::Part *> matchingBlocks =
        aux_meta().get_matching_io_parts_of_rank(blockNameOrPattern, stk::topology::ELEMENT_RANK);

    if (matchingBlocks.empty())
      stk::RuntimeDoomedAdHoc() << "The part " << blockNameOrPattern << " is listed in the finite element model specification, but does not match any element blocks in the input mesh.\n";

    for (const stk::mesh::Part * block : matchingBlocks)
    {
      STK_ThrowAssert(block != nullptr);

      for (const std::string & phaseNameOrPattern : phaseNamesOrPatterns)
      {
        const std::vector<unsigned> matchingPhaseIndices =
            get_matching_phase_indices(
                phaseNameOrPattern,
                block->name(),
                phaseNameToIndexMap);

        auto & blockPhases = my_mesh_block_phases[block->mesh_meta_data_ordinal()];

        for (const unsigned phaseIndex : matchingPhaseIndices)
          blockPhases.insert(phaseIndex);
      }
    }
  }
}

//--------------------------------------------------------------------------------

std::vector<stk::mesh::Part*>
Phase_Support::get_blocks_decomposed_by_surface(const std::vector<unsigned> & surfacePhases) const
{
  std::vector<stk::mesh::Part*> blocksDecomposedBySurface;
  for (const auto & entry : my_mesh_block_phases)
  {
    const auto & blockPhases = entry.second;
    if (does_block_have_any_of_surface_phases(blockPhases, surfacePhases))
    {
      stk::mesh::Part & block = meta().get_part(entry.first);
      blocksDecomposedBySurface.push_back(&block);
    }
  }
  return blocksDecomposedBySurface;
}

//--------------------------------------------------------------------------------

void
Phase_Support::setup_phases()
{
  DecompositionPackage decomps;
  for (auto && [blockOrdinal, blockPhaseIndices] : my_mesh_block_phases)
  {
    PhaseVec blockPhases;
    for (auto && phaseIndex : blockPhaseIndices)
      blockPhases.push_back(myMeshPhases[phaseIndex]);

    decomps.add_levelset_decomposition(stk::mesh::PartVector{&meta().get_part(blockOrdinal)}, blockPhases);
  }

  if (!validate_surface_phase_connectivity(decomps))
    return;

  decompose_blocks(decomps);
  build_decomposed_block_surface_connectivity();
}

//--------------------------------------------------------------------------------

stk::mesh::Part *
Phase_Support::find_conforming_part_if_decomposed_by_phase(const stk::mesh::Part & part, const PhaseTag & phase) const
{
  unsigned conformingPartOrdinal = 0;
  if (myPhasePartInfo.find_conforming_part(
      part.mesh_meta_data_ordinal(),
      phase,
      conformingPartOrdinal))
  {
    return &meta().get_part(conformingPartOrdinal);
  }

  return nullptr;
}

stk::mesh::Part &
Phase_Support::get_conforming_part(const stk::mesh::Part & part, const PhaseTag & phase) const
{
  stk::mesh::Part * conformingPartPtr = find_conforming_part_if_decomposed_by_phase(part, phase);
  STK_ThrowRequireMsg(nullptr != conformingPartPtr, "Part " << part.name() << " is not decomposed by " << phase);
  return *conformingPartPtr;
}

//--------------------------------------------------------------------------------

bool
Phase_Support::is_conformal(const stk::mesh::Part & part) const
{
  return myPhasePartInfo.is_conforming_part(part.mesh_meta_data_ordinal());
}

//--------------------------------------------------------------------------------

bool
Phase_Support::is_nonconformal(const stk::mesh::Part & part) const
{
  return myPhasePartInfo.is_nonconforming_part(part.mesh_meta_data_ordinal());
}

//--------------------------------------------------------------------------------

bool
Phase_Support::is_interface(const stk::mesh::Part & part) const
{
  return myPhasePartInfo.is_interface_part(part.mesh_meta_data_ordinal()) || myPhasePartInfo.is_interface_superset_part(part.mesh_meta_data_ordinal());
}

//--------------------------------------------------------------------------------

stk::mesh::Part & Phase_Support::get_nonconforming_part(const stk::mesh::Part & part) const
{
  // Will throw if not decomposed or interface part
  return meta().get_part(myPhasePartInfo.get_nonconforming_part(part.mesh_meta_data_ordinal()));
}

//--------------------------------------------------------------------------------

stk::mesh::Part & Phase_Support::find_original_part(const stk::mesh::Part & part) const
{
  STK_ThrowAssertMsg(!myPhasePartInfo.is_interface_superset_part(part.mesh_meta_data_ordinal()), "Should not call find_original_part with interface superset part.");

  if (myPhasePartInfo.is_decomposed(part.mesh_meta_data_ordinal()) || myPhasePartInfo.is_interface_part(part.mesh_meta_data_ordinal()))
    return meta().get_part(myPhasePartInfo.get_original_part(part.mesh_meta_data_ordinal()));
  return meta().get_part(part.mesh_meta_data_ordinal()); // hidden const_cast
}

//--------------------------------------------------------------------------------

const stk::mesh::Part *
Phase_Support::find_interface_part(const stk::mesh::Part & vol0, const stk::mesh::Part & vol1) const
{
  unsigned interfacePartOrd = 0;
  if (myPhasePartInfo.find_interface_part(
          vol0.mesh_meta_data_ordinal(),
          vol1.mesh_meta_data_ordinal(),
          interfacePartOrd))
  {
    return &meta().get_part(interfacePartOrd);
  }

  return nullptr;
}

//--------------------------------------------------------------------------------

const PhaseTag &
Phase_Support::get_iopart_phase(const stk::mesh::Part & part) const
{
  return myPhasePartInfo.get_conforming_part_phase(part.mesh_meta_data_ordinal());
}

void Phase_Support::set_all_decomposed_blocks_selector()
{
  stk::mesh::PartVector decomposedParts;
  for (auto * part : meta().get_mesh_parts())
    if (part->primary_entity_rank() == stk::topology::ELEMENT_RANK && is_decomposed(*part))
      decomposedParts.push_back(part);
  all_decomposed_blocks_selector = stk::mesh::selectUnion(decomposedParts);
}

bool Phase_Support::is_decomposed(const stk::mesh::Part & part) const
{
  // This includes original, conforming, and nonconforming parts (but not interface parts)
  return myPhasePartInfo.is_decomposed(part.mesh_meta_data_ordinal());
}

//--------------------------------------------------------------------------------

void Phase_Support::register_blocks_for_level_set(const Surface_Identifier levelSetIdentifier,
                                                  const std::vector<stk::mesh::Part *> & blocks_decomposed_by_ls)
{
  for (auto * blockPtr : blocks_decomposed_by_ls)
    myLsUsedByParts[levelSetIdentifier].insert(blockPtr);
}

void Phase_Support::register_blocks_for_surface(const Surface_Identifier & surfID)
{
  const std::vector<unsigned> surfPhases = get_indices_of_phases_that_depend_on_surface(has_one_levelset_per_phase(), myMeshPhases, surfID);
  const std::vector<stk::mesh::Part*> decomposedBlocks = get_blocks_decomposed_by_surface(surfPhases);
  register_blocks_for_level_set(surfID, decomposedBlocks);
}

//--------------------------------------------------------------------------------
bool Phase_Support::level_set_is_used_by_nonconforming_part(const Surface_Identifier levelSetIdentifier,
                                                           const stk::mesh::Part & nonconformingPart) const
{
  // Look up the original block of the nonconforming part, so this works whether the level set was
  // registered before or after the blocks were decomposed.
  if (!is_nonconformal(nonconformingPart))
    return false;
  const IOPartSet & blocks = myLsUsedByParts.at(levelSetIdentifier);
  return blocks.count(&find_original_part(nonconformingPart)) > 0;
}

} // namespace krino
