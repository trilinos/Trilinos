// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include <Akri_Phase_Support.hpp>

#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_io/IossBridge.hpp>

#include <Akri_AuxMetaData.hpp>
#include <Akri_CDFEM_Support.hpp>
#include <Akri_DiagWriter.hpp>
#include <Akri_LevelSet.hpp>

#include <stk_util/environment/RuntimeDoomed.hpp>
#include <Ioss_ElementTopology.h>
#include <stk_util/util/tokenize.hpp>

namespace 
{
stk::mesh::Part * get_parent_sideset(const stk::mesh::Part & part)
{
  stk::mesh::Part * result = nullptr;
  if(part.primary_entity_rank() == part.mesh_meta_data().side_rank())
  {
    for (auto superset : part.supersets())
    {
      if(stk::io::is_part_io_part(superset))
      {
        ThrowRequireMsg(result == nullptr, "krino::Akri_Phase_Support: Part has more than 1 parent IO part");
        result = superset;
      }
    }
  }
  return result;
}

std::string generate_ioss_sideblock_name(const std::string & base_name, 
  const std::string & block_or_topo_name, const std::string & topo_name)
{
  std::vector<std::string> tokens;
  stk::util::tokenize(base_name, "_", tokens);
  std::string result;
  if(tokens.size() == 2 && tokens.back().find_first_not_of("0123456789") == std::string::npos)
  {
    //Canonically named part
    result = tokens[0] + "_" + block_or_topo_name + "_" + topo_name + "_" + tokens[1];
  }
  else 
  {
    //Non-canonically named part
    result = base_name + "_" + block_or_topo_name + "_" + topo_name;
  }
  return result;
}

std::string build_part_name(const krino::Phase_Support & ps, 
  const stk::mesh::Part & part, const krino::PhaseTag & ls_phase_tag, 
  const std::string & suffix)
{
  stk::mesh::Part * parent_part = get_parent_sideset(part);
  if(parent_part == nullptr)
  {
    return part.name() + suffix;
  }

  std::string parent_part_name = parent_part->name();
  std::string io_part_name = part.name();
  std::set<stk::mesh::PartOrdinal> touching_block_ordinals;
  ps.get_input_blocks_touching_surface(ps.get_input_block_surface_connectivity(), part.mesh_meta_data_ordinal(), touching_block_ordinals);
  

  ThrowRequireMsg(touching_block_ordinals.size() > 0, 
      "krino::Akri_Phase_Support: Side block must be touching at least 1 block");
  stk::topology block_topo = part.mesh_meta_data().get_part(*touching_block_ordinals.begin()).topology();
  for(auto b : touching_block_ordinals)
  {
    ThrowRequireMsg(block_topo == part.mesh_meta_data().get_part(b).topology(), 
      "krino::Akri_Phase_Support: Touching blocks do not all have same topology");
  }
  Ioss::ElementTopology *ioss_side_topo = Ioss::ElementTopology::factory(part.topology().name());
  ThrowRequireMsg(ioss_side_topo != nullptr, 
      "krino::Akri_Phase_Support: IOSS topology factory must return a topology type");
  Ioss::ElementTopology *ioss_block_topo = Ioss::ElementTopology::factory(block_topo.name());
  ThrowRequireMsg(ioss_block_topo != nullptr, 
      "krino::Akri_Phase_Support: IOSS topology factory must return a topology type");

  //Figure out if this side block was split by element block or topology. Ideally need a
  //way to do this without string comparisons
  std::string block_or_topo_name;
  if(io_part_name.find("_"+ioss_block_topo->name()+"_") != std::string::npos)
  {
    //Split by topology
    block_or_topo_name = ioss_block_topo->name();
  }
  else 
  {
    //Split by element block
    ThrowRequireMsg(touching_block_ordinals.size() == 1, 
      "krino::Akri_Phase_Support: Side blocks split by element block should touch exactly 1 block");
    stk::mesh::Part & touching_block = part.mesh_meta_data().get_part(*touching_block_ordinals.begin());
    auto conformal_block = ps.find_conformal_io_part(touching_block, ls_phase_tag);
    if(conformal_block == &touching_block)
    {
      auto decomp_block = part.mesh_meta_data().get_part(touching_block.name() + suffix);
      if(decomp_block != nullptr) conformal_block = decomp_block;
    }
    block_or_topo_name = conformal_block->name();
  }

  std::string base_name = parent_part_name + suffix;
  //Here we can plug in the IOSS function from Greg once it is pushed
  return generate_ioss_sideblock_name(base_name, block_or_topo_name, ioss_side_topo->name());
}
}

namespace krino{

bool Phase_Support::oneLevelSetPerPhase = false;
std::map<std::string,PhaseVec> Phase_Support::the_mesh_phases;
std::map<std::string,PartnamePhasenameMap> Phase_Support::the_mesh_block_phases_by_name;

void Phase_Support::check_phase_parts() const
{
  const LevelSetManager & region_ls = LevelSetManager::get(my_meta);
  bool error = false;
  for (auto&& phase_part : my_phase_parts)
  {
    const PhaseTag & part_phase = (phase_part.is_interface()) ? phase_part.get_touching_phase() : phase_part.get_phase();
    const std::set<LS_SideTag> & phase_ls_sides = part_phase.ls_sides();
    for (std::set<LS_SideTag>::const_iterator side_it = phase_ls_sides.begin(); side_it != phase_ls_sides.end(); ++side_it)
    {
      const std::string phase_ls_name = LevelSet::get_name(side_it->get_ls_identifier());
      if (!region_ls.has_levelSet(phase_ls_name))
      {
        error = true;
        stk::mesh::Part & conformal_part = my_meta.get_part(phase_part.get_conformal_part_ordinal());
        krinolog << "Error: Phase \"" << conformal_part.name() << "\" uses level set \"" << phase_ls_name << "\" but no such level set exists." << stk::diag::dendl;
      }
    }
  }
  ThrowErrorMsgIf(error, "Error: Phases are not defined correctly.");
}

Phase_Support::Phase_Support(stk::mesh::MetaData & meta)
: my_meta(meta),
  my_aux_meta(AuxMetaData::get(meta)),
  nonconformalLsMapsAreFilled_(false) {}

Phase_Support &
Phase_Support::get(stk::mesh::MetaData & meta)
{
  Phase_Support * support = const_cast<Phase_Support *>(meta.get_attribute<Phase_Support>());
  if (nullptr == support)
  {
    support = new Phase_Support(meta);
    ThrowAssert(nullptr != support);
    meta.declare_attribute_with_delete<Phase_Support>(support);
  }
  return *support;
}

Phase_Support &
Phase_Support::get(const stk::mesh::MetaData & meta)
{
  Phase_Support * support = const_cast<Phase_Support *>(meta.get_attribute<Phase_Support>());
  ThrowRequireMsg(nullptr != support, "No Phase_Support found for MetaData.");
  return *support;
}

bool
Phase_Support::exists_and_has_phases_defined(const stk::mesh::MetaData & meta)
{
  Phase_Support * support = const_cast<Phase_Support *>(meta.get_attribute<Phase_Support>());
  return support != nullptr && support->phases_defined();
}

stk::mesh::Selector Phase_Support::get_all_conformal_surfaces_selector() const
{
  std::vector<stk::mesh::Part *> surfaceParts;
  for (auto && part : get_conformal_parts())
    if(is_interface(part))
      surfaceParts.push_back(part);
  return selectUnion(surfaceParts);
}

void
Phase_Support::addPhasePart(stk::mesh::Part & io_part, PhasePartSet & phase_parts, const NamedPhase & ls_phase)
{
  const std::string & phase_name = ls_phase.name();
  const std::string nonconf_suffix = "_nonconformal";

  // The alive phase for death problems just uses the io part as the conformal part
  const std::string suffix = phase_name == "" ? "" : "_" + phase_name;

  std::string phase_part_name = build_part_name(*this, io_part, ls_phase.tag(), suffix);
  std::string nonconf_name = build_part_name(*this, io_part, PhaseTag(), nonconf_suffix);

  stk::mesh::Part & conformal_io_part = my_aux_meta.declare_io_part(phase_part_name, io_part.primary_entity_rank());
  ThrowAssert(phase_name != "" || conformal_io_part.mesh_meta_data_ordinal() == io_part.mesh_meta_data_ordinal());
  std::string topology_name = "INVALID_TOPOLOGY";
  if (stk::topology::INVALID_TOPOLOGY != io_part.topology())
  {
    stk::mesh::set_topology(conformal_io_part, io_part.topology());
    topology_name = conformal_io_part.topology().name();
  }

  if(krinolog.shouldPrint(LOG_PARTS))
  {
    krinolog << "Created conformal IO part " << conformal_io_part.name() << " with topology " << topology_name << stk::diag::dendl;
    krinolog << "  Adding phase part: conformal_part_name=" << conformal_io_part.name()
           << ", nonconformal_part_name=" << nonconf_name
           << ", original_mesh_part_name=" << io_part.name();
    krinolog << " " << ls_phase.tag();
    krinolog << stk::diag::dendl;
  }


  stk::mesh::Part & nonconf_part = my_aux_meta.get_part(nonconf_name);
  phase_parts.insert(PhasePartTag(conformal_io_part.mesh_meta_data_ordinal(), nonconf_part.mesh_meta_data_ordinal(), io_part.mesh_meta_data_ordinal(), ls_phase.tag()));

  all_decomposed_blocks_selector |= conformal_io_part;
  part_is_conformal_map[&conformal_io_part] = true;
  part_to_nonconformal_part_map[&conformal_io_part] = &nonconf_part;
  part_to_phase_map[&conformal_io_part] = ls_phase.tag();
  PhaseTagToPartMap & phase_to_conformal_map = nonconformal_to_phase_conformal_map[&nonconf_part];
  phase_to_conformal_map[ls_phase.tag()] = &conformal_io_part;
}

void
Phase_Support::create_nonconformal_parts(const PartSet & decomposed_ioparts)
{
  const std::string nonconformal_part_suffix = "_nonconformal";
  for(PartSet::const_iterator it = decomposed_ioparts.begin(); it != decomposed_ioparts.end(); ++it)
  {
    const stk::mesh::Part & iopart = my_aux_meta.get_part((*it)->name());

    std::string nonconformal_part_name = build_part_name(*this, iopart, PhaseTag(), nonconformal_part_suffix);

    if(!my_aux_meta.has_part(nonconformal_part_name))
    {
      const bool restartOnlyIOPart = true;
      stk::mesh::Part & nonconformal_part =
          (iopart.topology() == stk::topology::INVALID_TOPOLOGY) ?
          my_aux_meta.declare_io_part(nonconformal_part_name, iopart.primary_entity_rank(), restartOnlyIOPart)    :
          my_aux_meta.declare_io_part_with_topology(nonconformal_part_name, iopart.topology(), restartOnlyIOPart);
      if(krinolog.shouldPrint(LOG_PARTS)) krinolog << "Created nonconformal part " << nonconformal_part_name << stk::diag::dendl;

      all_decomposed_blocks_selector |= iopart;
      all_decomposed_blocks_selector |= nonconformal_part;
      part_to_nonconformal_part_map[&iopart] = &nonconformal_part;
      part_is_nonconformal_map[&nonconformal_part] = true;
    }
  }
}

stk::mesh::PartVector
Phase_Support::get_nonconformal_parts() const
{
  stk::mesh::PartVector result;
  for(const auto & entry : part_is_nonconformal_map)
    if(entry.second)
      result.push_back(const_cast<stk::mesh::Part *>(entry.first));
  return result;
}

stk::mesh::PartVector
Phase_Support::get_nonconformal_parts_of_rank(const stk::mesh::EntityRank rank) const
{
  stk::mesh::PartVector result;
  for(const auto & entry : part_is_nonconformal_map)
    if(entry.second && entry.first->primary_entity_rank() == rank)
      result.push_back(const_cast<stk::mesh::Part *>(entry.first));
  return result;
}

stk::mesh::PartVector
Phase_Support::get_conformal_parts() const
{
  stk::mesh::PartVector result;
  for(const auto & entry : part_is_conformal_map)
  {
    if(entry.second) result.push_back(const_cast<stk::mesh::Part *>(entry.first));
  }
  return result;
}

std::vector<unsigned>
Phase_Support::get_level_set_phases(const PhaseVec & mesh_phases, const LevelSet & levelSet)
{
  // gather up the set of phases that depend on this levelSet
  PhaseTag ls_neg_phase;
  PhaseTag ls_pos_phase;
  ls_neg_phase.add(levelSet.get_identifier(),-1);
  ls_pos_phase.add(levelSet.get_identifier(),+1);
  std::vector<unsigned> ls_phases;
  for(unsigned phase_index=0; phase_index<mesh_phases.size(); ++phase_index)
  {
    const NamedPhase & phase = mesh_phases[phase_index];
    if (has_one_levelset_per_phase() || phase.tag().contain(ls_neg_phase) || phase.tag().contain(ls_pos_phase))
    {
      ls_phases.push_back(phase_index);
    }
  }
  return ls_phases;
}

Phase_Support::PartSet
Phase_Support::get_blocks_and_touching_surfaces(const stk::mesh::MetaData & mesh_meta,
  const stk::mesh::PartVector& input_blocks,
  const Block_Surface_Connectivity & input_block_surface_info)
{
  // For each block that is input_blocks, insert into the set that is returned.
  // Also add surfaces that touch input_blocks to the set.
  PartSet blocks_and_touching_sides;
  for (auto && block_ptr : input_blocks)
  {
    blocks_and_touching_sides.insert(block_ptr);
    std::set<stk::mesh::PartOrdinal> touching_surface_ordinals;
    get_input_surfaces_touching_block(input_block_surface_info, block_ptr->mesh_meta_data_ordinal(), touching_surface_ordinals);
    for (auto && surf_ordinal : touching_surface_ordinals)
    {
      stk::mesh::Part & surf_part = mesh_meta.get_part(surf_ordinal);
      blocks_and_touching_sides.insert(&surf_part);

      // Add all subset parts to the set of parts
      if (!(surf_part.subsets().empty()))
      {
        for (stk::mesh::PartVector::const_iterator subset = surf_part.subsets().begin(); subset != surf_part.subsets().end(); ++subset)
        {
          blocks_and_touching_sides.insert(*subset);
        }
      }
    }
  }
  return blocks_and_touching_sides;
}

void
Phase_Support::create_phase_parts(const PhaseVec& ls_phases, const PartSet& decomposed_ioparts)
{
  for (auto && ls_phase : ls_phases)
  {
    for (auto && io_part : decomposed_ioparts)
    {
      if(get_parent_sideset(*io_part) == nullptr)
        addPhasePart(*io_part, my_phase_parts, ls_phase);
    }
  }
}

void
Phase_Support::subset_and_alias_surface_phase_parts(const PhaseVec& ls_phases,
    const PartSet& decomposed_ioparts)
{
  std::set<stk::mesh::PartOrdinal> touching_block_ordinals;

  for (auto && io_part : decomposed_ioparts)
  {
    if (!(io_part->subsets().empty()))
    {
      for (auto && ls_phase_entry : ls_phases)
      {
        const PhaseTag & ls_phase = ls_phase_entry.tag();
        //FIXME: remove const_cast's
        stk::mesh::Part * conformal_iopart = const_cast<stk::mesh::Part *>(find_conformal_io_part(*io_part, ls_phase));
        ThrowRequire(NULL != conformal_iopart);

        stk::mesh::Part * nonconformal_iopart = const_cast<stk::mesh::Part *>(find_nonconformal_part(*io_part));
        ThrowRequire(NULL != nonconformal_iopart);

        for (stk::mesh::PartVector::const_iterator subset = io_part->subsets().begin();
            subset != io_part->subsets().end(); ++subset)
        {
          stk::mesh::Part * io_part_subset = *subset;
          ThrowRequire(NULL != io_part_subset);

          addPhasePart(*io_part_subset, my_phase_parts, ls_phase_entry);
          stk::mesh::Part * conformal_iopart_subset = const_cast<stk::mesh::Part *>(find_conformal_io_part(*io_part_subset, ls_phase));
          ThrowRequire(NULL != conformal_iopart_subset);

          if(krinolog.shouldPrint(LOG_PARTS)) krinolog << "Adding " << conformal_iopart_subset->name() << " as subset of " << conformal_iopart->name() << stk::diag::dendl;
          my_meta.declare_part_subset(*conformal_iopart, *conformal_iopart_subset);

          stk::mesh::Part * nonconformal_iopart_subset = const_cast<stk::mesh::Part *>(find_nonconformal_part(*io_part_subset));
          ThrowRequire(NULL != nonconformal_iopart_subset);

          if(krinolog.shouldPrint(LOG_PARTS)) krinolog << "Adding " << nonconformal_iopart_subset->name() << " as subset of " << nonconformal_iopart->name() << stk::diag::dendl;
          my_meta.declare_part_subset(*nonconformal_iopart, *nonconformal_iopart_subset);

          get_input_blocks_touching_surface(my_input_block_surface_connectivity, io_part_subset->mesh_meta_data_ordinal(), touching_block_ordinals);
          for (auto && touching_block_ordinal : touching_block_ordinals)
          {
            const std::string conformal_part_alias = conformal_iopart->name() + "_" + my_meta.get_part(touching_block_ordinal).name();
            if(krinolog.shouldPrint(LOG_PARTS)) krinolog << "Adding alias " << conformal_part_alias << " for conformal part " << conformal_iopart_subset->name() << stk::diag::dendl;
            my_aux_meta.define_part_alias(*conformal_iopart_subset, conformal_part_alias);
          }
        }
      }
    }
  }
}

void
Phase_Support::build_decomposed_block_surface_connectivity()
{
  for (auto && part : my_meta.get_mesh_parts())
  {
    if (part->primary_entity_rank() != my_meta.side_rank()) continue;
    const PhasePartTag * phase_part = find_conformal_phase_part(*part);
    if (nullptr == phase_part) continue;
    stk::mesh::Part & orig_part = my_meta.get_part(phase_part->get_original_part_ordinal());
    if (orig_part == my_meta.universal_part()) continue;

    if (phase_part->is_interface())
    {
      const stk::mesh::Part * conformal_touching_block = find_conformal_io_part(orig_part, phase_part->get_touching_phase());
      ThrowRequire(conformal_touching_block);
      if(krinolog.shouldPrint(LOG_PARTS)) krinolog << "Surface " << part->name() << " touches block " << conformal_touching_block->name() << "\n";
      std::vector<const stk::mesh::Part*> touching_blocks = my_meta.get_blocks_touching_surface(part);
      if (std::find(touching_blocks.begin(), touching_blocks.end(), conformal_touching_block) == touching_blocks.end())
      {
        touching_blocks.push_back(conformal_touching_block);
      }
      my_meta.set_surface_to_block_mapping(part, touching_blocks);
    }
    else
    {
      std::set<stk::mesh::PartOrdinal> touching_block_ordinals;
      get_input_blocks_touching_surface(my_input_block_surface_connectivity, orig_part.mesh_meta_data_ordinal(), touching_block_ordinals);

      for (auto && touching_block_ordinal : touching_block_ordinals)
      {
        stk::mesh::Part & touching_block = my_meta.get_part(touching_block_ordinal);
        const stk::mesh::Part * conformal_touching_block = find_conformal_io_part(touching_block, phase_part->get_phase());
        ThrowRequire(conformal_touching_block);
        if(krinolog.shouldPrint(LOG_PARTS)) krinolog << "Surface " << part->name() << " touches block " << conformal_touching_block->name() << "\n";
        std::vector<const stk::mesh::Part*> touching_blocks = my_meta.get_blocks_touching_surface(part);
        if (std::find(touching_blocks.begin(), touching_blocks.end(), conformal_touching_block) == touching_blocks.end())
        {
          touching_blocks.push_back(conformal_touching_block);
        }
        my_meta.set_surface_to_block_mapping(part, touching_blocks);
      }
    }
  }
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
    ThrowRequire((int) interface_phases.size() == num_interfaces);

    for (int i = 0; i < num_interfaces; ++i)
    {
      const int name_compare = ls_phases[interface_phases[i].first].name().compare(ls_phases[interface_phases[i].second].name());
      ThrowRequire(name_compare != 0);
      // form phase intersection
      const auto & phase0 = (name_compare < 0) ? ls_phases[interface_phases[i].first] : ls_phases[interface_phases[i].second];
      const auto & phase1 = (name_compare < 0) ? ls_phases[interface_phases[i].second] : ls_phases[interface_phases[i].first];
      std::string superset_interface_phase_name = phase0.name() + "_" + phase1.name();
      const std::string superset_phase_part_name = interface_name_gen.interface_superset_name(superset_interface_phase_name);

      // assumes all sides have same topology
      stk::mesh::Part &superset_phase_part = my_aux_meta.declare_io_part(superset_phase_part_name, my_meta.side_rank());

      if(krinolog.shouldPrint(LOG_PARTS)) krinolog << "Adding interfacial phase part: conformal_part_name=" << superset_phase_part.name() << "\n";

      for (auto && phase_pair : {std::make_pair(phase0,phase1), std::make_pair(phase1,phase0)})
      {
        const auto & touching_phase = phase_pair.first.tag();
        const auto & opposite_phase = phase_pair.second.tag();

        my_phase_parts.insert(PhasePartTag(
                superset_phase_part.mesh_meta_data_ordinal(), my_meta.universal_part().mesh_meta_data_ordinal(),
                my_meta.universal_part().mesh_meta_data_ordinal(), touching_phase,
                opposite_phase));

        for (auto * io_part : decomposed_ioparts)
        {
          ThrowRequire(NULL != io_part);
          // only handle blocks for now -> creating interfacial surface "phases"
          if (io_part->primary_entity_rank() == stk::topology::ELEMENT_RANK)
          {
            std::string subset_interface_phase_name = phase_pair.first.name() + "_" + phase_pair.second.name();
            const std::string subset_phase_part_name = interface_name_gen.interface_subset_name(io_part->name(), subset_interface_phase_name);
            const bool restartOnlyIOPart = true; // Subset part is only output for restart
            stk::mesh::Part &subset_phase_part = my_aux_meta.declare_io_part_with_topology(subset_phase_part_name, io_part->topology().side_topology(), restartOnlyIOPart);
            my_meta.declare_part_subset(superset_phase_part, subset_phase_part);

            const stk::mesh::Part * const nonconf_part = find_nonconformal_part(*io_part);
            ThrowRequire(NULL != nonconf_part);

            if(krinolog.shouldPrint(LOG_PARTS))
            {
              krinolog << "  Adding subset interfacial phase part: conformal_part_name=" << subset_phase_part.name()
                  << ", nonconformal_part_name=" << nonconf_part->name();
              krinolog << ", touching_phase=" << touching_phase
                       << ", opposite_phase=" << opposite_phase;
              krinolog << "\n";
            }

            my_phase_parts.insert(PhasePartTag(
                subset_phase_part.mesh_meta_data_ordinal(), nonconf_part->mesh_meta_data_ordinal(),
                io_part->mesh_meta_data_ordinal(), touching_phase,
                opposite_phase));

            all_decomposed_blocks_selector |= subset_phase_part;
            part_is_conformal_map[&subset_phase_part] = true;
            part_to_nonconformal_part_map[&subset_phase_part] = nonconf_part;
            auto & phase_to_conformal_map = nonconformal_to_phase_conformal_map[nonconf_part];
            auto & touching_phase_part = phase_to_conformal_map[touching_phase];
            auto & opposite_phase_part = phase_to_conformal_map[opposite_phase];
            auto key = std::make_pair(touching_phase_part, opposite_phase_part);
            volume_to_interface_parts_map[key] = &subset_phase_part;
          }
        }
      }
    }
  }
}

void
Phase_Support::decompose_blocks(std::vector<std::tuple<stk::mesh::PartVector, 
  std::shared_ptr<Interface_Name_Generator>, PhaseVec>> ls_sets)
{
  stk::topology topo;
  for(auto && ls_set : ls_sets)
  {
    stk::mesh::PartVector & blocks_to_decompose = std::get<0>(ls_set);
    PhaseVec & ls_phases = std::get<2>(ls_set);
    if(ls_phases.empty()) continue;
    for(auto && decompose_block : blocks_to_decompose)
    {
      topo = decompose_block->topology();
      if(topo != stk::topology::TET_4 && topo != stk::topology::TET_10
        && topo != stk::topology::TRI_3_2D && topo != stk::topology::TRI_6_2D)
      {
        stk::RuntimeDoomedAdHoc() << "Currently only Tetrahedron 4 or 10 (in 3D) and Triangle 3 or 6 (in 2D) element types are supported when using CDFEM.\n";
        return;
      }
    }
  }

  std::vector<PartSet> part_set_vec(ls_sets.size());
  for(unsigned i=0; i<ls_sets.size(); i++)
  {
    auto & ls_set = ls_sets[i];
    stk::mesh::PartVector & blocks_to_decompose = std::get<0>(ls_set);
    PhaseVec & ls_phases = std::get<2>(ls_set);
    if(std::get<2>(ls_set).empty()) continue;
    part_set_vec[i] = get_blocks_and_touching_surfaces(my_meta, blocks_to_decompose, my_input_block_surface_connectivity);

    create_nonconformal_parts(part_set_vec[i]);

    create_phase_parts(ls_phases, part_set_vec[i]);
  }

  for(unsigned i=0; i<ls_sets.size(); i++)
  {
    auto & ls_set = ls_sets[i];
    std::shared_ptr<Interface_Name_Generator> name_gen = std::get<1>(ls_set);
    PhaseVec & ls_phases = std::get<2>(ls_set);
    if(ls_phases.empty()) continue;

    subset_and_alias_surface_phase_parts(ls_phases, part_set_vec[i]);

    create_interface_phase_parts(ls_phases, part_set_vec[i], *name_gen);
  }
}

//--------------------------------------------------------------------------------

void
Phase_Support::determine_block_phases(const std::set<std::string> & FEmodel_block_names, const krino::PhaseVec & FEmodel_phases)
{
  if (FEmodel_phases.empty()) return;

  ThrowAssertMsg(my_mesh_block_phases.empty(), "determine_block_phases should only be called once per mesh");

  for (auto && part : my_meta.get_parts())
  {
    if (part->primary_entity_rank() == stk::topology::ELEMENT_RANK)
    {
      auto & block_phases = my_mesh_block_phases[part->mesh_meta_data_ordinal()];
      for(unsigned phase_index=0; phase_index<FEmodel_phases.size(); ++phase_index)
      {
        const std::string & phase_name = FEmodel_phases[phase_index].name();
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

void
Phase_Support::determine_block_phases(const PhaseVec & FEmodel_phases, const PartnamePhasenameMap & FEmodel_block_phase_names)
{
  if (FEmodel_phases.empty()) return;

  ThrowAssertMsg(my_mesh_block_phases.empty(), "determine_block_phases should only be called once per mesh");

  if (FEmodel_block_phase_names.empty())
  {
    stk::RuntimeDoomedAdHoc() << "Failure in determine_block_phases. Were the subdomains specified for the decomposed blocks in the finite element model specification?\n";
    return;
  }

  std::map<std::string, unsigned> phase_name_to_index_map;
  for(unsigned phase_index=0; phase_index<FEmodel_phases.size(); ++phase_index)
  {
    const std::string & phase_name = FEmodel_phases[phase_index].name();
    phase_name_to_index_map[phase_name] = phase_index;
  }

  for (auto && map_entry : FEmodel_block_phase_names)
  {
    if (my_aux_meta.has_part(map_entry.first))
    {
      stk::mesh::Part & block = my_aux_meta.get_part(map_entry.first);
      auto & block_phases = my_mesh_block_phases[block.mesh_meta_data_ordinal()];
      const auto & block_phase_names = map_entry.second;
      for(auto && block_phase_name : block_phase_names)
      {
        auto map_iter = phase_name_to_index_map.find(block_phase_name);
        if (map_iter == phase_name_to_index_map.end())
        {
          stk::RuntimeDoomedAdHoc() << "Failed to find phase " << block_phase_name << ", specified for block " << block.name() << " in the list of all phases.\n";
        }
        else
        {
          block_phases.insert(map_iter->second);
        }
      }
    }
    else
    {
      stk::RuntimeDoomedAdHoc() << "The block " << map_entry.first << " is listed in the finite element model specification, but is not in the input mesh.\n";
    }
  }
}

//--------------------------------------------------------------------------------

std::vector<stk::mesh::Part*>
Phase_Support::get_blocks_decomposed_by_levelset(const std::vector<unsigned> ls_phases) const
{
  std::vector<stk::mesh::Part*> blocks_decomposed_by_ls;
  for (auto && map_it : my_mesh_block_phases)
  {
    stk::mesh::Part & FEmodel_block = my_meta.get_part(map_it.first);
    const auto & block_phases = map_it.second;
    for(auto ls_phase : ls_phases)
    {
      if (block_phases.find(ls_phase) != block_phases.end())
      {
        blocks_decomposed_by_ls.push_back(&FEmodel_block);
        break;
      }
    }
  }

  return blocks_decomposed_by_ls;
}

//--------------------------------------------------------------------------------

void
Phase_Support::setup_phases(const krino::PhaseVec & FEmodel_phases)
{
  std::vector<std::tuple<stk::mesh::PartVector, std::shared_ptr<Interface_Name_Generator>, PhaseVec>>
    ls_sets;
  auto interface_name_gen = std::shared_ptr<Interface_Name_Generator>(new LS_Name_Generator());
  for (auto && map_it : my_mesh_block_phases)
  {
    stk::mesh::Part & FEmodel_block = my_meta.get_part(map_it.first);
    const auto & block_phase_ordinals = map_it.second;

    PhaseVec block_phases;
    for(auto && block_phase_ordinal : block_phase_ordinals)
    {
      block_phases.push_back(FEmodel_phases[block_phase_ordinal]);
    }
    ls_sets.push_back(std::make_tuple(stk::mesh::PartVector{&FEmodel_block},interface_name_gen,block_phases));
  }
  decompose_blocks(ls_sets);
  build_decomposed_block_surface_connectivity();
}

//--------------------------------------------------------------------------------

void
Phase_Support::get_iopart_roots(const stk::mesh::Part & iopart, std::vector<const stk::mesh::Part *> & subsets)
{ /* %TRACE% */  /* %TRACE% */
  if ( iopart.subsets().empty() )
  {
    subsets.push_back(&iopart);
  }
  else
  {
    for (stk::mesh::PartVector::const_iterator subset = iopart.subsets().begin(); subset != iopart.subsets().end() ; ++subset )
    {
      ThrowRequire(stk::io::is_part_io_part(**subset));
      get_iopart_roots(**subset, subsets);
    }
  }
}
//--------------------------------------------------------------------------------
void
Phase_Support::get_blocks_touching_surface(const std::string & surface_name, std::vector<std::string> & block_names)
{ /* %TRACE% */  /* %TRACE% */

  std::set<stk::mesh::PartOrdinal> block_ordinal_set;
  for (auto&& block_ordinal : block_ordinal_set)
  {
    block_names.push_back(my_meta.get_part(block_ordinal).name());
  }
}
//--------------------------------------------------------------------------------
void
Phase_Support::get_input_surfaces_touching_block(const Block_Surface_Connectivity & input_block_surface_connectivity,
    const stk::mesh::PartOrdinal block_ordinal, std::set<stk::mesh::PartOrdinal> & surface_ordinals)
{
  input_block_surface_connectivity.get_surfaces_touching_block(block_ordinal, surface_ordinals);
}

//--------------------------------------------------------------------------------
void
Phase_Support::get_input_blocks_touching_surface(const Block_Surface_Connectivity & input_block_surface_connectivity,
    const stk::mesh::PartOrdinal surfaceOrdinal, std::set<stk::mesh::PartOrdinal> & blockOrdinals) const
{
  input_block_surface_connectivity.get_blocks_touching_surface(surfaceOrdinal, blockOrdinals);
}
//--------------------------------------------------------------------------------
const stk::mesh::Part *
Phase_Support::find_conformal_io_part(const stk::mesh::Part & io_part, const PhaseTag & phase) const
{
  const stk::mesh::Part * const nonconf_part = find_nonconformal_part(io_part);
  auto entry = nonconformal_to_phase_conformal_map.find(nonconf_part);
  if(entry == nonconformal_to_phase_conformal_map.end()) return &io_part;
  PhaseTagToPartMap & tag_to_conformal_map = (*entry).second;

  const auto find_part = tag_to_conformal_map.find(phase);
  if(find_part != tag_to_conformal_map.end())
  {
    auto & tag_part_pair = *find_part;
    return tag_part_pair.second;
  }

  // This search is to handle the case where the phase passed in contains additional LS side tags that the
  // conformal part phase tag doesn't care about. For example if the conformal phase is defined as "where LS1 is negative"
  // it should match a phase tag that has (0,1) (1,-1) even though the tags are not ==
  // TODO: If we knew about all the LS's in the problem when decompose_blocks was called we could populate all the necessary entries
  // then and wouldn't need this search.
  for(auto && tag_part_pair : tag_to_conformal_map)
  {
    if(phase.contain(tag_part_pair.first))
    {
      tag_to_conformal_map[phase] = tag_part_pair.second;
      return tag_part_pair.second;
    }
  }

  return &io_part;
}

//--------------------------------------------------------------------------------

bool
Phase_Support::is_conformal(const stk::mesh::Part * io_part) const
{
  PartToBoolMap::const_iterator entry = part_is_conformal_map.find(io_part);
  if (entry != part_is_conformal_map.end()) return (entry->second);
  return false;
}

//--------------------------------------------------------------------------------

bool
Phase_Support::is_nonconformal(const stk::mesh::Part * io_part) const
{
  PartToBoolMap::const_iterator entry = part_is_nonconformal_map.find(io_part);
  if (entry != part_is_nonconformal_map.end()) return (entry->second);
  return false;
}

//--------------------------------------------------------------------------------

bool
Phase_Support::is_interface(const stk::mesh::Part * io_part) const
{
  ThrowAssert(io_part);
  auto phase_part = find_conformal_phase_part(*io_part);
  return phase_part && phase_part->is_interface();
}

//--------------------------------------------------------------------------------

const stk::mesh::Part *
Phase_Support::find_nonconformal_part(const stk::mesh::Part & io_part) const
{
  PartToPartMap::const_iterator entry = part_to_nonconformal_part_map.find(&io_part);
  if (entry != part_to_nonconformal_part_map.end()) return (entry->second);
  return &io_part;
}

//--------------------------------------------------------------------------------

const stk::mesh::Part *
Phase_Support::find_interface_part(const stk::mesh::Part & vol0, const stk::mesh::Part & vol1) const
{
  auto key = std::make_pair(&vol0, &vol1);
  auto find_it = volume_to_interface_parts_map.find(key);
  if(find_it != volume_to_interface_parts_map.end()) return find_it->second;

  return nullptr;
}

//--------------------------------------------------------------------------------

const PhaseTag &
Phase_Support::get_iopart_phase(const stk::mesh::Part & io_part) const
{
  PartToPhaseTagMap::const_iterator entry = part_to_phase_map.find(&io_part);
  if (entry != part_to_phase_map.end()) return (entry->second);

  static const PhaseTag empty_phase;
  return empty_phase;
}

//--------------------------------------------------------------------------------

void Phase_Support::register_blocks_for_level_set(LevelSet * levelSet,
                                                  const std::vector<stk::mesh::Part *> & blocks_decomposed_by_ls)
{
  nonconformalLsMapsAreFilled_ = false;

  // Loop over block names
  for (auto && block_ptr : blocks_decomposed_by_ls)
  {
    // Add direct relationship between this block and level set
    lsUsedByParts_[levelSet].insert(block_ptr);

    // Now get surfaces touching this block
    std::set<stk::mesh::PartOrdinal> surfaceOrdinals;
    get_input_surfaces_touching_block(my_input_block_surface_connectivity, block_ptr->mesh_meta_data_ordinal(), surfaceOrdinals);
    for (auto && surfaceOrdinal : surfaceOrdinals)
    {
      // For each surface, add IO Part/Level Set pairing to maps
      stk::mesh::Part & surfaceIOPart = my_meta.get_part(surfaceOrdinal);
      lsUsedByParts_[levelSet].insert(&surfaceIOPart);
    }
  }
}


//--------------------------------------------------------------------------------
bool Phase_Support::level_set_is_used_by_nonconformal_part(const LevelSet * levelSet,
                                                           const stk::mesh::Part * const ioPart) const
{
  if (!nonconformalLsMapsAreFilled_)
    fill_nonconformal_level_set_maps();

  IOPartSet & ioPartSet = lsUsedByNonconformalParts_[levelSet];
  bool contains = (ioPartSet.find(ioPart) != ioPartSet.end());
  return contains;
}

//--------------------------------------------------------------------------------
void Phase_Support::fill_nonconformal_level_set_maps() const
{
  if (nonconformalLsMapsAreFilled_ == true) return;
  nonconformalLsMapsAreFilled_ = true;

  // In this method we duplicate the level-set-to-part maps, but use
  // the nonconformal versions of each IOPart instead

  // lsUsedByNonconformalParts_
  lsUsedByNonconformalParts_.clear();
  std::map<LevelSet *, IOPartSet>::const_iterator iterLSMap;
  for (iterLSMap = lsUsedByParts_.begin(); iterLSMap != lsUsedByParts_.end(); ++iterLSMap)
  {
    const LevelSet * levelSet = iterLSMap->first;
    const IOPartSet & ioPartSet = iterLSMap->second;
    IOPartSet::const_iterator iterPart;
    for (iterPart = ioPartSet.begin(); iterPart != ioPartSet.end(); ++iterPart)
    {
      const stk::mesh::Part * ioPart = *iterPart;
      const stk::mesh::Part * nonconformalIOPart = find_nonconformal_part(*ioPart);
      lsUsedByNonconformalParts_[levelSet].insert(nonconformalIOPart);
    }
  }
}

const PhasePartTag * Phase_Support::find_conformal_phase_part(const stk::mesh::Part & conformal_part) const
{
  const PhasePartTag * result = NULL;

  const unsigned conformal_part_ordinal = conformal_part.mesh_meta_data_ordinal();

  // Take advantage of the fact that PhasePartTags are uniquely identified by the conformal_part_ordinal
  // and therefore the set of phase parts is sorted only based on that to use the binary search
  // my_phase_parts.find().
  PhaseTag dummy_tag;
  dummy_tag.add(LevelSet_Identifier(0), 1);
  PhasePartTag search_tag(conformal_part_ordinal, 0, 0, dummy_tag);

  auto it =  my_phase_parts.find(search_tag);
  if(it != my_phase_parts.end()) result = &(*it);

  return result;
}

//--------------------------------------------------------------------------------
CDFEM_Irreversible_Phase_Support &
CDFEM_Irreversible_Phase_Support::get(stk::mesh::MetaData & meta)
{
  CDFEM_Irreversible_Phase_Support * support = const_cast<CDFEM_Irreversible_Phase_Support *>(meta.get_attribute<CDFEM_Irreversible_Phase_Support>());
  if (NULL == support)
  {
    support = new CDFEM_Irreversible_Phase_Support();
    meta.declare_attribute_with_delete<CDFEM_Irreversible_Phase_Support>(support);
  }
  return *support;
}
//--------------------------------------------------------------------------
CDFEM_Inequality_Spec * CDFEM_Irreversible_Phase_Support::add_death_spec(const std::string & death_name, bool is_death)
{
  // CDFEM death needs the decomposition to occur at the end of the time step, irreversible phase change
  // needs it to occur at the start of the time step. For now throw if the user tries to do both in the same problem.
  if(is_death)
  {
    has_death = true;
  }
  else
  {
    has_irreversible_phase_change = true;
  }
  ThrowInvalidArgMsgIf(has_death && has_irreversible_phase_change,
    "Cannot have both CDFEM death and CDFEM irreversible phase change in the same problem.");

  for (auto&& death_spec : my_death_specs)
  {
    ThrowInvalidArgMsgIf(death_spec.name() == death_name,
      "Only one CDFEM death specification with the same name is allowed per mesh (" << death_name << "). ");
  }
  my_death_specs.push_back(CDFEM_Inequality_Spec(death_name));
  return &my_death_specs.back();
}

//--------------------------------------------------------------------------

CDFEM_Inequality_Spec::CDFEM_Inequality_Spec(const std::string & name_)
   : my_name(name_),
     my_criterion_compare_type(INEQUALITY_CRITERION_TYPE_UNDEFINED),
     my_ls(NULL)
{}

CDFEM_Inequality_Spec::InequalityCriterionType
CDFEM_Inequality_Spec::int_to_inequality_criterion_type(const int inequality_criterion)
{
  /* %TRACE[ON]% */ Trace trace__("krino::CDFEM_Death_Spec::int_to_death_criterion_compare_type(const int death_criterion_compare)"); /* %TRACE% */

  switch (inequality_criterion)
  {
  case int(INEQUALITY_CRITERION_TYPE_UNDEFINED):
    return INEQUALITY_CRITERION_TYPE_UNDEFINED ;
  case int(LESS_THAN                             ):
    return LESS_THAN                              ;
  case int(GREATER_THAN                          ):
    return GREATER_THAN                           ;
  default:
    ThrowRuntimeError("Bad integer passed.  Could not convert passed in integer "
                      << "death_criterion_compare = " << inequality_criterion
                      << " into a CDFEM_Death_Spec::DeathCriterionCompareType enum.  "
                      << "The code and the Sierra XML database may be out of sync." << std::endl
                      << StackTrace);
  }
}

void CDFEM_Inequality_Spec::add_element_volume_name (const std::string & a_element_volume_name)
{
  /* %TRACE[NONE]% */  /* %TRACE% */

  // Note: It is ok to add multiple volume names to vector, thus no return status given
  my_element_volume_names.push_back(a_element_volume_name);
}

bool CDFEM_Inequality_Spec::set_threshold_variable_name (const std::string & a_threshold_variable_name)
{
  /* %TRACE[ON]% */ Trace trace__("krino::CDFEM_Death_Spec::set_death_variable_name(const std::string & a_death_variable_name)"); /* %TRACE% */

  // Set return status to false if value has already been set
  bool return_status = (0 == my_threshold_variable_name.length());

  if (return_status)
    my_threshold_variable_name = a_threshold_variable_name;

  return return_status;
}

bool CDFEM_Inequality_Spec::set_threshold_value (const double a_threshold_value)
{
  /* %TRACE[ON]% */ Trace trace__("krino::CDFEM_Death_Spec::set_threshold_value(const double a_threshold_value)"); /* %TRACE% */

  // Set return status to false if value has already been set
  bool return_status = (std::numeric_limits<double>::max() != a_threshold_value);

  if (return_status)
    my_threshold_value = a_threshold_value;

  return return_status;
}

bool CDFEM_Inequality_Spec::set_criterion_compare_type
  (const CDFEM_Inequality_Spec::InequalityCriterionType & a_criterion_compare_type)
{
  /* %TRACE[ON]% */ Trace trace__("krino::CDFEM_Death_Spec::set_criterion_compare_type(const CDFEM_Death_Spec::DeathCriterionCompareType & a_criterion_compare_type)"); /* %TRACE% */

  // Set return status to false if value has already been set
  bool return_status = (CDFEM_Inequality_Spec::INEQUALITY_CRITERION_TYPE_UNDEFINED == my_criterion_compare_type);

  if (return_status)
    my_criterion_compare_type = a_criterion_compare_type;

  return return_status;
}

void CDFEM_Inequality_Spec::sanity_check(const std::vector<std::string> mesh_elem_blocks) const
{
  /* %TRACE[ON]% */ Trace trace__("krino::CDFEM_Death_Spec::sanity_check()"); /* %TRACE% */

  ThrowErrorMsgIf(CDFEM_Inequality_Spec::INEQUALITY_CRITERION_TYPE_UNDEFINED == my_criterion_compare_type,
   "CDFEM Death is not properly setup.  Was the death criterion specified?");

  const std::vector<std::string> & volume_names = get_element_volume_names();

  for (auto && block_name : volume_names)
  {
    ThrowErrorMsgIf(
        std::find(mesh_elem_blocks.begin(), mesh_elem_blocks.end(), block_name) == mesh_elem_blocks.end(),
        "Could not find an element volume named '"
        << block_name << "' that was specified in the "
        << "cdfem death command block for the cdfem death named '"
        << name() << "'.");
  }
}

void CDFEM_Inequality_Spec::create_levelset(stk::mesh::MetaData & meta, stk::diag::Timer & parent_timer)
{
  /* %TRACE[ON]% */ Trace trace__("krino::CDFEM_Death_Spec::create_levelset()"); /* %TRACE% */

  ThrowInvalidArgMsgIf(
      LevelSetManager::get(meta).has_levelSet(my_name),
      "Region already has a LevelSet named " << my_name);

  my_ls = &LevelSet::build(meta, my_name, parent_timer);
  my_ls->set_isovar(get_threshold_variable_name(), get_threshold_value());

  int deactivated_phase_sign = ( get_criterion_compare_type() == CDFEM_Inequality_Spec::LESS_THAN ) ? -1 : 1;
  my_deactivated_phase.add(my_ls->get_identifier(), deactivated_phase_sign);
  my_active_phase.add(my_ls->get_identifier(), -1 * deactivated_phase_sign);
}

} // namespace krino
