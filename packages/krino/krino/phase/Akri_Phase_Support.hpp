// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#ifndef Akri_Phase_Support_h
#define Akri_Phase_Support_h
//
#include <stk_mesh/base/Part.hpp>
#include <stk_mesh/base/Selector.hpp>
#include <stk_mesh/base/Types.hpp>
#include <stk_util/util/ReportHandler.hpp>

#include <Akri_Interface_Name_Generator.hpp>
#include <Akri_IO_Helpers.hpp>
#include <Akri_PhaseInput.hpp>
#include <Akri_PhasePartInfo.hpp>
#include <Akri_PhaseTag.hpp>
#include <Akri_Surface_Identifier.hpp>

#include <map>
#include <memory>
#include <set>
#include <string>
#include <tuple>
#include <vector>

namespace stk { namespace mesh { class MetaData; } }

namespace krino {

class AuxMetaData;
struct LS_Field; // TEMPORARY: see get_levelset_fields() below.

typedef std::set< const stk::mesh::Part * > IOPartSet;
typedef std::map<stk::mesh::PartOrdinal,std::set<unsigned>> PartPhaseMap;

class DecompositionPackage
{
public:
  unsigned size() const { return myDecompositions.size(); }
  const stk::mesh::PartVector & get_blocks_to_decompose(unsigned i) const { return myDecompositions[i].blocks; }
  const PhaseVec & get_phases(unsigned i) const { return myDecompositions[i].phases; }
  const Interface_Name_Generator & get_interface_name_generator(unsigned i) const { return *myDecompositions[i].interfaceNameGenerator; }
  // User-facing description of the decomposition, for error messages.
  const std::string & get_description(unsigned i) const { return myDecompositions[i].description; }

  void add_levelset_decomposition(const stk::mesh::PartVector & blocksToDecompose, const PhaseVec & namedPhases)
  {
    add_decomposition_if_not_empty(blocksToDecompose, namedPhases, std::make_shared<LS_Name_Generator>(), "level set decomposition");
  }
  void add_death_decomposition(const stk::mesh::PartVector & blocksToDecompose, const PhaseVec & namedPhases, const std::string & deathName)
  {
    add_decomposition_if_not_empty(blocksToDecompose, namedPhases, std::make_shared<Death_Name_Generator>(deathName), "CDFEM death " + deathName);
  }

private:
  struct Decomposition
  {
    stk::mesh::PartVector blocks;
    PhaseVec phases;
    std::shared_ptr<Interface_Name_Generator> interfaceNameGenerator;
    std::string description;
  };

  void add_decomposition_if_not_empty(const stk::mesh::PartVector & blocksToDecompose, const PhaseVec & namedPhases, const std::shared_ptr<Interface_Name_Generator> & nameGenerator, const std::string & description)
  {
    if (!blocksToDecompose.empty() && !namedPhases.empty())
      myDecompositions.push_back(Decomposition{blocksToDecompose, namedPhases, nameGenerator, description});
  }
  std::vector<Decomposition> myDecompositions;
};

class Phase_Support {
public:
  Phase_Support (const Phase_Support&) = delete;
  Phase_Support& operator= (const Phase_Support&) = delete;

  typedef std::set<unsigned> PartSet; // part ordinals; std::set orders by ordinal, matching the old stk::mesh::PartLess

  static bool exists_and_has_phases_defined(const stk::mesh::MetaData & meta);
  static Phase_Support & get(const stk::mesh::MetaData & meta);
  static Phase_Support & get(stk::mesh::MetaData & meta); // for unit test usage where the Phase_Support is created with the MetaData

  std::vector<unsigned> get_negative_levelset_block_ordinals(const Surface_Identifier levelSetIdentifier) const;
  std::vector<unsigned> get_negative_levelset_interface_ordinals(const Surface_Identifier levelSetIdentifier) const;
  stk::mesh::Selector get_negative_levelset_interface_selector(const Surface_Identifier levelSetIdentifier) const;
  stk::mesh::Selector get_negative_levelset_block_selector(const Surface_Identifier levelSetIdentifier) const;

  std::vector<unsigned> get_levelset_decomposed_block_ordinals(const Surface_Identifier levelSetIdentifier) const;
  stk::mesh::Selector get_levelset_decomposed_blocks_selector(const Surface_Identifier levelSetIdentifier) const;

  bool phases_defined() const { return !myPhasePartInfo.empty(); }
  bool is_cdfem_use_case() const;
  void force_cdfem_use_case_for_minimal_unit_tests();

  // Reports (RuntimeDoomedAdHoc) each sideset touching blocks that are decomposed with different
  // criteria (phase tags) for a phase with the same name, since the sideset's phase parts are shared.
  // Returns true if there are none.  Call at parse time, before decompose_blocks, which throws on
  // such sidesets.
  bool validate_surface_phase_connectivity(const DecompositionPackage & decomps) const;
  void decompose_blocks(const DecompositionPackage & decomps);
  std::vector<stk::mesh::Part*> get_blocks_decomposed_by_surface(const std::vector<unsigned> & surfacePhases) const;
  void setup_phases();
  void set_input_block_surface_connectivity(const Block_Surface_Connectivity & input_block_surface_info) { my_input_block_surface_connectivity = input_block_surface_info; }
  const Block_Surface_Connectivity & get_input_block_surface_connectivity() const {return my_input_block_surface_connectivity;}
  void build_decomposed_block_surface_connectivity();

  std::map<int,int> build_phase_conversion_map(const PhaseTag & fromPhase, const PhaseTag & toPhase) const;
  stk::mesh::Selector get_selector_of_blocks_with_phase(const PhaseTag & phase) const;

  bool is_nonconformal(const stk::mesh::Part & io_part) const;
  bool is_conformal(const stk::mesh::Part & io_part) const;
  bool is_interface(const stk::mesh::Part & io_part) const;
  stk::mesh::Part * find_conforming_part_if_decomposed_by_phase(const stk::mesh::Part & part, const PhaseTag & phase) const;
  stk::mesh::Part & get_conforming_part(const stk::mesh::Part & part, const PhaseTag & phase) const;
  stk::mesh::Part & get_nonconforming_part(const stk::mesh::Part & part) const;
  stk::mesh::Part & find_original_part(const stk::mesh::Part & part) const;
  const stk::mesh::Part * find_interface_part(const stk::mesh::Part & vol0, const stk::mesh::Part & vol1) const;
  const PhaseTag & get_iopart_phase(const stk::mesh::Part & io_part) const;

  void set_all_decomposed_blocks_selector();
  const stk::mesh::Selector & get_all_decomposed_blocks_selector() const { return all_decomposed_blocks_selector; }
  bool is_decomposed(const stk::mesh::Part & part) const;

  void register_blocks_for_level_set(const Surface_Identifier levelSetIdentifier,
      const std::vector<stk::mesh::Part *> & blocks_decomposed_by_ls);
  void register_blocks_for_surface(const Surface_Identifier & surfID);
  stk::mesh::Selector get_all_interface_surfaces_selector() const;

  // True if nonconformingPart is the nonconforming part of a block registered for the level set.
  bool level_set_is_used_by_nonconforming_part(const Surface_Identifier levelSetIdentifier, const stk::mesh::Part & nonconformingPart) const;

  bool has_one_levelset_per_phase() const { return oneLevelSetPerPhase; }
  void set_one_levelset_per_phase(const bool val) { oneLevelSetPerPhase = val; }

  // Sets the phases parsed for this mesh's finite element model.  Must be called before any phases are set.
  void set_phase_input(const PhaseInput & phaseInput);
  const PhaseVec & get_mesh_phases() const { return myMeshPhases; }
  const PartnamePhasenameMap & get_block_phases_by_name() const { return myMeshBlockPhasesByName; }

  stk::mesh::PartVector get_nonconformal_parts() const;
  stk::mesh::PartVector get_nonconformal_parts_of_rank(const stk::mesh::EntityRank rank) const;
  stk::mesh::PartVector get_conformal_parts(const bool includeInterfaceParts) const;
  stk::mesh::PartVector get_conformal_parts_of_rank(const stk::mesh::EntityRank rank, const bool includeInterfaceParts) const;
  stk::mesh::PartVector get_conformal_parts_of_rank(const stk::mesh::EntityRank rank) const; // Not allowed for side_rank

  void determine_block_phases(const std::set<std::string> & FEmodel_block_names);
  void determine_block_phases();

  // ==========================================================================
  // TEMPORARY COMPATIBILITY SHIM -- DO NOT BUILD ON THIS.  REMOVE ME.
  // --------------------------------------------------------------------------
  // get_levelset_fields() used to be a static member of Phase_Support.  During
  // the krino library modularization it was extracted into a FREE function,
  //     krino::get_levelset_fields(const stk::mesh::MetaData &)
  // declared in Akri_LevelSet_Phase_Utils.hpp and defined in
  // Akri_LevelSet_Phase_Utils.cpp (both in krino_interface_geometry).
  //
  // External callers (platoengine's third_party_integration/krino, e.g.
  // Utilities.cpp) still spell it the old way, "::krino::Phase_Support::
  // get_levelset_fields(meta)", so their build breaks with:
  //   no member named 'get_levelset_fields' in 'krino::Phase_Support'
  //
  // This static member is re-added ONLY to unbreak those callers while they
  // migrate to the free function.  Note the deliberate layering trick that
  // keeps this legal:
  //   * krino_phase (this header) is BELOW krino_interface_geometry, so it
  //     cannot see the real LS_Field definition (Akri_LS_Field.hpp) or the
  //     free function.  We therefore only forward-declare LS_Field here; a
  //     by-value return type in a declaration needs no complete type.
  //   * The DEFINITION lives in Akri_LevelSet_Phase_Utils.cpp
  //     (krino_interface_geometry), which already sees both this class and the
  //     complete LS_Field, and simply forwards to the free function.  So the
  //     symbol resolves out of krino_interface_geometry and no build-system /
  //     dependency wiring changes are needed -- consumers already link
  //     krino_interface_geometry (via the krino_lib umbrella).
  //
  // The correct permanent fix is for callers to use the free function directly
  // (the compiler even suggests it): ::krino::get_levelset_fields(meta).
  // ==========================================================================
  static std::vector<LS_Field> get_levelset_fields(const stk::mesh::MetaData & meta);

  const NamedPhase & get_phase_with_name(const std::string & phaseName) const;
  void set_small_cluster_cleanup(const std::vector<std::tuple<std::string,std::string,unsigned>> & smallClusterCleanup);
  const auto & get_small_cluster_cleanup() const { return mySmallClusterCleanup; }

private:
  explicit Phase_Support(stk::mesh::MetaData & meta);

  const stk::mesh::MetaData & meta() const { STK_ThrowAssertMsg(myMeta, "MetaDeta not yet set on Phase_Support"); return *myMeta; }
  stk::mesh::MetaData & meta() { STK_ThrowAssertMsg(myMeta, "MetaDeta not yet set on Phase_Support"); return *myMeta; }
  const AuxMetaData & aux_meta() const { STK_ThrowAssertMsg(myAuxMeta, "AuxMetaData not yet set on Phase_Support"); return *myAuxMeta; }
  AuxMetaData & aux_meta() { STK_ThrowAssertMsg(myAuxMeta, "AuxMetaData not yet set on Phase_Support"); return *myAuxMeta; }

  void update_touching_parts_for_phase_part(const stk::mesh::Part & origPart, const stk::mesh::Part & phasePart, const PhaseTag & phase);
  // Phases the decomposed blocks and the top-level surfaces (one-sided sidesets and supersets) touching them.
  void create_phase_parts(const PhaseVec& phases, const PartSet& decomposedParts);
  // Phases the subsets of each superset, declaring the subset relations and aliases.
  void subset_and_alias_surface_phase_parts(const PhaseVec& phases, const PartSet& decomposedParts);
  std::string get_nonconforming_part_name(const stk::mesh::Part & originalPart) const;
  // The nonconforming part created by create_nonconforming_parts for an original part (error if missing).
  stk::mesh::Part & get_nonconforming_part_for_original_part(const stk::mesh::Part & originalPart) const;
  void create_nonconforming_parts(const PartSet & decomposed_ioparts);
  std::string get_phase_part_name(const stk::mesh::Part & originalPart, const NamedPhase & phase) const;
  // Declares the part for originalPart in the given phase (originalPart itself for an unnamed phase).
  stk::mesh::Part & declare_phase_part(const stk::mesh::Part & originalPart, const NamedPhase & phase);
  void add_conforming_phase_part(stk::mesh::Part & originalPart, const NamedPhase & phase);
  // Top-level sidesets (not subsets) that touch, directly or through a subset, a block in decomps.
  std::vector<unsigned> get_sidesets_touching_decomposed_blocks(const DecompositionPackage & decomps) const;
  // Blocks touching the sideset or any of its subsets.
  std::set<stk::mesh::PartOrdinal> get_blocks_touching_sideset(const stk::mesh::Part & sideset) const;
  // Sidesets touching blocks that are decomposed with different criteria for a phase with the same name.
  std::vector<unsigned> find_sidesets_with_conflicting_phase_criteria(const DecompositionPackage & decomps) const;
  void create_interface_phase_parts(const PhaseVec& ls_phases, const PartSet& decomposed_ioparts, const Interface_Name_Generator& interface_name_gen);

  void add_interface_part(const stk::mesh::Part & interfacePart,
      const stk::mesh::Part & nonconformingPart,
      const PhaseTag & touchingPhase,
      const PhaseTag & oppositePhase);
  void add_interface_superset_part(const stk::mesh::Part & interfacePart);

private:
  stk::mesh::MetaData * myMeta;
  AuxMetaData * myAuxMeta;
  Block_Surface_Connectivity my_input_block_surface_connectivity;
  PartPhaseMap my_mesh_block_phases;
  std::map<Surface_Identifier, IOPartSet> myLsUsedByParts; // blocks decomposed by each level set
  bool oneLevelSetPerPhase;
  PhaseVec myMeshPhases;
  PartnamePhasenameMap myMeshBlockPhasesByName;

  stk::mesh::Selector all_decomposed_blocks_selector;
  PhasePartInfo myPhasePartInfo;
  std::vector<std::tuple<NamedPhase,NamedPhase,unsigned>> mySmallClusterCleanup;
};

} // namespace krino

#endif // Akri_Phase_Support_h
