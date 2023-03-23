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

#include <Akri_Interface_Name_Generator.hpp>
#include <Akri_IO_Helpers.hpp>
#include <Akri_PhaseTag.hpp>

#include <map>

#include "Akri_FieldRef.hpp"
namespace stk { namespace mesh { class BulkData; } }
namespace stk { namespace mesh { class MetaData; } }
namespace stk { namespace diag { class Timer; } }

namespace krino {

class AuxMetaData;
class LevelSet;
class CDFEM_Inequality_Spec;

typedef std::set< const stk::mesh::Part * > IOPartSet;
typedef std::set< LevelSet * > LevelSetSet;
typedef std::map<stk::mesh::PartOrdinal,std::set<unsigned>> PartPhaseMap;
typedef std::map<std::string,std::vector<std::string>> PartnamePhasenameMap;

struct LS_Field
{
  LS_Field(const std::string & name_, const Surface_Identifier & identifier_, const FieldRef isovar_, const double isoval_, const LevelSet * const ptr_, const CDFEM_Inequality_Spec * const deathPtr_ = nullptr)
    : name(name_), identifier(identifier_), isovar(isovar_), isoval(isoval_), ptr(ptr_), deathPtr(deathPtr_) {
    ThrowRequireMsg(isovar_.valid(), "Invalid field " + isovar_.name() + " used in CDFEM initialization");
  }

  // Constructor just for unit tests
  LS_Field(const std::string & name_, const Surface_Identifier & identifier_)
    : name(name_), identifier(identifier_), isoval(0), ptr(nullptr), deathPtr(nullptr) {
  }

  std::string name;
  Surface_Identifier identifier;
  FieldRef isovar;
  double isoval;
  const LevelSet * ptr;
  const CDFEM_Inequality_Spec * deathPtr;
};

class Phase_Support {
public:
  Phase_Support (const Phase_Support&) = delete;
  Phase_Support& operator= (const Phase_Support&) = delete;

  typedef std::set<stk::mesh::Part*,stk::mesh::PartLess> PartSet;

  static bool exists_and_has_phases_defined(const stk::mesh::MetaData & meta);
  static Phase_Support & get_or_create(const std::string & FEModelName);
  static void associate_FEModel_and_metadata(const std::string & FEModelName, stk::mesh::MetaData & meta); // for typical parsing scenario where Phase_Support is originally created with the FEModel and now is being associated with MetaData
  static Phase_Support & get(const stk::mesh::MetaData & meta);
  static Phase_Support & get(stk::mesh::MetaData & meta); // for unit test usage where the Phase_Support is created with the MetaData

  static std::vector<LS_Field> get_levelset_fields(const stk::mesh::MetaData & meta);
  static void check_isovariable_field_existence_on_decomposed_blocks(const stk::mesh::MetaData & meta, const bool conformal_parts_require_field);
  static void check_isovariable_field_existence_on_decomposed_blocks(const stk::mesh::MetaData & meta, const std::vector<LS_Field> & lsFields, const bool conformal_parts_require_field);

  void get_blocks_touching_surface(const std::string & surface_name, std::vector<std::string> & block_names);

  std::vector<unsigned> get_negative_levelset_interface_ordinals(const Surface_Identifier levelSetIdentifier) const;

  void check_phase_parts() const;

  bool phases_defined() const { return !my_phase_parts.empty(); }

  void decompose_blocks(std::vector<std::tuple<stk::mesh::PartVector, std::shared_ptr<Interface_Name_Generator>, PhaseVec>> ls_sets);
  std::vector<stk::mesh::Part*> get_blocks_decomposed_by_levelset(const std::vector<unsigned> ls_phases) const;
  std::vector<stk::mesh::Part*> get_blocks_decomposed_by_bounding_surface(const unsigned surfaceID) const;
  void setup_phases();
  void set_input_block_surface_connectivity(const Block_Surface_Connectivity & input_block_surface_info) { my_input_block_surface_connectivity = input_block_surface_info; }
  const Block_Surface_Connectivity & get_input_block_surface_connectivity() const {return my_input_block_surface_connectivity;}
  void build_decomposed_block_surface_connectivity();

  bool is_nonconformal(const stk::mesh::Part * io_part) const;
  bool is_conformal(const stk::mesh::Part * io_part) const;
  bool is_interface(const stk::mesh::Part * io_part) const;
  const stk::mesh::Part * find_conformal_io_part(const stk::mesh::Part & io_part, const PhaseTag & phase) const;
  const stk::mesh::Part * find_nonconformal_part(const stk::mesh::Part & io_part) const;
  const stk::mesh::Part * find_interface_part(const stk::mesh::Part & vol0, const stk::mesh::Part & vol1) const;
  const PhaseTag & get_iopart_phase(const stk::mesh::Part & io_part) const;

  void add_decomposed_part(const stk::mesh::Part & part) { all_decomposed_blocks_selector |= part; }
  const stk::mesh::Selector & get_all_decomposed_blocks_selector() const { return all_decomposed_blocks_selector; }

  stk::mesh::Selector get_interface_part_selector(const LS_Field & ls_field);
  void register_blocks_for_level_set(const Surface_Identifier levelSetIdentifier,
      const std::vector<stk::mesh::Part *> & blocks_decomposed_by_ls);
  stk::mesh::Selector get_all_conformal_surfaces_selector() const;

  bool level_set_is_used_by_nonconformal_part(const Surface_Identifier levelSetIdentifier, const stk::mesh::Part * const ioPart) const;

  bool has_one_levelset_per_phase() const { return oneLevelSetPerPhase; }
  void set_one_levelset_per_phase(const bool val) { oneLevelSetPerPhase = val; }

  const PhaseVec & get_mesh_phases() const { return myMeshPhases; }
  PhaseVec & get_mesh_phases() { return myMeshPhases; }
  const PartnamePhasenameMap & get_block_phases_by_name() const { return myMeshBlockPhasesByName; }
  PartnamePhasenameMap & get_block_phases_by_name() { return myMeshBlockPhasesByName; }

  stk::mesh::PartVector get_nonconformal_parts() const;
  stk::mesh::PartVector get_nonconformal_parts_of_rank(const stk::mesh::EntityRank rank) const;
  stk::mesh::PartVector get_conformal_parts() const;
  stk::mesh::PartVector get_conformal_parts_of_rank(const stk::mesh::EntityRank rank) const;

  void determine_block_phases(const std::set<std::string> & FEmodel_block_names);
  void determine_block_phases();
  static PartSet get_blocks_and_touching_surfaces(const stk::mesh::MetaData & mesh_meta, const stk::mesh::PartVector& input_blocks, const Block_Surface_Connectivity & input_block_surface_info);

  static std::vector<unsigned> get_level_set_phases(const bool oneLSPerPhase, const PhaseVec & mesh_phases, const LevelSet & levelSet);

private:
  Phase_Support();
  Phase_Support(stk::mesh::MetaData & meta);

  const stk::mesh::MetaData & meta() const { ThrowAssertMsg(myMeta, "MetaDeta not yet set on Phase_Support"); return *myMeta; }
  stk::mesh::MetaData & meta() { ThrowAssertMsg(myMeta, "MetaDeta not yet set on Phase_Support"); return *myMeta; }
  const AuxMetaData & aux_meta() const { ThrowAssertMsg(myAuxMeta, "AuxMetaData not yet set on Phase_Support"); return *myAuxMeta; }
  AuxMetaData & aux_meta() { ThrowAssertMsg(myAuxMeta, "AuxMetaData not yet set on Phase_Support"); return *myAuxMeta; }

  void update_touching_parts_for_phase_part(const stk::mesh::Part & origPart, const stk::mesh::Part & phasePart, const PhaseTag & phase);
  const PhasePartTag * find_conformal_phase_part(const stk::mesh::Part & conformal_part) const;
  void create_nonconformal_parts(const PartSet & decomposed_ioparts);
  void addPhasePart(stk::mesh::Part & io_part, PhasePartSet & phase_parts, const NamedPhase & ls_phase);
  void create_phase_parts(const PhaseVec& ls_phases, const PartSet& decomposed_ioparts);
  void subset_and_alias_surface_phase_parts(const PhaseVec& ls_phases, const PartSet& decomposed_ioparts);
  void create_interface_phase_parts(const PhaseVec& ls_phases, const PartSet& decomposed_ioparts, const Interface_Name_Generator& interface_name_gen);
  void get_iopart_roots(const stk::mesh::Part & iopart, std::vector<const stk::mesh::Part *> & subsets);
  void fill_nonconformal_level_set_maps() const;

private:
  static std::map<std::string,std::unique_ptr<Phase_Support>> theModeltoPhaseSupportMap;
  stk::mesh::MetaData * myMeta;
  AuxMetaData * myAuxMeta;
  PhasePartSet my_phase_parts;
  Block_Surface_Connectivity my_input_block_surface_connectivity;
  PartPhaseMap my_mesh_block_phases;
  std::map<Surface_Identifier, IOPartSet> lsUsedByParts_;
  mutable bool nonconformalLsMapsAreFilled_;
  bool oneLevelSetPerPhase;
  PhaseVec myMeshPhases;
  PartnamePhasenameMap myMeshBlockPhasesByName;

  mutable std::map<const Surface_Identifier, IOPartSet> lsUsedByNonconformalParts_;

  typedef std::map< const stk::mesh::Part *, bool, stk::mesh::PartLess> PartToBoolMap;
  typedef std::map< const stk::mesh::Part *, const stk::mesh::Part *, stk::mesh::PartLess > PartToPartMap;
  typedef std::map<PhaseTag, const stk::mesh::Part *> PhaseTagToPartMap;
  typedef std::map<const stk::mesh::Part *, PhaseTagToPartMap, stk::mesh::PartLess> PartToPhaseTagToPartMap;
  typedef std::map<const stk::mesh::Part *,PhaseTag, stk::mesh::PartLess> PartToPhaseTagMap;
  typedef std::map<std::pair<const stk::mesh::Part *, const stk::mesh::Part *>,
      const stk::mesh::Part *> VolumePartsToInterfacePartMap;
  mutable PartToPhaseTagToPartMap nonconformal_to_phase_conformal_map;
  VolumePartsToInterfacePartMap volume_to_interface_parts_map;
  PartToBoolMap part_is_conformal_map;
  PartToBoolMap part_is_nonconformal_map;
  PartToPartMap part_to_nonconformal_part_map;
  PartToPhaseTagMap part_to_phase_map;

  stk::mesh::Selector all_decomposed_blocks_selector;
};

class CDFEM_Inequality_Spec {
public:

  CDFEM_Inequality_Spec(const std::string & name_);
  ~CDFEM_Inequality_Spec() {}
  const std::string & name() const {return my_name;}

  enum InequalityCriterionType
  {
    INEQUALITY_CRITERION_TYPE_UNDEFINED = 0,
    LESS_THAN                                 ,
    GREATER_THAN
  };

  static InequalityCriterionType int_to_inequality_criterion_type(const int inequality_criterion);

  void add_element_volume_name(const std::string & a_element_volue_name);
  bool set_threshold_variable_name(const std::string & a_threshold_variable_name);
  bool set_threshold_value(const double a_threshold_value);
  bool set_criterion_compare_type(const InequalityCriterionType & a_criterion_compare_type);

  void sanity_check(const std::vector<std::string> mesh_elem_blocks) const;
  const std::vector<std::string> & get_element_volume_names() const { return my_element_volume_names; }
  const std::string & get_threshold_variable_name () const { return my_threshold_variable_name; }
  double get_threshold_value () const { return my_threshold_value; }
  InequalityCriterionType get_criterion_compare_type() const { return my_criterion_compare_type; }
  const PhaseTag & get_deactivated_phase() const { return my_deactivated_phase; }
  const PhaseTag & get_active_phase() const { return my_active_phase; }
  void create_levelset(stk::mesh::MetaData & meta, stk::diag::Timer & parent_timer);
  const LevelSet & get_levelset() const { ThrowRequire(my_ls != NULL); return *my_ls; }
  LevelSet & get_levelset() { ThrowRequire(my_ls != NULL); return *my_ls; }

  // for unit testing
  void set_phases(const PhaseTag & active_phase, const PhaseTag & inactive_phase)
  {
    my_active_phase = active_phase;
    my_deactivated_phase = inactive_phase;
  }
protected:
  std::string my_name;
  std::vector<std::string> my_element_volume_names;
  std::string my_threshold_variable_name;
  double my_threshold_value;
  InequalityCriterionType my_criterion_compare_type;
  LevelSet * my_ls;
  PhaseTag my_deactivated_phase;
  PhaseTag my_active_phase;
};
typedef std::vector< CDFEM_Inequality_Spec > CDFEM_Inequality_Spec_Vec;

class CDFEM_Irreversible_Phase_Support {
public:
  ~CDFEM_Irreversible_Phase_Support() {}

  static CDFEM_Irreversible_Phase_Support * get_if_present(const stk::mesh::MetaData & meta);
  static CDFEM_Irreversible_Phase_Support & get(stk::mesh::MetaData & meta);
  static CDFEM_Irreversible_Phase_Support & get(const stk::mesh::MetaData & meta);

  CDFEM_Inequality_Spec * add_death_spec(const std::string & death_name, bool is_death);

  bool is_active() const {return has_death || has_irreversible_phase_change;}

  const CDFEM_Inequality_Spec_Vec & get_death_specs() const {return my_death_specs;}
  const CDFEM_Inequality_Spec * get_death_spec_for_ls(const LevelSet * ls) const;

  // For irreversible phase changes we decompose at the start of the time step, for death at the end.
  bool decompose_at_start_of_time_step() { ThrowAssert(has_irreversible_phase_change == !has_death); return has_irreversible_phase_change; }

private:
  CDFEM_Irreversible_Phase_Support() : has_death(false), has_irreversible_phase_change(false) {}
private:
  bool has_death;
  bool has_irreversible_phase_change;
  CDFEM_Inequality_Spec_Vec my_death_specs;
};

} // namespace krino

#endif // Akri_Phase_Support_h
