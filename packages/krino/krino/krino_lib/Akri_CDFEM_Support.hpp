// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#ifndef Akri_CDFEM_Support_h
#define Akri_CDFEM_Support_h
//
#include <Akri_LevelSet_Identifier.hpp>
#include <Akri_Phase_Support.hpp>
#include <Akri_FieldRef.hpp>
#include <Akri_AuxMetaData.hpp>
#include <Akri_CDFEM_Snapper.hpp>

#include <stk_util/diag/Timer.hpp>
#include <map>
#include <stk_mesh/base/Field.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/Selector.hpp>

namespace krino {

class LevelSet;
class CDFEM_Inequality_Spec;

struct LS_Field
{
  LS_Field(const std::string & name_, const LevelSet_Identifier & identifier_, const FieldRef isovar_, const double isoval_, const LevelSet * const ptr_)
    : name(name_), identifier(identifier_), isovar(isovar_), isoval(isoval_), ptr(ptr_) {
    ThrowRequireMsg(isovar_.valid(), "Invalid field " + isovar_.name() + " used in CDFEM initialization");
  }

  // Constructor just for unit tests
  LS_Field(const std::string & name_, const LevelSet_Identifier & identifier_)
    : name(name_), identifier(identifier_), isoval(0), ptr(NULL) {
  }

  std::string name;
  LevelSet_Identifier identifier;
  FieldRef isovar;
  double isoval;
  const LevelSet * ptr;
};

enum Prolongation_Model
{
  ALE_NEAREST_NODE=0,
  ALE_NEAREST_POINT,
  INTERPOLATION,
  MAX_PROLONGATION_MODEL
};

enum Edge_Interpolation_Model
{
  PIECEWISE_LINEAR=0,
  CONSTRAINED_LINEAR,
  MAX_EDGE_INTERPOLATION_MODEL
};

enum Edge_Degeneracy_Handling
{
  SNAP_TO_NODE=0,
  SNAP_TO_INTERFACE_WHEN_QUALITY_ALLOWS_THEN_SNAP_TO_NODE,
  MAX_EDGE_DEGENERACY_HANDLING_TYPE
};

enum Simplex_Generation_Method
{
  CUT_QUADS_BY_GLOBAL_IDENTIFIER=0,
  CUT_QUADS_BY_LARGEST_ANGLE,
  CUT_QUADS_BY_NEAREST_EDGE_CUT,
  MAX_SIMPLEX_GENERATION_METHOD
};

class CDFEM_Support {
public:

  CDFEM_Support(stk::mesh::MetaData & meta);
  ~CDFEM_Support() {}
  CDFEM_Support(CDFEM_Support const&) = delete;
  CDFEM_Support& operator=(CDFEM_Support const&) = delete;

  static CDFEM_Support & get(stk::mesh::MetaData & meta);
  static CDFEM_Support & get(const stk::mesh::MetaData & meta);
  static void use_constrained_edge_interpolation() { the_edge_interpolation_model = CONSTRAINED_LINEAR; }
  static Edge_Interpolation_Model get_edge_interpolation_model() { return the_edge_interpolation_model; }
  static bool use_nonconformal_adaptivity(stk::mesh::MetaData & meta) { CDFEM_Support & cdfem_support = get(meta); return cdfem_support.get_interface_maximum_refinement_level() > 0; }
  static std::string cdfem_mesh_displacements_field_name() { return "CDFEM_MESH_DISPLACEMENTS"; }

  static bool is_active(const stk::mesh::MetaData & meta)
  {
    return Phase_Support::exists_and_has_phases_defined(meta);
  }

  stk::mesh::MetaData & get_mesh_meta() { return my_meta; }
  const stk::mesh::MetaData & get_mesh_meta() const { return my_meta; }
  int num_ls_fields() const { return my_ls_fields.size(); }
  const LS_Field & ls_field(int i) const { return my_ls_fields[i]; }
  LS_Field & ls_field(int i) { return my_ls_fields[i]; }
  const std::vector<LS_Field> & ls_fields() const { return my_ls_fields; }
  Prolongation_Model get_prolongation_model() const { return my_prolongation_model; }
  void set_prolongation_model(const Prolongation_Model & model) { my_prolongation_model = model; }
  Simplex_Generation_Method get_simplex_generation_method() const { return my_simplex_generation_method; }
  void set_simplex_generation_method(const Simplex_Generation_Method & method);
  bool get_global_ids_are_parallel_consistent() const { return myGlobalIDsAreParallelConsistent; }
  void set_global_ids_are_NOT_parallel_consistent() { myGlobalIDsAreParallelConsistent = false; }
  const CDFEM_Inequality_Spec * get_death_spec(int ls_index) const { return my_death_specs[ls_index]; }
  void activate_interface_refinement(int minimum_level, int maximum_level);
  void activate_nonconformal_adaptivity(const int num_levels);
  int get_ls_index(const LevelSet * ls) const;

  void create_parts();

  void register_cdfem_mesh_displacements_field();
  void register_parent_node_ids_field();

  void setup_fields();
  void finalize_fields();
  void set_cdfem_displacement_field(const FieldRef field) { my_cdfem_displacements_field = field; }
  void set_cdfem_snap_displacement_field(const FieldRef field) { myCDFEMSnapDisplacementsField = field; }
  void add_ale_prolongation_field(const FieldRef field);
  void add_interpolation_field(const FieldRef field);
  // Just for unit test setup purposes
  void add_ls_field(const LS_Field & to_add, const CDFEM_Inequality_Spec * death_spec = nullptr) { my_ls_fields.push_back(to_add); my_death_specs.push_back(death_spec); }

  void set_coords_field(const FieldRef coords_field) { my_coords_field = coords_field; }
  const FieldRef get_coords_field() const { return my_coords_field; }
  const FieldRef get_cdfem_displacements_field() { return my_cdfem_displacements_field; }
  const FieldRef get_cdfem_snap_displacements_field() const { return myCDFEMSnapDisplacementsField; }
  const FieldSet & get_ale_prolongation_fields() const { return my_ale_prolongation_fields; }
  const FieldSet & get_interpolation_fields() const { return my_interpolation_fields; }
  const FieldSet & get_zeroed_fields() const { return my_zeroed_fields; }
  const FieldSet & get_element_fields() const { return my_element_fields; }

  bool add_initial_prolongation_field(const std::string & dest_field_name, const std::string & src_field_name);
  FieldRef get_initial_prolongation_field(const FieldRef field) const;

  stk::mesh::Part & get_parent_part() const { return *my_parent_part; }
  stk::mesh::Part & get_child_part() const { return *my_child_part; }
  stk::mesh::Part & get_internal_side_part() const { return *my_internal_side_part; }
  stk::mesh::Part & get_active_part() const { return my_aux_meta.active_part(); }
  stk::mesh::Part & get_universal_part() const { return my_meta.universal_part(); }
  stk::mesh::Part & get_locally_owned_part() const { return my_meta.locally_owned_part(); }
  stk::mesh::Part & get_globally_shared_part() const { return my_meta.globally_shared_part(); }

  stk::mesh::Selector get_post_cdfem_refinement_selector() const;

  stk::mesh::Part & get_child_edge_node_part() const { return *my_child_edge_node_part; }
  FieldRef get_parent_node_ids_field() const { return my_parent_node_ids_field; }

  void activate_fully_coupled_cdfem() { my_fully_coupled_cdfem = true; }
  bool fully_coupled_cdfem() const { return my_fully_coupled_cdfem; }

  void activate_nonconformal_adapt_target_count(uint64_t val);
  uint64_t get_nonconformal_adapt_target_count() const { return my_nonconformal_adapt_target_element_count; }
  int get_interface_minimum_refinement_level() const { return my_interface_minimum_refinement_level; }
  int get_interface_maximum_refinement_level() const { return my_interface_maximum_refinement_level; }
  void set_post_adapt_refinement_levels(int levels) { my_post_adapt_uniform_refinement_levels = levels; }
  int get_post_adapt_refinement_levels() const { return my_post_adapt_uniform_refinement_levels; }
  void set_post_cdfem_refinement_levels(int levels) { my_post_cdfem_refinement_levels = levels; }
  int get_post_cdfem_refinement_levels() const { return my_post_cdfem_refinement_levels; }
  void set_post_cdfem_refinement_blocks(const std::vector<std::string> & post_cdfem_refinement_blocks) { my_post_cdfem_refinement_blocks = post_cdfem_refinement_blocks; }
  void set_num_initial_decomposition_cycles(int num_initial_decomposition_cycles) { my_num_initial_decomposition_cycles = num_initial_decomposition_cycles; }
  // In case of both nonconformal adaptivity, perform at least 2 rounds of level set initialization and decomposition.
  // This is to capture features that might be missed on the unrefined mesh.
  int get_num_initial_decomposition_cycles() const { return (my_num_initial_decomposition_cycles > 1) ? my_num_initial_decomposition_cycles : ((my_interface_maximum_refinement_level > 0) ? 2 : 1); }
  const std::string & get_nonconformal_adapt_marker_name() const { return my_nonconformal_adapt_marker_name; }
  const std::string & get_nonconformal_adapt_indicator_name() const { return my_nonconformal_adapt_indicator_name; }
  void set_nonconformal_hadapt(const std::function<void(const std::string &, int)> & hadapt) { my_nonconformal_hadapt = hadapt; }
  const std::function<void(const std::string &, int)> & get_nonconformal_hadapt() const { return my_nonconformal_hadapt; }

  bool is_ale_prolongation_field(const FieldRef field) const
  {
    return (my_ale_prolongation_fields.find(field) != my_ale_prolongation_fields.end());
  }
  bool is_interpolation_field(const FieldRef field) const
  {
    return (my_interpolation_fields.find(field) != my_interpolation_fields.end());
  }

  void use_nonconformal_element_size(bool flag) { my_flag_use_nonconformal_element_size = flag; }
  bool use_nonconformal_element_size() const { return my_flag_use_nonconformal_element_size; }

  Edge_Degeneracy_Handling get_cdfem_edge_degeneracy_handling() const { return my_cdfem_edge_degeneracy_handling; }
  void set_cdfem_edge_degeneracy_handling( const Edge_Degeneracy_Handling type ) { my_cdfem_edge_degeneracy_handling = type; }

  stk::diag::Timer & get_timer_cdfem() const { return my_timer_cdfem; }
  stk::diag::Timer & get_timer_adapt() const { return my_timer_adapt; }

  void set_cdfem_edge_tol( const double tol ) { my_cdfem_snapper.set_edge_tolerance(tol); }
  const CDFEM_Snapper & get_snapper() const { return my_cdfem_snapper; }
  const double & get_cdfem_dof_edge_tol() const { return my_cdfem_dof_edge_tol; }
  void set_cdfem_dof_edge_tol( const double tol ) { my_cdfem_dof_edge_tol = tol; }
  bool use_internal_face_stabilization() const { return my_internal_face_stabilization_multiplier > 0.0; }
  double get_internal_face_stabilization_multiplier() const { return my_internal_face_stabilization_multiplier; }
  void set_internal_face_stabilization_multiplier( const double mult ) { my_internal_face_stabilization_multiplier = mult; }
  bool get_use_hierarchical_dofs() const { return my_flag_use_hierarchical_dofs; }
  void set_use_hierarchical_dofs(bool flag) { my_flag_use_hierarchical_dofs = flag; }
  bool get_constrain_CDFEM_to_XFEM_space() const { return my_flag_constrain_CDFEM_to_XFEM_space; }
  void set_constrain_CDFEM_to_XFEM_space(bool flag) { my_flag_constrain_CDFEM_to_XFEM_space = flag; }

  void force_ale_prolongation_for_field(const std::string & field_name);

private:
  void setup_refinement_marker();

private:
  stk::mesh::MetaData & my_meta;
  AuxMetaData & my_aux_meta;

  FieldRef my_coords_field;
  std::vector<LS_Field> my_ls_fields;
  std::vector<const CDFEM_Inequality_Spec *> my_death_specs;
  FieldRef my_cdfem_displacements_field;
  FieldRef myCDFEMSnapDisplacementsField;
  Prolongation_Model my_prolongation_model;
  Simplex_Generation_Method my_simplex_generation_method;
  std::vector<std::string> my_force_ale_prolongation_fields;
  FieldSet my_ale_prolongation_fields;
  FieldSet my_interpolation_fields;
  FieldSet my_zeroed_fields;
  FieldSet my_element_fields;
  std::map<std::string, std::string> my_initial_prolongation_field_name_map;
  std::map<FieldRef, FieldRef> my_initial_prolongation_field_map;
  static Edge_Interpolation_Model the_edge_interpolation_model;
  stk::mesh::Part * my_parent_part;
  stk::mesh::Part * my_child_part;
  stk::mesh::Part * my_internal_side_part;
  stk::mesh::Part * my_child_edge_node_part;
  FieldRef my_parent_node_ids_field;
  bool my_fully_coupled_cdfem;
  int my_num_initial_decomposition_cycles;
  bool myGlobalIDsAreParallelConsistent;
  int my_interface_minimum_refinement_level;
  int my_interface_maximum_refinement_level;
  int my_post_adapt_uniform_refinement_levels;
  int my_post_cdfem_refinement_levels;
  std::vector<std::string> my_post_cdfem_refinement_blocks;
  uint64_t my_nonconformal_adapt_target_element_count;
  std::string my_nonconformal_adapt_marker_name;
  std::string my_nonconformal_adapt_indicator_name;
  std::function<void(const std::string &, int)> my_nonconformal_hadapt;
  Edge_Degeneracy_Handling my_cdfem_edge_degeneracy_handling;
  CDFEM_Snapper my_cdfem_snapper;
  double my_cdfem_dof_edge_tol;
  double my_internal_face_stabilization_multiplier;
  bool my_flag_use_hierarchical_dofs;
  bool my_flag_constrain_CDFEM_to_XFEM_space;
  bool my_flag_use_nonconformal_element_size;
  mutable stk::diag::Timer my_timer_cdfem;
  mutable stk::diag::Timer my_timer_adapt;
};

} // namespace krino

#endif // Akri_CDFEM_Support_h
