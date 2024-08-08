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
#include "Akri_CDFEM_Support.hpp"

#include <Akri_FieldRef.hpp>
#include <Akri_AuxMetaData.hpp>
#include <Akri_CDFEM_Snapper.hpp>

#include <stk_util/diag/Timer.hpp>
#include <map>
#include <stk_mesh/base/Field.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/Selector.hpp>

namespace krino {

class Phase_Support;
class RefinementInterface;

enum Prolongation_Model
{
  ALE_CLOSEST_POINT=0,
  ALE_CLOSEST_NODE,
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

enum Resnap_Method
{
  RESNAP_AFTER_USING_ALE_TO_UNSNAP=0,
  RESNAP_AFTER_USING_INTERPOLATION_TO_UNSNAP,
  RESNAP_USING_INTERPOLATION,
  RESNAP_USING_INTERFACE_ON_PREVIOUS_SNAPPED_MESH,
  MAX_RESNAP_METHOD
};

enum Simplex_Generation_Method
{
  CUT_QUADS_BY_GLOBAL_IDENTIFIER=0,
  CUT_QUADS_BY_LARGEST_ANGLE,
  CUT_QUADS_BY_NEAREST_EDGE_CUT,
  MAX_SIMPLEX_GENERATION_METHOD
};

enum Interface_CFL_Length_Scale
{
  CONSTANT_LENGTH_SCALE=0,
  LOCAL_LENGTH_SCALE,
  L1_NORM_LENGTH_SCALE,
  MAX_LENGTH_SCALE_TYPE
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
  static std::string cdfem_mesh_displacements_field_name() { return "CDFEM_MESH_DISPLACEMENTS"; }
  static std::string cdfem_snap_displacements_field_name() { return "CDFEM_SNAP_DISPLACEMENTS"; }

  static bool is_active(const stk::mesh::MetaData & meta);

  stk::mesh::MetaData & get_mesh_meta() { return my_meta; }
  const stk::mesh::MetaData & get_mesh_meta() const { return my_meta; }
  Prolongation_Model get_prolongation_model() const { return my_prolongation_model; }
  void set_prolongation_model(const Prolongation_Model & model) { my_prolongation_model = model; }
  Simplex_Generation_Method get_simplex_generation_method() const { return my_simplex_generation_method; }
  void set_simplex_generation_method(const Simplex_Generation_Method & method);
  bool get_global_ids_are_parallel_consistent() const { return myGlobalIDsAreParallelConsistent; }
  void set_global_ids_are_NOT_parallel_consistent() { myGlobalIDsAreParallelConsistent = false; }
  void set_snapping_sharp_feature_angle_in_degrees(const double snappingSharpFeatureAngleInDegrees) { mySnappingSharpFeatureAngleInDegrees = snappingSharpFeatureAngleInDegrees; }
  double get_snapping_sharp_feature_angle_in_degrees() const { return mySnappingSharpFeatureAngleInDegrees; }
  void set_is_transient(const bool isTransient) { myFlagIsTransient = isTransient; }
  bool get_is_transient() const { return myFlagIsTransient; }

  void create_parts();

  void register_cdfem_mesh_displacements_field();
  void register_cdfem_snap_displacements_field();
  void register_parent_node_ids_field();

  void setup_fields();
  void finalize_fields();
  void set_cdfem_displacement_field(const FieldRef field) { my_cdfem_displacements_field = field; }
  void set_cdfem_snap_displacement_field(const FieldRef field) { myCDFEMSnapDisplacementsField = field; }
  void add_ale_prolongation_field(const FieldRef field);
  void add_interpolation_field(const FieldRef field);
  void add_edge_interpolation_field(const FieldRef field);

  void set_coords_field(const FieldRef coords_field) { my_coords_field = coords_field; }
  const FieldRef get_coords_field() const { return my_coords_field; }
  const FieldRef get_cdfem_displacements_field() { return my_cdfem_displacements_field; }
  const FieldRef get_cdfem_snap_displacements_field() const { return myCDFEMSnapDisplacementsField; }
  const FieldSet & get_ale_prolongation_fields() const { return my_ale_prolongation_fields; }
  const FieldSet & get_interpolation_fields() const { return my_interpolation_fields; }
  const FieldSet & get_edge_interpolation_fields() const { return my_edge_interpolation_fields; }
  const FieldSet & get_zeroed_fields() const { return my_zeroed_fields; }
  const FieldSet & get_element_fields() const { return my_element_fields; }
  const FieldSet & get_snap_fields() const { return mySnapFields; }
  const FieldSet & get_levelset_fields() const { return myLevelSetFields; }

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

  stk::mesh::Part & get_child_node_part() const { return *myChildNodePart; }
  FieldRef get_parent_node_ids_field() const { return myParentNodeIdsField; }
  FieldRef get_parent_node_weights_field() const { return myParentNodeWtsField; }

  void activate_fully_coupled_cdfem() { my_fully_coupled_cdfem = true; }
  bool fully_coupled_cdfem() const { return my_fully_coupled_cdfem; }

  void set_post_cdfem_refinement_levels(int levels) { my_post_cdfem_refinement_levels = levels; }
  int get_post_cdfem_refinement_levels() const { return my_post_cdfem_refinement_levels; }
  void set_post_cdfem_refinement_blocks(const std::vector<std::string> & post_cdfem_refinement_blocks) { my_post_cdfem_refinement_blocks = post_cdfem_refinement_blocks; }
  void set_num_initial_decomposition_cycles(int num_initial_decomposition_cycles) { my_num_initial_decomposition_cycles = num_initial_decomposition_cycles; }
  int get_num_initial_decomposition_cycles() const { return (my_num_initial_decomposition_cycles > 1) ? my_num_initial_decomposition_cycles : 1; }

  bool is_ale_prolongation_field(const FieldRef field) const
  {
    return (my_ale_prolongation_fields.find(field) != my_ale_prolongation_fields.end());
  }
  bool is_interpolation_field(const FieldRef field) const
  {
    return (my_interpolation_fields.find(field) != my_interpolation_fields.end());
  }
  bool is_edge_interpolation_field(const FieldRef field) const
  {
    return (my_edge_interpolation_fields.find(field) != my_edge_interpolation_fields.end());
  }

  void use_nonconformal_element_size(bool flag) { my_flag_use_nonconformal_element_size = flag; }
  bool use_nonconformal_element_size() const { return my_flag_use_nonconformal_element_size; }

  void set_use_facets_instead_of_levelset_fields(bool flag) { myFlagUseFacetsInsteadOfLsFields = flag; }
  bool use_facets_instead_of_levelset_fields() const { return myFlagUseFacetsInsteadOfLsFields; }

  Edge_Degeneracy_Handling get_cdfem_edge_degeneracy_handling() const { return my_cdfem_edge_degeneracy_handling; }
  void set_cdfem_edge_degeneracy_handling( const Edge_Degeneracy_Handling type );

  stk::diag::Timer & get_timer_cdfem() const { return my_timer_cdfem; }

  void set_cdfem_edge_tol( const double tol ) { my_cdfem_snapper.set_edge_tolerance(tol); }
  const CDFEM_Snapper & get_snapper() const { return my_cdfem_snapper; }
  double get_cdfem_dof_edge_tol() const { return my_cdfem_dof_edge_tol; }
  void set_cdfem_dof_edge_tol( const double tol ) { my_cdfem_dof_edge_tol = tol; }
  double get_max_edge_snap() const { return myMaxEdgeSnap; }
  void set_max_edge_snap( const double snap ) { myMaxEdgeSnap = snap; }
  bool use_internal_face_stabilization() const { return my_internal_face_stabilization_multiplier > 0.0; }
  double get_internal_face_stabilization_multiplier() const { return my_internal_face_stabilization_multiplier; }
  void set_internal_face_stabilization_multiplier( const double mult ) { my_internal_face_stabilization_multiplier = mult; }
  bool get_use_hierarchical_dofs() const { return my_flag_use_hierarchical_dofs; }
  void set_use_hierarchical_dofs(bool flag) { my_flag_use_hierarchical_dofs = flag; }
  bool get_constrain_CDFEM_to_XFEM_space() const { return my_flag_constrain_CDFEM_to_XFEM_space; }
  void set_constrain_CDFEM_to_XFEM_space(bool flag) { my_flag_constrain_CDFEM_to_XFEM_space = flag; }
  void set_resnap_method(Resnap_Method resnapMethod) { myResnapMethod = resnapMethod; }
  Resnap_Method get_resnap_method() const { return myResnapMethod; }
  void set_perform_volume_correction_after_levelset_solve(bool flag) { myFlagPerformVolumeCorrectionAfterLsSolve = flag; }
  bool get_perform_volume_correction_after_levelset_solve() const { return myFlagPerformVolumeCorrectionAfterLsSolve; }

  void set_constant_length_scale_for_interface_CFL(double lengthScale) { myConstantLengthScaleForInterfaceCFL = lengthScale; }
  double get_constant_length_scale_for_interface_CFL() const { return myConstantLengthScaleForInterfaceCFL; }
  void set_length_scale_type_for_interface_CFL(Interface_CFL_Length_Scale lengthScaleType) { myLengthScaleTypeForInterfaceCFL = lengthScaleType; }
  Interface_CFL_Length_Scale get_length_scale_type_for_interface_CFL() const { return myLengthScaleTypeForInterfaceCFL; }
  bool get_use_velocity_to_evaluate_interface_CFL() const { return myFlagUseVelocityToEvaluateInterfaceCFL; }
  void set_use_velocity_to_evaluate_interface_CFL(bool flag) { myFlagUseVelocityToEvaluateInterfaceCFL = flag; }

  void force_ale_prolongation_for_field(const std::string & field_name);

  void add_user_prolongation_field(const std::vector<std::string> & fields)
  {
    my_user_prolongation_field_names.insert(fields.begin(), fields.end());
  }
  
  const std::set<std::string> & get_user_prolongation_field() const
  {
    return my_user_prolongation_field_names;
  }


private:
  void set_snap_fields();

private:
  stk::mesh::MetaData & my_meta;
  AuxMetaData & my_aux_meta;

  FieldRef my_coords_field;
  FieldRef my_cdfem_displacements_field;
  FieldRef myCDFEMSnapDisplacementsField;
  Prolongation_Model my_prolongation_model;
  Simplex_Generation_Method my_simplex_generation_method;
  std::vector<std::string> my_force_ale_prolongation_fields;
  FieldSet my_ale_prolongation_fields;
  FieldSet my_interpolation_fields;
  FieldSet my_edge_interpolation_fields;
  FieldSet my_zeroed_fields;
  FieldSet my_element_fields;
  FieldSet mySnapFields;
  FieldSet myLevelSetFields;
  std::map<std::string, std::string> my_initial_prolongation_field_name_map;
  std::map<FieldRef, FieldRef> my_initial_prolongation_field_map;
  static Edge_Interpolation_Model the_edge_interpolation_model;
  stk::mesh::Part * my_parent_part;
  stk::mesh::Part * my_child_part;
  stk::mesh::Part * my_internal_side_part;
  stk::mesh::Part * myChildNodePart;
  FieldRef myParentNodeIdsField;
  FieldRef myParentNodeWtsField;
  bool my_fully_coupled_cdfem;
  int my_num_initial_decomposition_cycles;
  bool myGlobalIDsAreParallelConsistent;
  int my_post_cdfem_refinement_levels;
  std::vector<std::string> my_post_cdfem_refinement_blocks;
  Edge_Degeneracy_Handling my_cdfem_edge_degeneracy_handling;
  CDFEM_Snapper my_cdfem_snapper;
  double my_cdfem_dof_edge_tol;
  double myMaxEdgeSnap{1.0};
  double my_internal_face_stabilization_multiplier;
  double mySnappingSharpFeatureAngleInDegrees;
  Interface_CFL_Length_Scale myLengthScaleTypeForInterfaceCFL;
  double myConstantLengthScaleForInterfaceCFL;
  bool my_flag_use_hierarchical_dofs;
  bool my_flag_constrain_CDFEM_to_XFEM_space;
  bool my_flag_use_nonconformal_element_size;
  bool myFlagUseFacetsInsteadOfLsFields{false};
  bool myFlagUseVelocityToEvaluateInterfaceCFL;
  Resnap_Method myResnapMethod{RESNAP_AFTER_USING_ALE_TO_UNSNAP};
  bool myFlagPerformVolumeCorrectionAfterLsSolve{false};
  bool myFlagIsTransient{false};
  mutable stk::diag::Timer my_timer_cdfem;
  std::set<std::string> my_user_prolongation_field_names;
};

} // namespace krino

#endif // Akri_CDFEM_Support_h
