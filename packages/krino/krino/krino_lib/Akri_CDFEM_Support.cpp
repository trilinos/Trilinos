// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include <Akri_CDFEM_Support.hpp>
#include <Akri_LevelSet.hpp>
#include <Akri_DiagWriter.hpp>
#include <stk_util/diag/Timer.hpp>
#include <stk_util/environment/RuntimeDoomed.hpp>

namespace krino{

Edge_Interpolation_Model CDFEM_Support::the_edge_interpolation_model = PIECEWISE_LINEAR;

CDFEM_Support &
CDFEM_Support::get(stk::mesh::MetaData & meta)
{
  CDFEM_Support * support = const_cast<CDFEM_Support *>(meta.get_attribute<CDFEM_Support>());
  if (NULL == support)
  {
    support = new CDFEM_Support(meta);
    meta.declare_attribute_with_delete<CDFEM_Support>(support);
  }
  return *support;
}

CDFEM_Support &
CDFEM_Support::get(const stk::mesh::MetaData & meta)
{
  CDFEM_Support * support = const_cast<CDFEM_Support *>(meta.get_attribute<CDFEM_Support>());
  ThrowRequireMsg(nullptr != support, "Could not find CDFEM_Support attribute on MetaData.");
  return *support;
}

CDFEM_Support::CDFEM_Support(stk::mesh::MetaData & meta)
  : my_meta(meta),
    my_aux_meta(AuxMetaData::get(meta)),
    my_simplex_generation_method(CUT_QUADS_BY_LARGEST_ANGLE),
    my_fully_coupled_cdfem(false),
    my_num_initial_decomposition_cycles(1),
    myGlobalIDsAreParallelConsistent(true),
    my_interface_minimum_refinement_level(0),
    my_interface_maximum_refinement_level(0),
    my_post_adapt_uniform_refinement_levels(0),
    my_post_cdfem_refinement_levels(0),
    my_nonconformal_adapt_target_element_count(0),
    my_cdfem_edge_degeneracy_handling(SNAP_TO_NODE),
    my_cdfem_snapper(),
    my_cdfem_dof_edge_tol(0.0),
    my_internal_face_stabilization_multiplier(0.0),
    my_flag_use_hierarchical_dofs(false),
    my_flag_constrain_CDFEM_to_XFEM_space(false),
    my_flag_use_nonconformal_element_size(true),
    my_timer_cdfem("CDFEM", sierra::Diag::sierraTimer()),
    my_timer_adapt("Nonconformal Adapt", my_timer_cdfem)
{
  my_prolongation_model = ALE_NEAREST_POINT;

  if (3 == my_meta.spatial_dimension())
    my_simplex_generation_method = CUT_QUADS_BY_NEAREST_EDGE_CUT;

  create_parts();
}

void CDFEM_Support::create_parts()
{
  my_parent_part = &my_meta.declare_part("CDFEM_PARENT_CONTEXT_BIT");
  my_child_part = &my_meta.declare_part("CDFEM_CHILD_CONTEXT_BIT");
  my_internal_side_part = &my_meta.declare_part("CDFEM_INTERNAL_SIDE");

  if (my_aux_meta.using_fmwk())
  {
    const bool restartOnlyIOPart = true;
    my_child_edge_node_part = &my_aux_meta.declare_io_part_with_topology("CDFEM_EDGE_NODE_2_PARENTS", stk::topology::NODE, restartOnlyIOPart);
  }
  else
  {
    // Currently no need to output nodeset for krino usage
    my_child_edge_node_part = &my_meta.declare_part_with_topology("CDFEM_EDGE_NODE_2_PARENTS", stk::topology::NODE);
  }
}

void CDFEM_Support::register_cdfem_mesh_displacements_field()
{
  ThrowRequireMsg(my_meta.spatial_dimension() > 1, "Spatial dimension must be set and equal to 2 or 3.");

  const FieldType & vec_type = (my_meta.spatial_dimension() == 3) ? FieldType::VECTOR_3D : FieldType::VECTOR_2D;
  FieldRef cdfem_disp_field = my_aux_meta.register_field(cdfem_mesh_displacements_field_name(), vec_type, stk::topology::NODE_RANK, 2, 1, get_universal_part());
  set_cdfem_displacement_field(cdfem_disp_field);
}

void CDFEM_Support::register_parent_node_ids_field()
{
  FieldType id_field_type = (my_aux_meta.get_assert_32bit_flag()) ? FieldType::UNSIGNED_INTEGER : FieldType::UNSIGNED_INTEGER_64;
  my_parent_node_ids_field = my_aux_meta.register_field("CDFEM_2_PARENT_NODE_IDS",
      id_field_type, stk::topology::NODE_RANK, 1,
      2, *my_child_edge_node_part);
}

stk::mesh::Selector
CDFEM_Support::get_post_cdfem_refinement_selector() const
{
  stk::mesh::Selector selector = my_aux_meta.active_part();
  if (!my_post_cdfem_refinement_blocks.empty())
  {
    stk::mesh::PartVector blocks;
    for (auto && block_name : my_post_cdfem_refinement_blocks)
    {
      if (my_aux_meta.has_part(block_name))
      {
        stk::mesh::Part & block = my_aux_meta.get_part(block_name);
        blocks.push_back(&block);
      }
      else
      {
        stk::RuntimeDoomedAdHoc() << "post_cdfem_refinement_block " << block_name << " not found.\n";
      }
    }
    selector &= stk::mesh::selectUnion(blocks);
  }
  return selector;
}

void
CDFEM_Support::setup_fields()
{
  my_coords_field = LevelSet::get_current_coordinates(my_meta);

  Phase_Support::get(my_meta).check_phase_parts();

  my_ls_fields.clear();
  const LevelSetManager & region_ls = LevelSetManager::get(my_meta);
  for (auto&& ls : region_ls)
  {
    LS_Field calc_ls_field(ls->name(),ls->get_identifier(),ls->get_isovar_field(),ls->get_isoval(),ls.get());
    my_ls_fields.push_back(calc_ls_field);
  }

  const unsigned num_ls = region_ls.numberLevelSets();
  my_death_specs.resize(num_ls, nullptr);
  for (unsigned i=0; i<num_ls; ++i)
  {
    const LevelSet & ls = LevelSetManager::get(my_meta).levelSet(i);

    if ( CDFEM_Irreversible_Phase_Support::get(my_meta).is_active() )
    {
      const CDFEM_Inequality_Spec_Vec & death_specs = CDFEM_Irreversible_Phase_Support::get(my_meta).get_death_specs();

      for (auto && death_spec : death_specs)
      {
        if (&death_spec.get_levelset() == &ls)
        {
          my_death_specs[i] = &death_spec;
        }
      }
    }
  }
}

void
CDFEM_Support::add_ale_prolongation_field(const FieldRef field)
{
  ThrowAssert(field.valid());
  ThrowRequireMsg(!is_interpolation_field(field), "Cannot add " << field.name() << " as ALE prolongation field because it is already an interpolation field.");
  for ( unsigned is = 0; is < field.number_of_states(); ++is )
  {
    const stk::mesh::FieldState state = static_cast<stk::mesh::FieldState>(is);
    my_ale_prolongation_fields.insert(field.field_state(state));
  }
}

void
CDFEM_Support::add_interpolation_field(const FieldRef field)
{
  ThrowAssert(field.valid());
  ThrowRequireMsg(!is_ale_prolongation_field(field), "Cannot add " << field.name() << " as interpolation field because it is already an ALE prolongation field.");
  for ( unsigned is = 0; is < field.number_of_states(); ++is )
  {
    const stk::mesh::FieldState state = static_cast<stk::mesh::FieldState>(is);
    my_interpolation_fields.insert(field.field_state(state));
  }
}

void
CDFEM_Support::finalize_fields()
{
  krinolog << "Finalizing field prolongation strategies for CDFEM." << stk::diag::push << stk::diag::dendl;
  for ( auto && field_ptr : my_meta.get_fields() )
  {
    const FieldRef field(field_ptr);
    if( !field.type_is<double>() || field.field_state(stk::mesh::StateNew) == my_cdfem_displacements_field ) continue;

    if( field.entity_rank()==stk::topology::ELEMENT_RANK &&
        !stk::equal_case(field.name(), "transition_element") &&
        !stk::equal_case(field.name(), "transition_element_3") &&
        !stk::equal_case(field.name(), "parent_element") )
    {
      my_element_fields.insert(field);
    }

    if( field.entity_rank()!=stk::topology::NODE_RANK ) continue;

    const auto& initial_prolong_field_name_entry = my_initial_prolongation_field_name_map.find(field.name());
    if (initial_prolong_field_name_entry != my_initial_prolongation_field_name_map.end())
    {
      const std::string & src_field_name = initial_prolong_field_name_entry->second;
      ThrowErrorMsgIf(!my_aux_meta.has_field(stk::topology::NODE_RANK, src_field_name),
          "Error: Could not find initial prolongation field with name " << src_field_name);

      // If the src field does not have the desired state, use StateNone (which is the same as StateNew).
      FieldRef src_field = my_aux_meta.get_field(stk::topology::NODE_RANK, src_field_name);
      stk::mesh::FieldState src_state = (field.state() < src_field.number_of_states()) ? field.state() : stk::mesh::StateNone;
      src_field = src_field.field_state(src_state);

      my_initial_prolongation_field_map[field] = src_field;
      krinolog << "Added "
          << src_field.name() << " (" << state_string(src_state) << ") as initial prolongation field for "
          << field.name()  << " (" << state_string(field.state()) << ")" << stk::diag::dendl;
    }

    if (std::find(my_force_ale_prolongation_fields.begin(),
            my_force_ale_prolongation_fields.end(),
            field.name()) != my_force_ale_prolongation_fields.end())
    {
      add_ale_prolongation_field(field);
    }

    if(is_interpolation_field(field))
    {
      krinolog << field.name() << " will use interpolation." << stk::diag::dendl;
      continue;
    }
    else if(is_ale_prolongation_field(field))
    {
      krinolog << field.name() << " will use ALE prolongation." << stk::diag::dendl;
      continue;
    }
    else if (field.name() == "node_registry")
    {
      krinolog << field.name() << " will not be modified." << stk::diag::dendl;
      continue;
    }

    krinolog << field.name() << " will be zeroed." << stk::diag::dendl;
    my_zeroed_fields.insert(field);
  }
  krinolog << stk::diag::pop << stk::diag::dendl;
}

void CDFEM_Support::force_ale_prolongation_for_field(const std::string & field_name)
{
  my_force_ale_prolongation_fields.push_back(field_name);
}

bool
CDFEM_Support::add_initial_prolongation_field(const std::string & dest_field_name, const std::string & src_field_name)
{
  if (my_initial_prolongation_field_name_map.find(dest_field_name) != my_initial_prolongation_field_name_map.end())
  {
    return false;
  }
  my_initial_prolongation_field_name_map[dest_field_name] = src_field_name;
  return true;
}

FieldRef
CDFEM_Support::get_initial_prolongation_field(const FieldRef field) const
{
  auto iter = my_initial_prolongation_field_map.find(field);
  if (iter == my_initial_prolongation_field_map.end())
    return FieldRef();
  return iter->second;
}

void
CDFEM_Support::set_simplex_generation_method(const Simplex_Generation_Method & method)
{
  ThrowAssert(method < MAX_SIMPLEX_GENERATION_METHOD);
  ThrowRequireMsg(3 == my_meta.spatial_dimension() || method != CUT_QUADS_BY_NEAREST_EDGE_CUT, "Simplex generation method CUT_QUADS_BY_NEAREST_EDGE_CUT only supported in 3d.");
  my_simplex_generation_method = method;
}

void
CDFEM_Support::activate_interface_refinement(int minimumLevel, int maximumLevel)
{
  /* %TRACE% */ Traceback trace__("krino::CDFEM_Support::activate_interface_refinement(int minimum_level, int maximum_level)"); /* %TRACE% */

  ThrowRequireMsg(my_interface_minimum_refinement_level == 0 && my_interface_maximum_refinement_level == 0,
      "Interface refinement levels should only be specified once.");
  ThrowRequireMsg(maximumLevel >= minimumLevel || maximumLevel == 0,
      "Maximum interface refinement level must be greater than or equal to the minimum interface refinement level or left unspecified.");
  if (maximumLevel == 0) maximumLevel = minimumLevel;

  my_interface_minimum_refinement_level = minimumLevel;
  my_interface_maximum_refinement_level = maximumLevel;

  setup_refinement_marker();

  if (maximumLevel > 0)
    set_global_ids_are_NOT_parallel_consistent();
}

void
CDFEM_Support::activate_nonconformal_adaptivity(const int numLevels)
{
  /* %TRACE% */ Traceback trace__("krino::CDFEM_Support::activate_nonconformal_adaptivity(const int num_levels)"); /* %TRACE% */

  if (numLevels < my_interface_maximum_refinement_level)
  {
    krinolog << "Ignoring request to activate " << numLevels << " of CDFEM nonconformal adaptivity because a maximum of " << my_interface_maximum_refinement_level << " have already been activated." << stk::diag::dendl;
    return;
  }

  my_interface_minimum_refinement_level = numLevels;
  my_interface_maximum_refinement_level = numLevels;

  setup_refinement_marker();

  if (numLevels > 0)
    set_global_ids_are_NOT_parallel_consistent();
}

void
CDFEM_Support::setup_refinement_marker()
{
  /* %TRACE% */ Traceback trace__("krino::CDFEM_Support::activate_nonconformal_adaptivity(const int num_levels)"); /* %TRACE% */

  my_nonconformal_adapt_marker_name = "CDFEM_NONCONFORMAL_MARKER";

  my_aux_meta.register_field(my_nonconformal_adapt_marker_name, FieldType::INTEGER, stk::topology::ELEMENT_RANK, 1, 1, get_universal_part());
  my_aux_meta.register_field(my_nonconformal_adapt_marker_name, FieldType::INTEGER, stk::topology::NODE_RANK, 1, 1, get_universal_part());
}

int
CDFEM_Support::get_ls_index(const LevelSet * ls) const
{
  for ( unsigned ls_index = 0; ls_index < my_ls_fields.size(); ++ls_index )
  {
    if (ls == ls_field(ls_index).ptr) return ls_index;
  }
  ThrowRuntimeError("Invalid level set field index");
}

void
CDFEM_Support::activate_nonconformal_adapt_target_count(const uint64_t target_count)
{
  /* %TRACE% */ Traceback trace__("CDFEM_Support::activate_nonconformal_adapt_target_count(const uint64_t target_count)"); /* %TRACE% */

  my_nonconformal_adapt_target_element_count = target_count;
  my_nonconformal_adapt_indicator_name = "CDFEM_ADAPTIVITY_ERROR_INDICATOR";

  my_aux_meta.register_field(my_nonconformal_adapt_indicator_name,
      FieldType::REAL,
      stk::topology::ELEMENT_RANK,
      1,
      1,
      get_universal_part());
}

} // namespace krino
