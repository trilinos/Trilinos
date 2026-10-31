// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include <Akri_LevelSet_Phase_Utils.hpp>

#include <limits>
#include <ostream>
#include <sstream>
#include <stdexcept>

#include <stk_mesh/base/MetaData.hpp>
#include <stk_util/diag/Timer.hpp>
#include <stk_util/util/ReportHandler.hpp>

#include <Akri_AuxMetaData.hpp>
#include <Akri_DiagWriter.hpp>
#include <Akri_LevelSet.hpp>
#include <Akri_Phase_Support.hpp>
#include <Akri_Surface_Manager.hpp>

namespace krino {

//--------------------------------------------------------------------------------
CDFEM_Irreversible_Phase_Support *
CDFEM_Irreversible_Phase_Support::get_if_present(const stk::mesh::MetaData & meta)
{
  CDFEM_Irreversible_Phase_Support * mgr = const_cast<CDFEM_Irreversible_Phase_Support *>(meta.get_attribute<CDFEM_Irreversible_Phase_Support>());
  return mgr;
}
//--------------------------------------------------------------------------------
CDFEM_Irreversible_Phase_Support &
CDFEM_Irreversible_Phase_Support::get(const stk::mesh::MetaData & meta)
{
  CDFEM_Irreversible_Phase_Support * support = get_if_present(meta);
  STK_ThrowRequireMsg(nullptr != support, "Could not find CDFEM_Irreversible_Phase_Support attribute on MetaData.");
  return *support;
}
//--------------------------------------------------------------------------------
CDFEM_Irreversible_Phase_Support &
CDFEM_Irreversible_Phase_Support::get(stk::mesh::MetaData & meta)
{
  CDFEM_Irreversible_Phase_Support * support = get_if_present(meta);
  if (nullptr == support)
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
  STK_ThrowInvalidArgMsgIf(has_death && has_irreversible_phase_change,
    "Cannot have both CDFEM death and CDFEM irreversible phase change in the same problem.");

  for (auto&& death_spec : my_death_specs)
  {
    STK_ThrowInvalidArgMsgIf(death_spec.name() == death_name,
      "Only one CDFEM death specification with the same name is allowed per mesh (" << death_name << "). ");
  }
  my_death_specs.push_back(CDFEM_Inequality_Spec(death_name));
  return &my_death_specs.back();
}
//--------------------------------------------------------------------------
const CDFEM_Inequality_Spec * CDFEM_Irreversible_Phase_Support::get_death_spec_for_ls(const LevelSet * ls) const
{
  for (auto && deathSpec : get_death_specs())
    if (&deathSpec.get_levelset() == ls)
      return &deathSpec;
  return nullptr;
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
    STK_ThrowErrorMsg("Bad integer passed.  Could not convert passed in integer "
                      << "death_criterion_compare = " << inequality_criterion
                      << " into a CDFEM_Death_Spec::DeathCriterionCompareType enum.  "
                      << "The code and the Sierra XML database may be out of sync." << std::endl
                      << StackTrace);
  }
  return INEQUALITY_CRITERION_TYPE_UNDEFINED; // unreachable; STK_ThrowErrorMsg does not return
}

void CDFEM_Inequality_Spec::add_element_block(stk::mesh::Part * block)
{
  // Note: It is ok to add multiple volume names to vector, thus no return status given
  myElementBlocks.push_back(block);
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

void CDFEM_Inequality_Spec::create_levelset(stk::mesh::MetaData & meta, stk::diag::Timer & parent_timer)
{
  /* %TRACE[ON]% */ Trace trace__("krino::CDFEM_Death_Spec::create_levelset()"); /* %TRACE% */

  STK_ThrowInvalidArgMsgIf(
      Surface_Manager::get(meta).has_levelset(my_name),
      "Region already has a LevelSet named " << my_name);

  my_ls = &LevelSet::build(meta, my_name, parent_timer);
  my_ls->set_tracked_isovar(get_threshold_variable_name(), get_threshold_value());

  int deactivated_phase_sign = ( get_criterion_compare_type() == CDFEM_Inequality_Spec::LESS_THAN ) ? -1 : 1;
  my_deactivated_phase.add(my_ls->get_identifier(), deactivated_phase_sign);
  my_active_phase.add(my_ls->get_identifier(), -1 * deactivated_phase_sign);
}

std::vector<LS_Field> get_levelset_fields(const stk::mesh::MetaData & meta)
{
  const Surface_Manager & surfaceManager = Surface_Manager::get(meta);
  const CDFEM_Irreversible_Phase_Support * deathManager = CDFEM_Irreversible_Phase_Support::get_if_present(meta);

  std::vector<LS_Field> lsFields;
  for (auto&& ls : surfaceManager.get_levelsets())
  {
    const CDFEM_Inequality_Spec * deathSpec = deathManager ? deathManager->get_death_spec_for_ls(ls.get()) : nullptr;
    LS_Field LSField(ls->name(),ls->get_identifier(),ls->get_isovar_field(),ls->get_isoval(),ls.get(), deathSpec);
    lsFields.push_back(LSField);
  }
  return lsFields;
}

static bool field_exists_on_part(const FieldRef field, const stk::mesh::Part & part)
{
  stk::mesh::Selector field_selector = AuxMetaData::get(part.mesh_meta_data()).selectField(field, stk::topology::ELEMENT_RANK);
  return field_selector(part);
}

static std::ostream & dump_fields_on_part(std::ostream & os, const stk::mesh::Part & part)
{
  const auto & meta = part.mesh_meta_data();
  for(auto && field : meta.get_fields(stk::topology::NODE_RANK))
  {
    if(field_exists_on_part(*field, part))
      os << "  * " << field->name() << "\n";
  }
  return os;
}

void check_isovariable_field_existence_on_decomposed_blocks(const stk::mesh::MetaData & meta, const std::vector<LS_Field> & lsFields, const bool conformal_parts_require_field)
{
  const Phase_Support & phaseSupport = Phase_Support::get(meta);

  for(auto * part : meta.get_mesh_parts())
  {
    if(part->primary_entity_rank() != stk::topology::ELEMENT_RANK ||
        phaseSupport.is_nonconformal(*part) ||
        !phaseSupport.is_decomposed(*part))
    {
      continue;
    }

    const stk::mesh::Part & nonconformingPart = phaseSupport.get_nonconforming_part(*part);

    /* 3 cases we care about here:
     * 1) Initial IO part (i.e. block_1) - Must have field
     * 2) Active conformal part (either LS or active part for irreversible) - Must have field
     * 3) Deactivated conformal part (irreversible only) - Skip
     */

    // If part is the initial IO part the phase will be empty
    const auto & phase = phaseSupport.get_iopart_phase(*part);

    // Steady state problems that only do a single decomposition can
    // get away with only defining the isovariable on the initial IO part
    // so skip checking for the field in that case.
    if(!conformal_parts_require_field && !phase.empty()) continue;

    for (auto&& lsField : lsFields)
    {
      const CDFEM_Inequality_Spec * deathSpec = lsField.deathPtr;

      const bool is_dead_block = deathSpec && phase.contain(deathSpec->get_deactivated_phase());
      const bool LS_uses_part =
          phaseSupport.level_set_is_used_by_nonconforming_part(lsField.identifier, nonconformingPart);
      if(LS_uses_part && !is_dead_block)
      {
        const auto field = lsField.isovar;
        if(!field_exists_on_part(field, *part))
        {
          std::ostringstream err_msg;
          err_msg << "CDFEM isovariable field " << field.name() << " is not defined on "
              << part->name() << " that is supposed to be decomposed.\n"
              << "Available fields on " << part->name() << " are:\n";
          dump_fields_on_part(err_msg, *part);
          throw std::runtime_error(err_msg.str());
        }
      }
    }
  }
}

void check_isovariable_field_existence_on_decomposed_blocks(const stk::mesh::MetaData & meta, const bool conformal_parts_require_field)
{
  check_isovariable_field_existence_on_decomposed_blocks(meta, get_levelset_fields(meta), conformal_parts_require_field);
}

} // namespace krino
