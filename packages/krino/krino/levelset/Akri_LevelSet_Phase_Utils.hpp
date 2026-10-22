// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#ifndef KRINO_KRINO_LEVELSET_AKRI_LEVELSET_PHASE_UTILS_HPP_
#define KRINO_KRINO_LEVELSET_AKRI_LEVELSET_PHASE_UTILS_HPP_

#include <stk_mesh/base/Part.hpp>
#include <stk_util/util/ReportHandler.hpp>

#include <Akri_LS_Field.hpp>
#include <Akri_PhaseTag.hpp>

#include <string>
#include <vector>

namespace stk { namespace mesh { class MetaData; } }
namespace stk { namespace diag { class Timer; } }

namespace krino {

class LevelSet;

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

  void add_element_block(stk::mesh::Part * block);
  bool set_threshold_variable_name(const std::string & a_threshold_variable_name);
  bool set_threshold_value(const double a_threshold_value);
  bool set_criterion_compare_type(const InequalityCriterionType & a_criterion_compare_type);

  const std::vector<stk::mesh::Part*> & get_element_blocks() const { return myElementBlocks; }
  const std::string & get_threshold_variable_name () const { return my_threshold_variable_name; }
  double get_threshold_value () const { return my_threshold_value; }
  InequalityCriterionType get_criterion_compare_type() const { return my_criterion_compare_type; }
  const PhaseTag & get_deactivated_phase() const { return my_deactivated_phase; }
  const PhaseTag & get_active_phase() const { return my_active_phase; }
  void create_levelset(stk::mesh::MetaData & meta, stk::diag::Timer & parent_timer);
  const LevelSet & get_levelset() const { STK_ThrowRequire(my_ls != NULL); return *my_ls; }
  LevelSet & get_levelset() { STK_ThrowRequire(my_ls != NULL); return *my_ls; }

  // for unit testing
  void set_phases(const PhaseTag & active_phase, const PhaseTag & inactive_phase)
  {
    my_active_phase = active_phase;
    my_deactivated_phase = inactive_phase;
  }
protected:
  std::string my_name;
  std::vector<stk::mesh::Part *> myElementBlocks;
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
  bool decompose_at_start_of_time_step() { STK_ThrowAssert(has_irreversible_phase_change == !has_death); return has_irreversible_phase_change; }

private:
  CDFEM_Irreversible_Phase_Support() : has_death(false), has_irreversible_phase_change(false) {}
private:
  bool has_death;
  bool has_irreversible_phase_change;
  CDFEM_Inequality_Spec_Vec my_death_specs;
};

std::vector<LS_Field> get_levelset_fields(const stk::mesh::MetaData & meta);
void check_isovariable_field_existence_on_decomposed_blocks(const stk::mesh::MetaData & meta, const std::vector<LS_Field> & lsFields, const bool conformal_parts_require_field);
void check_isovariable_field_existence_on_decomposed_blocks(const stk::mesh::MetaData & meta, const bool conformal_parts_require_field);

} // namespace krino

#endif // KRINO_KRINO_LEVELSET_AKRI_LEVELSET_PHASE_UTILS_HPP_
