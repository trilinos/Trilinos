// @HEADER
// *****************************************************************************
//        Phalanx: A Partial Differential Equation Field Evaluation 
//       Kernel for Flexible Management of Complex Dependency Chains
//
// Copyright 2008 NTESS and the Phalanx contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef PHX_FIELD_TEST_EVALUATORS_HPP
#define PHX_FIELD_TEST_EVALUATORS_HPP

#include "Phalanx_Evaluator_Macros.hpp"
#include "Phalanx_Field.hpp"

namespace PHX {

  // Both required and dependent fields are all unmanaged.  Covers all
  // combinations of data types (static/dynamic, const/nonconst).
  PHX_EVALUATOR_CLASS(EvalUnmanaged)
  PHX::Field<double,2> a; // static evaluated
  PHX::Field<const double,2> b; // static dependent
  PHX::Field<double,2> c; // dynamic evalauted
  PHX::Field<const double,2> d; // dynamic dependent
  PHX_EVALUATOR_CLASS_END

  // Dummy to satisfy dependent unmanaged fields
  PHX_EVALUATOR_CLASS(EvalDummy)
  PHX::Field<double,2> b;
  PHX::Field<double,2> d;
  PHX_EVALUATOR_CLASS_END

}

#include "Field_TestEvaluators_Def.hpp"

#endif
