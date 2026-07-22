// @HEADER
// *****************************************************************************
//        Phalanx: A Partial Differential Equation Field Evaluation 
//       Kernel for Flexible Management of Complex Dependency Chains
//
// Copyright 2008 NTESS and the Phalanx contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef PHX_MDFIELD_TEST_EVALUATORS_HPP
#define PHX_MDFIELD_TEST_EVALUATORS_HPP

#include "Phalanx_Evaluator_Macros.hpp"
#include "Phalanx_MDField.hpp"

namespace PHX {

  // Both required and dependent fields are all unmanaged.  Covers all
  // combinations of data types (static/dynamic, const/nonconst).
  PHX_EVALUATOR_CLASS(EvalUnmanaged)
  PHX::MDField<double,CELL,BASIS> a; // static evaluated
  PHX::MDField<const double,CELL,BASIS> b; // static dependent
  PHX::MDField<double> c; // dynamic evalauted
  PHX::MDField<const double> d; // dynamic dependent
  PHX_EVALUATOR_CLASS_END

  // Dummy to satisfy dependent unmanaged fields
  PHX_EVALUATOR_CLASS(EvalDummy)
  PHX::MDField<double,CELL,BASIS> b;
  PHX::MDField<double> d;
  PHX_EVALUATOR_CLASS_END

}

#include "MDField_TestEvaluators_Def.hpp"

#endif
