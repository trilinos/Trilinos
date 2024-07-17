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

// Forward declaration
namespace Kokkos {
  template<typename DataT,typename... Props> class View;
}

namespace PHX {

  // Both required and dependent fields are all unmanaged.  Covers all
  // combinations of data types (static/dynamic, const/nonconst).
  PHX_EVALUATOR_CLASS(EvalUnmanaged)
    Tag<double> tag_a;
    Tag<double> tag_b;
    Tag<double> tag_c;
    Tag<double> tag_d;
    PHX::View<double**> a; // static evaluated
    PHX::View<const double**> b; // static dependent
    PHX::View<double**> c; // dynamic evalauted
    PHX::View<const double**> d; // dynamic dependent
  PHX_EVALUATOR_CLASS_END

  // Dummy to satisfy dependent unmanaged fields
  PHX_EVALUATOR_CLASS(EvalDummy)
    Tag<double> tag_b;
    Tag<double> tag_d;
    PHX::View<double**> b;
    PHX::View<double**> d;
  PHX_EVALUATOR_CLASS_END

}

#include "View_TestEvaluators_Def.hpp"

#endif
