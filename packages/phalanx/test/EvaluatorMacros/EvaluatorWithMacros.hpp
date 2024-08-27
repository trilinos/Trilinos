// @HEADER
// *****************************************************************************
//        Phalanx: A Partial Differential Equation Field Evaluation 
//       Kernel for Flexible Management of Complex Dependency Chains
//
// Copyright 2008 NTESS and the Phalanx contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef PHX_EVALUATOR_WITH_MACROS_HPP
#define PHX_EVALUATOR_WITH_MACROS_HPP

#include "Phalanx_Evaluator_Macros.hpp"

namespace PHX {

  // Macro with no pre/post evaluate methods
  PHX_EVALUATOR_CLASS(EvaluatorWithMacros1)
  public:
    void evaluates(const std::string& field_name);
    void depends(const std::string& field_name);
    void bindField(const PHX::FieldTag& ft, const std::any& f);
  PHX_EVALUATOR_CLASS_END
  
  // Macro with no pre/post evaluate methods
  PHX_EVALUATOR_CLASS_PP(EvaluatorWithMacros2)
  public:
    void evaluates(const std::string& field_name);
    void depends(const std::string& field_name);
    void bindField(const PHX::FieldTag& ft, const std::any& f);
  PHX_EVALUATOR_CLASS_END

}

#include "EvaluatorWithMacros_Def.hpp"

#endif
