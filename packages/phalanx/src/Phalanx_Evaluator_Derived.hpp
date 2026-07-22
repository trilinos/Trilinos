// @HEADER
// *****************************************************************************
//        Phalanx: A Partial Differential Equation Field Evaluation 
//       Kernel for Flexible Management of Complex Dependency Chains
//
// Copyright 2008 NTESS and the Phalanx contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef PHX_FIELD_EVALUATOR_DERIVED_HPP
#define PHX_FIELD_EVALUATOR_DERIVED_HPP

#include <vector>

#include "Phalanx_Evaluator_Base.hpp"
#include "Phalanx_Evaluator_Utilities.hpp"

namespace PHX {

  template<typename EvalT, typename Traits>
  class EvaluatorDerived : public PHX::EvaluatorBase<Traits> {
    
  public:
    
    EvaluatorDerived() {}

    virtual ~EvaluatorDerived() {}
    
  protected:
    
    PHX::EvaluatorUtilities<EvalT,Traits> utils;

  };

}

#endif
