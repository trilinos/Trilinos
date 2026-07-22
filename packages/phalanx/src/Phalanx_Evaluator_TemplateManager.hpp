// @HEADER
// *****************************************************************************
//        Phalanx: A Partial Differential Equation Field Evaluation 
//       Kernel for Flexible Management of Complex Dependency Chains
//
// Copyright 2008 NTESS and the Phalanx contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef PHX_FIELDEVALUATOR_TEMPLATE_MANAGER_HPP
#define PHX_FIELDEVALUATOR_TEMPLATE_MANAGER_HPP

#include "Phalanx_TemplateManager.hpp"
#include "Phalanx_Evaluator_Derived.hpp"

#include "Sacado_mpl_placeholders.hpp"
using namespace Sacado::mpl::placeholders;

namespace PHX {
  
  template<typename Traits>
  class Evaluator_TemplateManager : 
    public PHX::TemplateManager<typename Traits::EvalTypes,
				PHX::EvaluatorBase<Traits>,
				PHX::EvaluatorDerived<_,Traits> > {
    
  public:
    
    Evaluator_TemplateManager() {}
    
    ~Evaluator_TemplateManager() {}
    
  };

} 

#endif 
