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

#ifdef PHALANX_HIDE_DEPRECATED_CODE
  namespace detail {
    template<typename Traits>
    struct EvaluatorDerivedOf {
      template<typename EvalT> using type = PHX::EvaluatorDerived<EvalT,Traits>;
    };
  }
#endif
  
  template<typename Traits>
  class Evaluator_TemplateManager : 
#ifdef PHALANX_HIDE_DEPRECATED_CODE
    public PHX::TemplateManager<typename Traits::EvalTypes,
                                PHX::detail::EvaluatorDerivedOf<Traits>::template type,
                                PHX::EvaluatorBase<Traits> > {
#else
    public PHX::TemplateManager<typename Traits::EvalTypes,
				PHX::EvaluatorBase<Traits>,
				PHX::EvaluatorDerived<_,Traits> > {
#endif
    
  public:
    
    Evaluator_TemplateManager() {}
    
    ~Evaluator_TemplateManager() {}
    
  };

} 

#endif 
