// @HEADER
// *****************************************************************************
//        Phalanx: A Partial Differential Equation Field Evaluation 
//       Kernel for Flexible Management of Complex Dependency Chains
//
// Copyright 2008 NTESS and the Phalanx contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef PHX_FIELDEVALUATOR_TEMPLATE_BUILDER_HPP
#define PHX_FIELDEVALUATOR_TEMPLATE_BUILDER_HPP

#include <string>
#include "Sacado_mpl_apply.hpp"
#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"

namespace PHX {

  template<typename Traits, typename ObjectT>
  class Evaluator_TemplateBuilder {
    
    Teuchos::RCP<Teuchos::ParameterList> p;

  public:
    
    Evaluator_TemplateBuilder(const Teuchos::RCP<Teuchos::ParameterList>& param) :
      p(param) {}

    template <typename ScalarT>
    Teuchos::RCP< PHX::EvaluatorBase<Traits> > build() const {
      typedef typename Sacado::mpl::apply<ObjectT,ScalarT>::type type;
      return Teuchos::rcp( static_cast< PHX::EvaluatorBase<Traits>* > (new type(*p)) );
    }
    
  };
  
}

#endif 
