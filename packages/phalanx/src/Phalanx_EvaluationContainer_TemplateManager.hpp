// @HEADER
// *****************************************************************************
//        Phalanx: A Partial Differential Equation Field Evaluation 
//       Kernel for Flexible Management of Complex Dependency Chains
//
// Copyright 2008 NTESS and the Phalanx contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef PHX_SCALAR_CONTAINER_TEMPLATE_MANAGER_HPP
#define PHX_SCALAR_CONTAINER_TEMPLATE_MANAGER_HPP

#include "Phalanx_TemplateManager.hpp"
#include "Phalanx_EvaluationContainer.hpp"

#include "Sacado_mpl_placeholders.hpp"
using namespace Sacado::mpl::placeholders;

namespace PHX {

  template<typename Traits>
  class EvaluationContainer_TemplateManager : 
    public PHX::TemplateManager<typename Traits::EvalTypes,
				PHX::EvaluationContainerBase<Traits>,
				PHX::EvaluationContainer<_,Traits> > {

  public:

    EvaluationContainer_TemplateManager() {}

    ~EvaluationContainer_TemplateManager() {}

  };

} 

#endif 
