// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef PANZER_RESPONSE_EVALUATOR_FACTORY_TEMPLATE_MANAGER_H
#define PANZER_RESPONSE_EVALUATOR_FACTORY_TEMPLATE_MANAGER_H

#include "Phalanx_TemplateManager.hpp"

#include "Panzer_Traits.hpp"
#include "Panzer_Base.hpp"

#include "Panzer_ResponseEvaluatorFactory.hpp"

#include "Sacado_mpl_placeholders.hpp"

namespace panzer {

  template<typename Traits>
  class ResponseEvaluatorFactory_TemplateManager : 
    public PHX::TemplateManager<typename Traits::EvalTypes,
				panzer::ResponseEvaluatorFactoryBase,
				panzer::ResponseEvaluatorFactory<_> > {

  public:

    ResponseEvaluatorFactory_TemplateManager() {}

    ~ResponseEvaluatorFactory_TemplateManager() {}

  };

} 

#endif 
