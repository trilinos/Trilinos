// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef PANZER_MODEL_FACTORY_TEMPLATE_MANAGER_HPP
#define PANZER_MODEL_FACTORY_TEMPLATE_MANAGER_HPP

#include "Phalanx_TemplateManager.hpp"

#include "Panzer_Traits.hpp"
#include "Panzer_Base.hpp"
#include "Panzer_ClosureModel_Factory_Base.hpp"
#include "Panzer_ClosureModel_Factory.hpp"

#include "Sacado_mpl_placeholders.hpp"

namespace panzer {

  template<typename Traits>
  class ClosureModelFactory_TemplateManager : 
    public PHX::TemplateManager<typename Traits::EvalTypes,
				panzer::ClosureModelFactoryBase,
                                panzer::ClosureModelFactory<_> > {

  public:

    ClosureModelFactory_TemplateManager() {}

    ~ClosureModelFactory_TemplateManager() {}

  };

} 

#endif 
