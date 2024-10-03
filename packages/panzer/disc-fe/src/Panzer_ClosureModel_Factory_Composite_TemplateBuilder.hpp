// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef PANZER_CLOSURE_MODEL_FACTORY_COMPOSITE_TEMPLATE_BUILDER_HPP
#define PANZER_CLOSURE_MODEL_FACTORY_COMPOSITE_TEMPLATE_BUILDER_HPP

#include <string>
#include "Sacado_mpl_apply.hpp"
#include "Teuchos_RCP.hpp"
#include "Panzer_ClosureModel_Factory_Composite.hpp"

namespace panzer {

  class ClosureModelFactoryComposite_TemplateBuilder {

    std::vector<Teuchos::RCP<panzer::ClosureModelFactory_TemplateManager<panzer::Traits> > > m_factories;

  public:

    template <typename EvalT>
    Teuchos::RCP<panzer::ClosureModelFactoryBase> build() const {
      return Teuchos::rcp( static_cast<panzer::ClosureModelFactoryBase*>
			   (new panzer::ClosureModelFactoryComposite<EvalT>(m_factories)) );
    }
    
    void addFactory(const Teuchos::RCP<panzer::ClosureModelFactory_TemplateManager<panzer::Traits> >& factory)
    {
      m_factories.push_back(factory);
    }

  };
  
}

#endif 
