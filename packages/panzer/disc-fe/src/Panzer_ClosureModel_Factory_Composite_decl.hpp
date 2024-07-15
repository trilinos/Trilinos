// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef PANZER_CLOSURE_MODEL_FACTORY_COMPOSITE_DECL_HPP
#define PANZER_CLOSURE_MODEL_FACTORY_COMPOSITE_DECL_HPP

#include "Panzer_ClosureModel_Factory.hpp"
#include "Panzer_ClosureModel_Factory_TemplateManager.hpp"

namespace panzer {

  template<typename EvalT>
  class ClosureModelFactoryComposite : public panzer::ClosureModelFactory<EvalT> {

  public:

    ClosureModelFactoryComposite(const std::vector<Teuchos::RCP<panzer::ClosureModelFactory_TemplateManager<panzer::Traits> > >& factories);

    Teuchos::RCP< std::vector< Teuchos::RCP<PHX::Evaluator<panzer::Traits> > > >
    buildClosureModels(const std::string& model_id,
		       const Teuchos::ParameterList& models,
		       const panzer::FieldLayoutLibrary& fl,
		       const Teuchos::RCP<panzer::IntegrationRule>& ir,
		       const Teuchos::ParameterList& default_params,
		       const Teuchos::ParameterList& user_data,
		       const Teuchos::RCP<panzer::GlobalData>& global_data,
		       PHX::FieldManager<panzer::Traits>& fm) const;

  private:

    std::vector<Teuchos::RCP<panzer::ClosureModelFactory_TemplateManager<panzer::Traits> > > m_factories;

  };

}

#endif
