// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef USER_APP_CLOSURE_MODEL_FACTORY_HPP
#define USER_APP_CLOSURE_MODEL_FACTORY_HPP

#include "Panzer_ClosureModel_Factory.hpp"

namespace panzer {
  class InputEquationSet;

  template <typename> class LinearObjFactory;
}

namespace user_app {

  template<typename EvalT>
  class MyModelFactory : public panzer::ClosureModelFactory<EvalT> {
    Teuchos::RCP<const panzer::LinearObjFactory<panzer::Traits> > distr_param_lof;

  public:
    void setDistributedParameterLOF(const Teuchos::RCP<const panzer::LinearObjFactory<panzer::Traits> > & dpl)
    { distr_param_lof = dpl; }

    Teuchos::RCP< std::vector< Teuchos::RCP<PHX::Evaluator<panzer::Traits> > > >
    buildClosureModels(const std::string& model_id,
		       const Teuchos::ParameterList& models,
		       const panzer::FieldLayoutLibrary& fl,
		       const Teuchos::RCP<panzer::IntegrationRule>& ir,
		       const Teuchos::ParameterList& default_params,
		       const Teuchos::ParameterList& user_data,
		       const Teuchos::RCP<panzer::GlobalData>& global_data,
		       PHX::FieldManager<panzer::Traits>& fm) const;

  };

}

#include "user_app_ClosureModel_Factory_impl.hpp"

#endif
