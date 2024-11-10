// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef USER_APP_CLOSURE_MODEL_FACTORY_TEMPLATE_BUILDER_HPP
#define USER_APP_CLOSURE_MODEL_FACTORY_TEMPLATE_BUILDER_HPP

#include <string>
#include "Sacado_mpl_apply.hpp"
#include "Teuchos_RCP.hpp"
#include "Panzer_Base.hpp"
#include "user_app_ClosureModel_Factory.hpp"

namespace user_app {

  class MyModelFactory_TemplateBuilder {
    Teuchos::RCP<const panzer::LinearObjFactory<panzer::Traits> > distr_param_lof;

  public:

    void setDistributedParameterLOF(const Teuchos::RCP<const panzer::LinearObjFactory<panzer::Traits> > & dpl)
    { distr_param_lof = dpl; }
    
    template <typename EvalT>
    Teuchos::RCP<panzer::ClosureModelFactoryBase> build() const {
      Teuchos::RCP<user_app::MyModelFactory<EvalT> > closure_factory =
          Teuchos::rcp(new user_app::MyModelFactory<EvalT>) ;
      closure_factory->setDistributedParameterLOF(distr_param_lof);

      return Teuchos::rcp_static_cast<panzer::ClosureModelFactoryBase>(closure_factory);
    }
    
  };
  
}

#endif 
