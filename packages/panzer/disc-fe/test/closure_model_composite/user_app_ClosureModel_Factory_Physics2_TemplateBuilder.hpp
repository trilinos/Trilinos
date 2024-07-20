// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef USER_APP_CLOSURE_MODEL_FACTORY_PHYSICS_2_TEMPLATE_BUILDER_HPP
#define USER_APP_CLOSURE_MODEL_FACTORY_PHYSICS_2_TEMPLATE_BUILDER_HPP

#include <string>
#include "Sacado_mpl_apply.hpp"
#include "Teuchos_RCP.hpp"
#include "user_app_ClosureModel_Factory_Physics2.hpp"

namespace user_app {

  class MyModelFactory_Physics2_TemplateBuilder {

    bool m_throw_if_model_not_found;

  public:

    MyModelFactory_Physics2_TemplateBuilder(bool throw_if_model_not_found) :
      m_throw_if_model_not_found(throw_if_model_not_found) {}

    template <typename EvalT>
    Teuchos::RCP<panzer::ClosureModelFactoryBase> build() const {
      return Teuchos::rcp( static_cast<panzer::ClosureModelFactoryBase*>
			   (new user_app::MyModelFactory_Physics2<EvalT>(m_throw_if_model_not_found)) );
    }
    
  };
  
}

#endif 
