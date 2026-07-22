// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef USER_APP_STK_CLOSURE_MODEL_FACTORY_TEMPLATE_BUILDER_HPP
#define USER_APP_STK_CLOSURE_MODEL_FACTORY_TEMPLATE_BUILDER_HPP

#include <string>
#include "Sacado_mpl_apply.hpp"
#include "Teuchos_RCP.hpp"
#include "Panzer_Base.hpp"
#include "user_app_STKClosureModel_Factory.hpp"

namespace user_app {

  class STKModelFactory_TemplateBuilder {

  public:
    
    template <typename EvalT>
    Teuchos::RCP<panzer::ClosureModelFactoryBase> build() const {
      return Teuchos::rcp( static_cast<panzer::ClosureModelFactoryBase*>
			   (new user_app::STKModelFactory<EvalT>) );
    }
    
  };
  
}

#endif 
