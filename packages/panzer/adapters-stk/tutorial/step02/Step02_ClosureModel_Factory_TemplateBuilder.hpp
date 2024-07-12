// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef __Step02_ClosureModel_Factory_TemplateBuilder_hpp__
#define __Step02_ClosureModel_Factory_TemplateBuilder_hpp__

#include <string>
#include "Sacado_mpl_apply.hpp"
#include "Teuchos_RCP.hpp"

#include "Step02_ClosureModel_Factory.hpp"

namespace user_app {

class ClosureModelFactory_TemplateBuilder {
  Teuchos::RCP<const panzer::LinearObjFactory<panzer::Traits> > distr_param_lof;

public:

  template <typename EvalT>
  Teuchos::RCP<panzer::ClosureModelFactoryBase> build() const 
  {
    using Teuchos::RCP;
    using Teuchos::rcp;
 
    RCP<user_app::ClosureModelFactory<EvalT> > closure_factory = rcp(new user_app::ClosureModelFactory<EvalT>);

    return Teuchos::rcp_static_cast<panzer::ClosureModelFactoryBase>(closure_factory);
  }
    
};
  
}

#endif 
