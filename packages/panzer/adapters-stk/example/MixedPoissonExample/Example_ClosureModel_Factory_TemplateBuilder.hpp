// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef __Example_ClosureModel_Factory_TemplateBuilder_hpp__
#define __Example_ClosureModel_Factory_TemplateBuilder_hpp__

#include <string>
#include "Sacado_mpl_apply.hpp"
#include "Teuchos_RCP.hpp"
#include "Example_ClosureModel_Factory.hpp"

namespace Example {

class ClosureModelFactory_TemplateBuilder {
public:
    
   template <typename EvalT>
   Teuchos::RCP<panzer::ClosureModelFactoryBase> build() const 
   {
      return Teuchos::rcp( static_cast<panzer::ClosureModelFactoryBase*>(new Example::ModelFactory<EvalT>) );
   }
    
};
  
}

#endif 
