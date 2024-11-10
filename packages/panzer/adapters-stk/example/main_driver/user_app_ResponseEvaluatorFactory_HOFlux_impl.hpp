// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef __user_app_ResponseEvaluatorFactory_HOFlux_impl_hpp__
#define __user_app_ResponseEvaluatorFactory_HOFlux_impl_hpp__

#include "Panzer_Normals.hpp"
#include "Panzer_Sum.hpp"
#include "Panzer_DotProduct.hpp"
#include "Panzer_Integrator_Scalar.hpp"
#include "Panzer_SubcellSum.hpp"

namespace user_app {

template <typename EvalT,typename LO,typename GO>
void ResponseEvaluatorFactory_HOFlux<EvalT,LO,GO>::
buildAndRegisterEvaluators(const std::string & responseName,
                           PHX::FieldManager<panzer::Traits> & fm,
                           const panzer::PhysicsBlock & physicsBlock,
                           const Teuchos::ParameterList & user_data) const
{
  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::ParameterList; 

  // create basis 
  RCP<const panzer::FieldLibrary> fieldLib = physicsBlock.getFieldLibrary();
  RCP<const panzer::PureBasis> basis = fieldLib->lookupBasis("TEMPERATURE");

  {
    Teuchos::ParameterList pl;
    pl.set("Sum Name",responseName);
    pl.set("Field Name","RESIDUAL_TEMPERATURE");
    pl.set("Basis",basis);
    pl.set("Multiplier",1.0);
    
    Teuchos::RCP<PHX::Evaluator<panzer::Traits> > eval 
        = Teuchos::rcp(new panzer::SubcellSum<EvalT,panzer::Traits>(pl));
 
    this->template registerEvaluator<EvalT>(fm, eval);
  }

  panzer::ResponseEvaluatorFactory_Functional<EvalT,LO,GO>::buildAndRegisterEvaluators(responseName,fm,physicsBlock,user_data);
}

}

#endif
