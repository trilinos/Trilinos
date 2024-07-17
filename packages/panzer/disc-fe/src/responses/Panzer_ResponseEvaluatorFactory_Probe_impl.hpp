// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef __Panzer_ResponseEvaluatorFactory_Probe_impl_hpp__
#define __Panzer_ResponseEvaluatorFactory_Probe_impl_hpp__

#include <string>

#include "PanzerDiscFE_config.hpp"

#include "Panzer_IntegrationRule.hpp"
#include "Panzer_PhysicsBlock.hpp"
#include "Panzer_CellExtreme.hpp"
#include "Panzer_ResponseScatterEvaluator_Probe.hpp"
#include "Panzer_Response_Probe.hpp"

namespace panzer {

template <typename EvalT,typename LO,typename GO>
Teuchos::RCP<ResponseBase> ResponseEvaluatorFactory_Probe<EvalT,LO,GO>::
buildResponseObject(const std::string & responseName) const
{
  Teuchos::RCP<ResponseBase> response = Teuchos::rcp(new Response_Probe<EvalT>(responseName,comm_,linearObjFactory_));
  response->setRequiresDirichletAdjustment(applyDirichletToDerivative_);

  return response;
}

template <typename EvalT,typename LO,typename GO>
void ResponseEvaluatorFactory_Probe<EvalT,LO,GO>::
buildAndRegisterEvaluators(const std::string & responseName,
                           PHX::FieldManager<panzer::Traits> & fm,
                           const panzer::PhysicsBlock & physicsBlock,
                           const Teuchos::ParameterList & user_data) const
{
   using Teuchos::RCP;
   using Teuchos::rcp;

   // build scatter evaluator
   {
     Teuchos::RCP<ProbeScatterBase> scatterObj =
         (globalIndexer_!=Teuchos::null) ?  Teuchos::rcp(new ProbeScatter<LO,GO>(globalIndexer_)) : Teuchos::null;
     std::string field = (fieldName_=="" ? responseName : fieldName_);

     // Get basis and integration rule associated with field
     std::vector<panzer::StrPureBasisPair> blockFields = physicsBlock.getProvidedDOFs();
     RCP<const panzer::PureBasis> basis;
     for (auto&& v : blockFields) {
       if (v.first == field) {
         basis = v.second;
         break;
       }
     }
     RCP<panzer::IntegrationRule> ir = physicsBlock.getIntegrationRules().at(cubatureDegree_);

     // build useful evaluator
     Teuchos::RCP<PHX::Evaluator<panzer::Traits> > eval
       = Teuchos::rcp(new ResponseScatterEvaluator_Probe<EvalT,panzer::Traits,LO,GO>(responseName,
                                                                                     field,
                                                                                     fieldComponent_,
                                                                                     point_,
                                                                                     *ir,
                                                                                     basis,
                                                                                     globalIndexer_,
                                                                                     scatterObj));

     this->template registerEvaluator<EvalT>(fm, eval);

     // require last field
     fm.template requireField<EvalT>(*eval->evaluatedFields()[0]);
   }
}

template <typename EvalT,typename LO,typename GO>
bool ResponseEvaluatorFactory_Probe<EvalT,LO,GO>::
typeSupported() const
{
  if(PHX::print<EvalT>()==PHX::print<panzer::Traits::Residual>()  ||
     PHX::print<EvalT>()==PHX::print<panzer::Traits::Tangent>() ||
     PHX::print<EvalT>()==PHX::print<panzer::Traits::Jacobian>()
    )
    return true;

  return false;
}

}

#endif
