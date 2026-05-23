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
     // If blocked, pull out the local indexer for the particular field
     auto tmp_global_indexer = globalIndexer_;
     auto blocked_dof_manager = Teuchos::rcp_dynamic_cast<const panzer::BlockedDOFManager>(globalIndexer_);
     if (nonnull(blocked_dof_manager)) {
       int field_number = blocked_dof_manager->getFieldNum(this->fieldName_);
       int product_vector_block_index = blocked_dof_manager->getFieldBlock(field_number);
       tmp_global_indexer = (blocked_dof_manager->getFieldDOFManagers()[product_vector_block_index]);
     }

     Teuchos::RCP<ProbeScatterBase> scatterObj =
         (tmp_global_indexer != Teuchos::null) ?  Teuchos::rcp(new ProbeScatter<LO,GO>(tmp_global_indexer)) : Teuchos::null;
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
  // TODO BWR does this need to happen??
  if(   PHX::print<EvalT>()==PHX::print<panzer::Traits::Residual>()// ||
        //PHX::print<EvalT>()==PHX::print<panzer::Traits::Tangent>()
    )
    return true;

  if(PHX::print<EvalT>()==PHX::print<panzer::Traits::Jacobian>())
    return linearObjFactory_!=Teuchos::null;

  if(PHX::print<EvalT>()==PHX::print<panzer::Traits::Tangent>())
    return linearObjFactory_!=Teuchos::null;

#ifdef Panzer_BUILD_HESSIAN_SUPPORT
  if(PHX::print<EvalT>()==PHX::print<panzer::Traits::Hessian>()) {
    return linearObjFactory_!=Teuchos::null;
  }
#endif

  return false;

  // TODO BWR REMOVE ME IF WE KEEP ABOVE
//  if(PHX::print<EvalT>()==PHX::print<panzer::Traits::Residual>()  ||
//     PHX::print<EvalT>()==PHX::print<panzer::Traits::Tangent>() ||
//     PHX::print<EvalT>()==PHX::print<panzer::Traits::Jacobian>()
//    )
//    return true;
//
//  return false;
}

}

#endif
