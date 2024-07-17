// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef __Panzer_ResponseEvaluatorFactory_Functional_impl_hpp__
#define __Panzer_ResponseEvaluatorFactory_Functional_impl_hpp__

#include <string>

#include "PanzerDiscFE_config.hpp"

#include "Panzer_IntegrationRule.hpp"
#include "Panzer_PhysicsBlock.hpp"
#include "Panzer_Integrator_Scalar.hpp"
#include "Panzer_ResponseScatterEvaluator_Functional.hpp"
#include "Panzer_Response_Functional.hpp"
#include "Panzer_BlockedDOFManager.hpp"
#include "Panzer_GlobalIndexer_Utilities.hpp"

namespace panzer {

template <typename EvalT,typename LO,typename GO>
Teuchos::RCP<ResponseBase> ResponseEvaluatorFactory_Functional<EvalT,LO,GO>::
buildResponseObject(const std::string & responseName) const
{ 
  Teuchos::RCP<ResponseBase> response = Teuchos::rcp(new Response_Functional<EvalT>(responseName,comm_,linearObjFactory_)); 
  response->setRequiresDirichletAdjustment(applyDirichletToDerivative_);
 
  return response;
}

template <typename EvalT,typename LO,typename GO>
void ResponseEvaluatorFactory_Functional<EvalT,LO,GO>::
buildAndRegisterEvaluators(const std::string & responseName,
                           PHX::FieldManager<panzer::Traits> & fm,
                           const panzer::PhysicsBlock & physicsBlock,
                           const Teuchos::ParameterList & /* user_data */) const
{
   using Teuchos::RCP;
   using Teuchos::rcp;


   // build integration evaluator (integrate over element)
   if(requiresCellIntegral_) {
     std::string field = (quadPointField_=="" ? responseName : quadPointField_);

     // build integration rule to use in cell integral
     RCP<IntegrationRule> ir = rcp(new IntegrationRule(cubatureDegree_,physicsBlock.cellData()));

     Teuchos::ParameterList pl;
     pl.set("Integral Name", field);
     pl.set("Integrand Name",field);
     pl.set("IR",ir);

     Teuchos::RCP<PHX::Evaluator<panzer::Traits> > eval 
         = Teuchos::rcp(new Integrator_Scalar<EvalT,panzer::Traits>(pl));
 
     this->template registerEvaluator<EvalT>(fm, eval);
   }


   // build scatter evaluator
   {
     Teuchos::RCP<FunctionalScatterBase> scatterObj;
     if(linearObjFactory_!=Teuchos::null) {

        TEUCHOS_ASSERT(linearObjFactory_->getDomainGlobalIndexer()!=Teuchos::null);

        auto ugi = Teuchos::rcp_dynamic_cast<const GlobalIndexer>(linearObjFactory_->getDomainGlobalIndexer());
        auto bugi = Teuchos::rcp_dynamic_cast<const BlockedDOFManager>(linearObjFactory_->getDomainGlobalIndexer());

        if(ugi!=Teuchos::null) {
          std::vector<Teuchos::RCP<const GlobalIndexer> > ugis; 
          ugis.push_back(ugi);

          scatterObj = Teuchos::rcp(new FunctionalScatter<LO,GO>(ugis));
        }
        else if(bugi!=Teuchos::null) {
          scatterObj = Teuchos::rcp(new FunctionalScatter<LO,GO>(nc2c_vector(bugi->getFieldDOFManagers())));
        }
        else {
          TEUCHOS_ASSERT(false); // no real global indexer to use
        }
     }

     std::string field = (quadPointField_=="" ? responseName : quadPointField_);

     // build useful evaluator
     Teuchos::RCP<PHX::Evaluator<panzer::Traits> > eval 
         = Teuchos::rcp(new ResponseScatterEvaluator_Functional<EvalT,panzer::Traits>(field,responseName,physicsBlock.cellData(),scatterObj));

     this->template registerEvaluator<EvalT>(fm, eval);

     // require last field
     fm.template requireField<EvalT>(*eval->evaluatedFields()[0]);
   }
}

template <typename EvalT,typename LO,typename GO>
bool ResponseEvaluatorFactory_Functional<EvalT,LO,GO>::
typeSupported() const
{
  if(   PHX::print<EvalT>()==PHX::print<panzer::Traits::Residual>() ||
        PHX::print<EvalT>()==PHX::print<panzer::Traits::Tangent>()
    )
    return true;

  if(PHX::print<EvalT>()==PHX::print<panzer::Traits::Jacobian>())
    return linearObjFactory_!=Teuchos::null;

#ifdef Panzer_BUILD_HESSIAN_SUPPORT
  if(PHX::print<EvalT>()==PHX::print<panzer::Traits::Hessian>()) {
    return linearObjFactory_!=Teuchos::null;
  }
#endif

  return false;
}

}

#endif
