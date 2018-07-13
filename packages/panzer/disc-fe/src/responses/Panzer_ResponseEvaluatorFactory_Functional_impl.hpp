// @HEADER
// ***********************************************************************
//
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//                 Copyright (2011) Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Roger P. Pawlowski (rppawlo@sandia.gov) and
// Eric C. Cyr (eccyr@sandia.gov)
// ***********************************************************************
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
#include "Panzer_UniqueGlobalIndexer_Utilities.hpp"

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

        auto ugi = Teuchos::rcp_dynamic_cast<const UniqueGlobalIndexer<LO,GO> >(linearObjFactory_->getDomainGlobalIndexer());
        auto bugi = Teuchos::rcp_dynamic_cast<const BlockedDOFManager<LO,GO> >(linearObjFactory_->getDomainGlobalIndexer());

        if(ugi!=Teuchos::null) {
          std::vector<Teuchos::RCP<const UniqueGlobalIndexer<LO,GO> > > ugis; 
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
  if(   PHX::typeAsString<EvalT>()==PHX::typeAsString<panzer::Traits::Residual>() ||
        PHX::typeAsString<EvalT>()==PHX::typeAsString<panzer::Traits::Tangent>()
    )
    return true;

  if(PHX::typeAsString<EvalT>()==PHX::typeAsString<panzer::Traits::Jacobian>())
    return linearObjFactory_!=Teuchos::null;

#ifdef Panzer_BUILD_HESSIAN_SUPPORT
  if(PHX::typeAsString<EvalT>()==PHX::typeAsString<panzer::Traits::Hessian>()) {
    return linearObjFactory_!=Teuchos::null;
  }
#endif

  return false;
}

}

#endif
