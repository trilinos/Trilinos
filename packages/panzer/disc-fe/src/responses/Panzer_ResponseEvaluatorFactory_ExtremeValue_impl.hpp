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

#ifndef __Panzer_ResponseEvaluatorFactory_ExtremeValue_impl_hpp__
#define __Panzer_ResponseEvaluatorFactory_ExtremeValue_impl_hpp__

#include <string>

#include "PanzerDiscFE_config.hpp"

#include "Panzer_IntegrationRule.hpp"
#include "Panzer_PhysicsBlock.hpp"
#include "Panzer_CellExtreme.hpp"
#include "Panzer_ResponseScatterEvaluator_ExtremeValue.hpp"
#include "Panzer_Response_ExtremeValue.hpp"

namespace panzer {

template <typename EvalT,typename LO,typename GO>
Teuchos::RCP<ResponseBase> ResponseEvaluatorFactory_ExtremeValue<EvalT,LO,GO>::
buildResponseObject(const std::string & responseName) const
{ 
  Teuchos::RCP<ResponseBase> response = Teuchos::rcp(new Response_ExtremeValue<EvalT>(responseName,comm_,useMax_,linearObjFactory_)); 
  response->setRequiresDirichletAdjustment(applyDirichletToDerivative_);
 
  return response;
}

template <typename EvalT,typename LO,typename GO>
void ResponseEvaluatorFactory_ExtremeValue<EvalT,LO,GO>::
buildAndRegisterEvaluators(const std::string & responseName,
                           PHX::FieldManager<panzer::Traits> & fm,
                           const panzer::PhysicsBlock & physicsBlock,
                           const Teuchos::ParameterList & user_data) const
{
   using Teuchos::RCP;
   using Teuchos::rcp;


   // build integration evaluator (integrate over element)
   if(requiresCellExtreme_) {
     std::string field = (quadPointField_=="" ? responseName : quadPointField_);

     // build integration rule to use in cell integral
     RCP<IntegrationRule> ir = rcp(new IntegrationRule(cubatureDegree_,physicsBlock.cellData()));

     Teuchos::ParameterList pl;
     // add prefix_ to help identify
     pl.set("Extreme Name",prefix_+field);
     pl.set("Field Name",field);
     pl.set("IR",ir);
     pl.set("Use Max",useMax_);

     Teuchos::RCP<PHX::Evaluator<panzer::Traits> > eval 
         = Teuchos::rcp(new CellExtreme<EvalT,panzer::Traits>(pl));
 
     this->template registerEvaluator<EvalT>(fm, eval);
   }


   // build scatter evaluator
   {
     Teuchos::RCP<ExtremeValueScatterBase> scatterObj =
         (globalIndexer_!=Teuchos::null) ?  Teuchos::rcp(new ExtremeValueScatter<LO,GO>(globalIndexer_)) : Teuchos::null;
     std::string field = (quadPointField_=="" ? responseName : quadPointField_);
     field = prefix_+field; // add prefix to help identify

     // build useful evaluator
     Teuchos::RCP<PHX::Evaluator<panzer::Traits> > eval 
         = Teuchos::rcp(new ResponseScatterEvaluator_ExtremeValue<EvalT,panzer::Traits>(field,         
                                                                                        responseName,         
                                                                                        physicsBlock.cellData(),
                                                                                        useMax_,
                                                                                        scatterObj));

     this->template registerEvaluator<EvalT>(fm, eval);

     // require last field
     fm.template requireField<EvalT>(*eval->evaluatedFields()[0]);
   }
}

template <typename EvalT,typename LO,typename GO>
bool ResponseEvaluatorFactory_ExtremeValue<EvalT,LO,GO>::
typeSupported() const
{
  if(   PHX::typeAsString<EvalT>()==PHX::typeAsString<panzer::Traits::Residual>()  ||
        PHX::typeAsString<EvalT>()==PHX::typeAsString<panzer::Traits::Tangent>()
    )
    return true;

  if(PHX::typeAsString<EvalT>()==PHX::typeAsString<panzer::Traits::Jacobian>())
    return false;

  return false;
}

}

#endif
