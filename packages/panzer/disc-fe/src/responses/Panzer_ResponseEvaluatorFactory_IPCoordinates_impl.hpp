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

#ifndef __Panzer_ResponseEvaluatorFactory_IPCoordinates_impl_hpp__
#define __Panzer_ResponseEvaluatorFactory_IPCoordinates_impl_hpp__

#include <string>
#include <sstream>

#include "PanzerDiscFE_config.hpp"

#include "Panzer_IntegrationRule.hpp"
#include "Panzer_PhysicsBlock.hpp"
#include "Panzer_Integrator_Scalar.hpp"
#include "Panzer_ResponseScatterEvaluator_IPCoordinates.hpp"
#include "Panzer_Response_IPCoordinates.hpp"

namespace panzer {

template <typename EvalT>
Teuchos::RCP<ResponseBase> ResponseEvaluatorFactory_IPCoordinates<EvalT>::
buildResponseObject(const std::string & responseName,const std::vector<WorksetDescriptor> & wkstDesc) const
{ 
  // check that the input worksets constains only element blocks 
  bool failure = false;
  std::stringstream failureStrm;
  for(std::size_t i=0;i<wkstDesc.size();i++) {
    failure |= wkstDesc[i].useSideset();
    failureStrm << wkstDesc[i] << std::endl;
  }
  TEUCHOS_TEST_FOR_EXCEPTION(failure,std::runtime_error,
                             "REF_IPCoordinates::buildResponseObject: could not build using side set descriptors:\n"
                             << failureStrm.str());

  return Teuchos::rcp(new Response_IPCoordinates<EvalT>(responseName)); 
}

template <typename EvalT>
void ResponseEvaluatorFactory_IPCoordinates<EvalT>::
buildAndRegisterEvaluators(const std::string & responseName,
                           PHX::FieldManager<panzer::Traits> & fm,
                           const panzer::PhysicsBlock& /* physicsBlock */,
                           const Teuchos::ParameterList& /* user_data */) const
{
   using Teuchos::RCP;
   using Teuchos::rcp;

   // build scatter evaluator
   {
     // build useful evaluator
     Teuchos::RCP<PHX::Evaluator<panzer::Traits> > eval 
         = Teuchos::rcp(new ResponseScatterEvaluator_IPCoordinates<EvalT,panzer::Traits>(responseName,cubatureDegree_));

     this->template registerEvaluator<EvalT>(fm, eval);

     // require last field
     fm.template requireField<EvalT>(*eval->evaluatedFields()[0]);
   }
}

template <typename EvalT>
bool ResponseEvaluatorFactory_IPCoordinates<EvalT>::
typeSupported() const
{
   return false;
}

template < >
bool ResponseEvaluatorFactory_IPCoordinates<panzer::Traits::Residual>::
typeSupported() const
{
  return true;
}

}

#endif
