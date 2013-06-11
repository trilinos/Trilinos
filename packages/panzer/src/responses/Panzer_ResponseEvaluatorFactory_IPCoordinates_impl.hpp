#ifndef __Panzer_ResponseEvaluatorFactory_IPCoordinates_impl_hpp__
#define __Panzer_ResponseEvaluatorFactory_IPCoordinates_impl_hpp__

#include <string>
#include <sstream>

#include "Panzer_config.hpp"

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
                           const panzer::PhysicsBlock & physicsBlock,
                           const Teuchos::ParameterList & user_data) const
{
   using Teuchos::RCP;
   using Teuchos::rcp;

   // build scatter evaluator
   {
     // build useful evaluator
     Teuchos::RCP<PHX::Evaluator<panzer::Traits> > eval 
         = Teuchos::rcp(new ResponseScatterEvaluator_IPCoordinates<EvalT,panzer::Traits>(responseName,cubatureDegree_));

     fm.template registerEvaluator<EvalT>(eval);

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
