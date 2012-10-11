#ifndef __Panzer_ResponseEvaluatorFactory_Functional_impl_hpp__
#define __Panzer_ResponseEvaluatorFactory_Functional_impl_hpp__

#include <string>

#include "Panzer_config.hpp"
#include "Panzer_ResponseScatterEvaluator_Functional.hpp"
#include "Panzer_Response_Functional.hpp"

namespace panzer {

template <typename EvalT>
Teuchos::RCP<ResponseBase> ResponseEvaluatorFactory_Functional<EvalT>::
buildResponseObject(const std::string & responseName) const
{ return Teuchos::rcp(new Response_Functional<typename EvalT::ScalarT>(responseName)); }

template <typename EvalT>
Teuchos::RCP<const PHX::FieldTag> ResponseEvaluatorFactory_Functional<EvalT>::
buildAndRegisterEvaluators(const std::string & responseName,
                           PHX::FieldManager<panzer::Traits> & fm,
                           const panzer::PhysicsBlock & physicsBlock,
                           const Teuchos::ParameterList & user_data) const
{
   // build useful evaluator
   Teuchos::RCP<PHX::Evaluator<panzer::Traits> > eval 
       = Teuchos::rcp(new ResponseScatterEvaluator_Functional<EvalT,panzer::Traits>(responseName));

   fm.template registerEvaluator<EvalT>(eval);

   // require last field
   return eval->evaluatedFields()[0];
}

}

#endif
