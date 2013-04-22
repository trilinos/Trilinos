#ifndef __Panzer_ResponseEvaluatorFactory_Functional_impl_hpp__
#define __Panzer_ResponseEvaluatorFactory_Functional_impl_hpp__

#include <string>

#include "Panzer_config.hpp"

#include "Panzer_IntegrationRule.hpp"
#include "Panzer_PhysicsBlock.hpp"
#include "Panzer_Integrator_Scalar.hpp"
#include "Panzer_ResponseScatterEvaluator_Functional.hpp"
#include "Panzer_Response_Functional.hpp"

namespace panzer {

template <typename EvalT,typename LO,typename GO>
Teuchos::RCP<ResponseBase> ResponseEvaluatorFactory_Functional<EvalT,LO,GO>::
buildResponseObject(const std::string & responseName) const
{ return Teuchos::rcp(new Response_Functional<EvalT>(responseName,comm_,linearObjFactory_)); }

template <typename EvalT,typename LO,typename GO>
void ResponseEvaluatorFactory_Functional<EvalT,LO,GO>::
buildAndRegisterEvaluators(const std::string & responseName,
                           PHX::FieldManager<panzer::Traits> & fm,
                           const panzer::PhysicsBlock & physicsBlock,
                           const Teuchos::ParameterList & user_data) const
{
   using Teuchos::RCP;
   using Teuchos::rcp;


   // build integration evaluator (integrate over element)
   if(requiresCellIntegral_) {
     std::string field = (quadPointField_=="" ? responseName : quadPointField_);

     // build integration rule to use in cell integral
     RCP<IntegrationRule> ir = rcp(new IntegrationRule(cubatureDegree_,physicsBlock.cellData()));

     Teuchos::ParameterList pl;
     pl.set("Integral Name",field);
     pl.set("Integrand Name",field);
     pl.set("IR",ir);

     Teuchos::RCP<PHX::Evaluator<panzer::Traits> > eval 
         = Teuchos::rcp(new Integrator_Scalar<EvalT,panzer::Traits>(pl));
 
     fm.template registerEvaluator<EvalT>(eval);
   }


   // build scatter evaluator
   {
     Teuchos::RCP<FunctionalScatterBase> scatterObj =
         (globalIndexer_!=Teuchos::null) ?  Teuchos::rcp(new FunctionalScatter<LO,GO>(globalIndexer_)) : Teuchos::null;
     std::string field = (quadPointField_=="" ? responseName : quadPointField_);

     // build useful evaluator
     Teuchos::RCP<PHX::Evaluator<panzer::Traits> > eval 
         = Teuchos::rcp(new ResponseScatterEvaluator_Functional<EvalT,panzer::Traits>(field,responseName,physicsBlock.cellData(),scatterObj));

     fm.template registerEvaluator<EvalT>(eval);

     // require last field
     fm.template requireField<EvalT>(*eval->evaluatedFields()[0]);
   }
}

template <typename EvalT,typename LO,typename GO>
bool ResponseEvaluatorFactory_Functional<EvalT,LO,GO>::
typeSupported() const
{
  if(   PHX::TypeString<EvalT>::value==PHX::TypeString<panzer::Traits::Residual>::value 
#ifdef HAVE_STOKHOS
     || PHX::TypeString<EvalT>::value==PHX::TypeString<panzer::Traits::SGResidual>::value
#endif
    )
    return true;

  if(PHX::TypeString<EvalT>::value==PHX::TypeString<panzer::Traits::Jacobian>::value)
    return linearObjFactory_!=Teuchos::null;

  return false;
}

}

#endif
