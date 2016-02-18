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
     pl.set("Extreme Name",field);
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
