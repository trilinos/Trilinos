#ifndef UNIT_VALUE_EVALUATOR
#define UNIT_VALUE_EVALUATOR

#include <string>
#include "Phalanx_Evaluator_Macros.hpp"
#include "Phalanx_Field.hpp"
#include "Panzer_IntegrationRule.hpp"

namespace panzer {

PHX_EVALUATOR_CLASS(TestEvaluator)

  PHX::MDField<ScalarT,Cell> dogValues;
  PHX::MDField<ScalarT,Cell> hrsValues;

PHX_EVALUATOR_CLASS_END

//**********************************************************************
PHX_EVALUATOR_CTOR(TestEvaluator,p)
{
  // Read from parameters
  int worksetSize = p.get<int>("Workset Size");

  Teuchos::RCP<PHX::DataLayout> dl = Teuchos::rcp(new PHX::MDALayout<Cell>(worksetSize));
  // grab information from quadrature rule
  dogValues = PHX::MDField<ScalarT,Cell>("Dog", dl);
  hrsValues = PHX::MDField<ScalarT,Cell>("Horse", dl);

  this->addEvaluatedField(dogValues);
  this->addEvaluatedField(hrsValues);
  
  std::string n = "TestEvaluator";
  this->setName(n);
}

//**********************************************************************
PHX_POST_REGISTRATION_SETUP(TestEvaluator,sd,fm)
{
  this->utils.setFieldData(dogValues,fm);
  this->utils.setFieldData(hrsValues,fm);
}

//**********************************************************************
PHX_EVALUATE_FIELDS(TestEvaluator,workset)
{ 
   double extra = 0.0;
   if(workset.block_id=="block_1")
      extra = 44.3;

   for(int i=0;i<dogValues.dimension(0);i++) {
      dogValues(i) = double(i) + 1.0 + extra;
      hrsValues(i) = -double(i) - 5.5 + extra;
   }
}

//**********************************************************************

}

#endif
