#ifndef UNIT_VALUE_EVALUATOR
#define UNIT_VALUE_EVALUATOR

#include <string>
#include "Phalanx_Evaluator_Macros.hpp"
#include "Phalanx_Field.hpp"
#include "Panzer_IntegrationRule.hpp"

namespace panzer {

PHX_EVALUATOR_CLASS(UnitValueEvaluator)

  PHX::MDField<ScalarT,Cell,panzer::IP> unitValue;

PHX_EVALUATOR_CLASS_END

//**********************************************************************
PHX_EVALUATOR_CTOR(UnitValueEvaluator,p)
{
  // Read from parameters
  const std::string name = p.get<std::string>("Name");
  Teuchos::RCP<panzer::IntegrationRule> ir
     = p.get< Teuchos::RCP<panzer::IntegrationRule> >("IR");

  // grab information from quadrature rule
  unitValue = PHX::MDField<ScalarT,Cell,IP>(name, ir->dl_scalar);

  this->addEvaluatedField(unitValue);
  
  std::string n = "UnitValueEvaluator: " + name;
  this->setName(n);
}

//**********************************************************************
PHX_POST_REGISTRATION_SETUP(UnitValueEvaluator,sd,fm)
{
  this->utils.setFieldData(unitValue,fm);
}

//**********************************************************************
PHX_EVALUATE_FIELDS(UnitValueEvaluator,workset)
{ 
   for(int i=0;i<unitValue.size();i++)
      unitValue[i] = 1.0;
}

//**********************************************************************

}

#endif
