#include "Panzer_ScalarParameterEntry.hpp"
#include <cstddef>
#include <string>
#include <vector>

namespace panzer {

//**********************************************************************
PHX_EVALUATOR_CTOR(Parameter,p)
{
  std::string target_name = p.get<std::string>("Name");
  Teuchos::RCP<PHX::DataLayout> data_layout = 
    p.get< Teuchos::RCP<PHX::DataLayout> >("Data Layout");
  
  initial_value = ScalarT(p.get<double>("Value"));

  target_field = PHX::MDField<ScalarT, Cell, Point>(target_name, data_layout);
  
  this->addEvaluatedField(target_field);
 
  param = Teuchos::rcp(new panzer::ScalarParameterEntry<EvalT>);
  // need to register with parmeter library

  std::string n = "Parameter Evaluator";
  this->setName(n);
}

//**********************************************************************
PHX_POST_REGISTRATION_SETUP(Parameter,worksets,fm)
{
  this->utils.setFieldData(target_field,fm);
}

//**********************************************************************
PHX_EVALUATE_FIELDS(Parameter,workset)
{ 
  for (std::size_t cell = 0; cell < workset.num_cells; ++cell) {
    for (typename PHX::MDField<ScalarT, Cell, Point>::size_type pt = 0;
	 pt < target_field.dimension(1); ++pt) {
      target_field(cell,pt) = initial_value;
    }
  }

}

//**********************************************************************

}
