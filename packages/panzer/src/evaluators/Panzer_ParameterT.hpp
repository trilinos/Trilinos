#include "Panzer_ScalarParameterEntry.hpp"
#include <cstddef>
#include <string>
#include <vector>

namespace panzer {

//**********************************************************************
template<typename EvalT, typename Traits>
Parameter<EvalT, Traits>::
Parameter(const std::string name,
	  const Teuchos::RCP<PHX::DataLayout>& data_layout,
	  const double in_initial_value)
{ 
  initial_value = ScalarT(in_initial_value);

  target_field = PHX::MDField<ScalarT, Cell, Point>(name, data_layout);
  
  this->addEvaluatedField(target_field);
 
  param = Teuchos::rcp(new panzer::ScalarParameterEntry<EvalT>);
  // need to register with parmeter library

  std::string n = "Parameter Evaluator";
  this->setName(n);
}

//**********************************************************************
template<typename EvalT, typename Traits>
void Parameter<EvalT, Traits>::
postRegistrationSetup(typename Traits::SetupData worksets,
		      PHX::FieldManager<Traits>& fm)
{
  this->utils.setFieldData(target_field,fm);
}

//**********************************************************************
template<typename EvalT, typename Traits>
void Parameter<EvalT, Traits>::
evaluateFields(typename Traits::EvalData workset)
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
