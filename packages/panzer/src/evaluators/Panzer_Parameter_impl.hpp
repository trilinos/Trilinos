#ifndef PANZER_PARAMETER_IMPL_HPP
#define PANZER_PARAMETER_IMPL_HPP

#include "Panzer_config.hpp"
#include "Panzer_ScalarParameterEntry.hpp"
#include "Panzer_ParameterLibraryUtilities.hpp"
#include <cstddef>
#include <string>
#include <vector>

#ifdef HAVE_STOKHOS
#include "Panzer_SGUtilities.hpp"
#endif

namespace panzer {

//**********************************************************************
template<typename EvalT, typename Traits>
Parameter<EvalT, Traits>::
Parameter(const std::string name,
	  const Teuchos::RCP<PHX::DataLayout>& data_layout,
	  const double in_initial_value,
	  panzer::ParamLib& param_lib)
{ 
  initial_value = ScalarT(in_initial_value);

  target_field = PHX::MDField<ScalarT, Cell, Point>(name, data_layout);
  
  this->addEvaluatedField(target_field);
 
  param = panzer::createAndRegisterScalarParameter<EvalT>(name,param_lib);

  std::string n = "Parameter Evaluator";
  this->setName(n);
}

//**********************************************************************
#ifdef HAVE_STOKHOS

template<typename EvalT, typename Traits>
Parameter<EvalT, Traits>::
Parameter(const std::string name,
	  const Teuchos::RCP<PHX::DataLayout>& data_layout,
	  const std::vector<double> & in_initial_value,
          const Teuchos::RCP<Stokhos::OrthogPolyExpansion<int,double> > & expansion,
	  panzer::ParamLib& param_lib)
{ 
  // using expansion convert vector to a scalar value
  sg_utils::vectorToValue(in_initial_value,expansion,initial_value); 

  target_field = PHX::MDField<ScalarT, Cell, Point>(name, data_layout);
  
  this->addEvaluatedField(target_field);
 
  param = panzer::createAndRegisterScalarParameter<EvalT>(name,param_lib);

  std::string n = "Parameter Evaluator";
  this->setName(n);
}

#endif

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
      target_field(cell,pt) = param->getValue();
    }
  }

}

//**********************************************************************

}

#endif
