#ifndef PANZER_EVALUATOR_SCALAR_HPP
#define PANZER_EVALUATOR_SCALAR_HPP

#include <string>
#include "Panzer_Dimension.hpp"
#include "Phalanx_Evaluator_Macros.hpp"
#include "Phalanx_MDField.hpp"
#include "Intrepid_FieldContainer.hpp"

namespace panzer {
    
/** This integrates a scalar quanity over each cell.
  * It is useful for comptuing integral responses.

    <ParameterList name="Integral Name" type="string" value="<Name to give to the integral field>"/>
    <ParameterList name="Integrand Name" type="string" value="<Name of integrand>"/>
    <ParameterList name="IR" type="RCP<IntegrationRule>" value="<user specified IntegrationRule>"/>
    <ParameterList name="Multiplier" type="double" value="<Scaling factor, default=1>"/>
    <ParameterList name="Field Multipliers" type="RCP<vector<string> >" value="<Other scalar multiplier fields>"/>
  */
PHX_EVALUATOR_CLASS(Integrator_Scalar)
  
  PHX::MDField<ScalarT,Cell> integral;  // result
    
  PHX::MDField<ScalarT,Cell,IP> scalar; // function to be integrated

  std::vector<PHX::MDField<ScalarT,Cell,IP> > field_multipliers;

  std::size_t num_qp;
  std::size_t quad_index;
  int quad_order;

  double multiplier;

  Intrepid::FieldContainer<ScalarT> tmp;

public:
  // for testing purposes
  const PHX::FieldTag & getFieldTag() const 
  { return integral.fieldTag(); }

private:
  Teuchos::RCP<Teuchos::ParameterList> getValidParameters() const;

PHX_EVALUATOR_CLASS_END

}

#include "Panzer_Integrator_ScalarT.hpp"

#endif
