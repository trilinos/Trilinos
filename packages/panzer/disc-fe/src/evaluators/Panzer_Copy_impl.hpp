// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef PANZER_COPY_IMPL_HPP
#define PANZER_COPY_IMPL_HPP

#include <cstddef>
#include <string>
#include <vector>

namespace panzer {

//**********************************************************************
template<typename EvalT, typename Traits>
Copy<EvalT, Traits>::
Copy(
  const Teuchos::ParameterList& p)
{
  std::string input_name = p.get<std::string>("Source Name");
  std::string output_name = p.get<std::string>("Destination Name");
  Teuchos::RCP<PHX::DataLayout> data_layout = p.get< Teuchos::RCP<PHX::DataLayout> >("Data Layout");
  
  input  = PHX::MDField<const ScalarT>(input_name, data_layout);
  output = PHX::MDField<ScalarT>(output_name, data_layout);
  
  this->addDependentField(input);
  this->addEvaluatedField(output);
 
  std::string n = "Copy Evaluator: " + input_name + " => " + output_name;
  this->setName(n);
}

//**********************************************************************
template<typename EvalT, typename Traits>
void
Copy<EvalT, Traits>::
postRegistrationSetup(
  typename Traits::SetupData /* worksets */,
  PHX::FieldManager<Traits>& /* fm */)
{
  TEUCHOS_ASSERT(input.size()==output.size());
}

//**********************************************************************
template<typename EvalT, typename Traits>
void
Copy<EvalT, Traits>::
evaluateFields(
  typename Traits::EvalData  /* workset */)
{ 
  output.deep_copy(input);
}

//**********************************************************************

}

#endif
