// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef PANZER_PARAMETER_IMPL_HPP
#define PANZER_PARAMETER_IMPL_HPP

#include "PanzerDiscFE_config.hpp"
#include "Panzer_ScalarParameterEntry.hpp"
#include "Panzer_ParameterLibraryUtilities.hpp"
#include <cstddef>
#include <string>
#include <vector>

namespace panzer {

//**********************************************************************
template<typename EvalT, typename TRAITS>
Parameter<EvalT, TRAITS>::
Parameter(const std::string parameter_name,
	  const std::string field_name,
	  const Teuchos::RCP<PHX::DataLayout>& data_layout,
	  panzer::ParamLib& param_lib)
{
  target_field = PHX::MDField<ScalarT, Cell, Point>(field_name, data_layout);

  this->addEvaluatedField(target_field);

  //param = panzer::accessScalarParameter<EvalT>(parameter_name,param_lib);
  param = panzer::createAndRegisterScalarParameter<EvalT>(parameter_name,param_lib);
    // no initialization, this will be done by someone else (possibly the ME) later

  std::string n = "Parameter Evaluator";
  this->setName(n);
}

//**********************************************************************
template<typename EvalT, typename TRAITS>
void Parameter<EvalT, TRAITS>::
evaluateFields(typename TRAITS::EvalData workset)
{
  //std::cout << "Parameter::evalauteFields() ParamValue = " << param->getValue() << std::endl;
  auto param_val = param->getValue();
  auto target_field_v = target_field.get_static_view();
  auto target_field_h = Kokkos::create_mirror_view(target_field_v);
  
  for (int cell=0; cell < workset.num_cells; ++cell) {
    for (std::size_t pt=0; pt<target_field_v.extent(1); ++pt)
      target_field_h(cell,pt) = param_val;
  }
  Kokkos::deep_copy(target_field_v, target_field_h);

}

//**********************************************************************

}

#endif
