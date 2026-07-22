// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef PANZER_MULTIVARIATE_PARAMETER_IMPL_HPP
#define PANZER_MULTIVARIATE_PARAMETER_IMPL_HPP

#include "PanzerDiscFE_config.hpp"
#include "Panzer_ScalarParameterEntry.hpp"
#include "Panzer_ParameterLibraryUtilities.hpp"
#include <cstddef>
#include <string>
#include <vector>
#include <sstream>

namespace panzer {

//**********************************************************************
template<typename EvalT, typename TRAITS>
MultiVariateParameter<EvalT, TRAITS>::
MultiVariateParameter(const std::string parameter_name,
                      const int num_param,
                      const std::string field_name,
                      const Teuchos::RCP<PHX::DataLayout>& data_layout,
                      panzer::ParamLib& param_lib)
{
  target_field = PHX::MDField<ScalarT, Cell, Point>(field_name, data_layout);
  this->addEvaluatedField(target_field);

  param.resize(num_param);
  for (int i=0; i<num_param; ++i) {
    std::stringstream ss;
    ss << parameter_name << "_" << i;
    param[i] = panzer::createAndRegisterScalarParameter<EvalT>(
      ss.str(),param_lib);
  }

  // no initialization, this will be done by someone else (possibly the ME) later

  std::string n = "Multi-variate Parameter Evaluator";
  this->setName(n);
}

//**********************************************************************
template<typename EvalT, typename TRAITS>
void MultiVariateParameter<EvalT, TRAITS>::
evaluateFields(typename TRAITS::EvalData workset)
{
  ScalarT sum = 0;
  const int num_param = param.size();
  for (int i=0; i<num_param; ++i)
    sum += param[i]->getValue();

  for (index_t cell = 0; cell < workset.num_cells; ++cell) {
    for (typename PHX::MDField<ScalarT, Cell, Point>::size_type pt = 0;
         pt < target_field.extent(1); ++pt) {
      target_field(cell,pt) = sum;
    }
  }

}

//**********************************************************************

}

#endif
