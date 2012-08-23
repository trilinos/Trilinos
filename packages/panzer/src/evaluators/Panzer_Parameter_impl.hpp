// @HEADER
// ***********************************************************************
//
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//                 Copyright (2011) Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Roger P. Pawlowski (rppawlo@sandia.gov) and
// Eric C. Cyr (eccyr@sandia.gov)
// ***********************************************************************
// @HEADER

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
