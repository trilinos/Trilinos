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
