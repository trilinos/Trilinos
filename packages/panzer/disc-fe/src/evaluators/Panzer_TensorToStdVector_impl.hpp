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

#ifndef PANZER_TENSOR_TO_STD_VECTOR_IMPL_HPP
#define PANZER_TENSOR_TO_STD_VECTOR_IMPL_HPP

#include <string>

namespace panzer {

//**********************************************************************
template<typename EvalT, typename Traits>
TensorToStdVector<EvalT, Traits>::
TensorToStdVector(
  const Teuchos::ParameterList&  p)
{
  Teuchos::RCP<PHX::DataLayout> vector_dl = 
    p.get< Teuchos::RCP<PHX::DataLayout> >("Data Layout Vector");

  Teuchos::RCP<PHX::DataLayout> tensor_dl =
    p.get< Teuchos::RCP<PHX::DataLayout> >("Data Layout Tensor");

  const std::vector<std::string>& vector_names =
    *(p.get< Teuchos::RCP<const std::vector<std::string> > >("Vector Names"));

  vector_fields.resize(vector_names.size());
  for (std::size_t i=0; i < vector_names.size(); ++i)
    vector_fields[i] =
      PHX::MDField<ScalarT,Cell,Point,Dim>(vector_names[i], vector_dl);

  tensor_field =
    PHX::MDField<const ScalarT,Cell,Point,Dim,Dim>(p.get<std::string>
					 ("Tensor Name"), tensor_dl);

  this->addDependentField(tensor_field);
  
  for (std::size_t i=0; i < vector_fields.size(); ++i)
    this->addEvaluatedField(vector_fields[i]);
  
  std::string n = "TensorToStdVector: " + tensor_field.fieldTag().name();
  this->setName(n);
}

//**********************************************************************
template<typename EvalT, typename Traits>
void
TensorToStdVector<EvalT, Traits>::
postRegistrationSetup(
  typename Traits::SetupData  /* worksets */,
  PHX::FieldManager<Traits>&  fm)
{
  for (std::size_t i=0; i < vector_fields.size(); ++i)
    this->utils.setFieldData(vector_fields[i], fm);

  this->utils.setFieldData(tensor_field, fm);
}

//**********************************************************************
template<typename EvalT, typename Traits>
void
TensorToStdVector<EvalT, Traits>::
evaluateFields(
  typename Traits::EvalData  workset)
{ 

  typedef typename PHX::MDField<ScalarT,Cell,Point,Dim>::size_type size_type;

  // Loop over cells
  for (index_t cell = 0; cell < workset.num_cells; ++cell) {

    // Loop over points
    for (size_type pt = 0; pt < tensor_field.dimension(1); ++pt) {
      
      // Loop over vectors
      for (std::size_t vec = 0; vec < vector_fields.size(); ++vec) {

        // Loop over spatial dimensions
        for (std::size_t dim = 0; dim < tensor_field.dimension(2); ++dim) {
      
          vector_fields[vec](cell,pt,dim) = tensor_field(cell,pt,vec,dim);

        }
      }
    }
  }
}

//**********************************************************************

}

#endif
