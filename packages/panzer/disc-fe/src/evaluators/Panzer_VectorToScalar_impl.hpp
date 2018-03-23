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

#ifndef PANZER_VECTOR_TO_SCALAR_IMPL_HPP
#define PANZER_VECTOR_TO_SCALAR_IMPL_HPP

#include <string>

namespace panzer {

//**********************************************************************
template<typename EvalT, typename Traits>
VectorToScalar<EvalT, Traits>::
VectorToScalar(
  const Teuchos::ParameterList& p)
{
  Teuchos::RCP<PHX::DataLayout> scalar_dl = 
    p.get< Teuchos::RCP<PHX::DataLayout> >("Data Layout Scalar");

  Teuchos::RCP<PHX::DataLayout> vector_dl = 
    p.get< Teuchos::RCP<PHX::DataLayout> >("Data Layout Vector");

  const std::vector<std::string>& scalar_names = 
    *(p.get< Teuchos::RCP<const std::vector<std::string> > >("Scalar Names"));

  scalar_fields.resize(scalar_names.size());
  for (std::size_t i=0; i < scalar_names.size(); ++i)
    scalar_fields[i] = 
      PHX::MDField<ScalarT,Cell,Point>(scalar_names[i], scalar_dl);

  vector_field = 
    PHX::MDField<const ScalarT,Cell,Point,Dim>(p.get<std::string>
					 ("Vector Name"), vector_dl);

  this->addDependentField(vector_field);
  
  for (std::size_t i=0; i < scalar_fields.size(); ++i)
    this->addEvaluatedField(scalar_fields[i]);
  
  std::string n = "VectorToScalar: " + vector_field.fieldTag().name();
  this->setName(n);
}

//**********************************************************************

template<typename EvalT, typename Traits>				\
VectorToScalar<EvalT,Traits>::
VectorToScalar(const PHX::FieldTag & input,
               const std::vector<PHX::Tag<ScalarT>> & output)
{
  // setup the fields
  vector_field = input;

  scalar_fields.resize(output.size());
  for(std::size_t i=0;i<output.size();i++) 
    scalar_fields[i] = output[i];

  // add dependent/evaluate fields
  this->addDependentField(vector_field);
  
  for (std::size_t i=0; i < scalar_fields.size(); ++i)
    this->addEvaluatedField(scalar_fields[i]);
  
  // name array
  std::string n = "VectorToScalar: " + vector_field.fieldTag().name();
  this->setName(n);
}

//**********************************************************************
template<typename EvalT, typename Traits>
void
VectorToScalar<EvalT, Traits>::
postRegistrationSetup(
  typename Traits::SetupData  /* worksets */,
  PHX::FieldManager<Traits>&  fm)
{
  for (std::size_t i=0; i < scalar_fields.size(); ++i)
    this->utils.setFieldData(scalar_fields[i],fm);

  this->utils.setFieldData(vector_field,fm);
}

//**********************************************************************
template<typename EvalT, typename Traits>
void
VectorToScalar<EvalT, Traits>::
evaluateFields(
  typename Traits::EvalData workset)
{ 

  typedef typename PHX::MDField<ScalarT,Cell,Point>::size_type size_type;

  // Need local copies for cuda, *this is not usable
  auto local_vector_field = vector_field;
  // Loop over scalars
  for (std::size_t sc = 0; sc < scalar_fields.size(); ++sc) {
    auto local_scalar_field = scalar_fields[sc];
    // Loop over cells
    Kokkos::parallel_for(workset.num_cells, KOKKOS_LAMBDA (const index_t cell) {
      // Loop over points
      for (size_type pt = 0; pt < local_vector_field.dimension(1); ++pt) {
        local_scalar_field(cell,pt) = local_vector_field(cell,pt,sc);
      }
    });
  }
}

//**********************************************************************

}

#endif
