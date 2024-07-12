// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
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
evaluateFields(
  typename Traits::EvalData  workset)
{ 

  typedef typename PHX::MDField<ScalarT,Cell,Point,Dim>::size_type size_type;

  // Loop over cells
  for (index_t cell = 0; cell < workset.num_cells; ++cell) {

    // Loop over points
    for (size_type pt = 0; pt < tensor_field.extent(1); ++pt) {
      
      // Loop over vectors
      for (std::size_t vec = 0; vec < vector_fields.size(); ++vec) {

        // Loop over spatial dimensions
        for (std::size_t dim = 0; dim < tensor_field.extent(2); ++dim) {
      
          vector_fields[vec](cell,pt,dim) = tensor_field(cell,pt,vec,dim);

        }
      }
    }
  }
}

//**********************************************************************

}

#endif
