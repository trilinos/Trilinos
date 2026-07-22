// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
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
evaluateFields(
  typename Traits::EvalData workset)
{ 

  // Iteration bounds
  const int num_points = vector_field.extent_int(1);
  const int num_scalars = std::min(static_cast<int>(scalar_fields.size()),
                                   vector_field.extent_int(2));

  // Need local copies for cuda, *this is not usable
  auto local_vector_field = vector_field;

  // We parallelize over each scalar field
  for (int sc = 0; sc < num_scalars; ++sc) {
    auto local_scalar_field = scalar_fields[sc];
    Kokkos::parallel_for(workset.num_cells, KOKKOS_LAMBDA (const int cell) {
      for (int pt = 0; pt < num_points; ++pt)
        local_scalar_field(cell,pt) = local_vector_field(cell,pt,sc);
    });
  }

  // If there are remaining fields, set them to zero
  for(unsigned int sc = num_scalars; sc < scalar_fields.size(); ++sc)
    Kokkos::deep_copy(scalar_fields[sc].get_static_view(), 0.);
}

//**********************************************************************

}

#endif
