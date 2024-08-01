// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef PANZER_SCALAR_TO_VECTOR_IMPL_HPP
#define PANZER_SCALAR_TO_VECTOR_IMPL_HPP

#include <string>

namespace panzer {

//**********************************************************************
template<typename EvalT, typename Traits>
ScalarToVector<EvalT, Traits>::
ScalarToVector(
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
      PHX::MDField<const ScalarT,Cell,Point>(scalar_names[i], scalar_dl);

  vector_field = 
    PHX::MDField<ScalarT,Cell,Point,Dim>(p.get<std::string>
					 ("Vector Name"), vector_dl);

  for (std::size_t i=0; i < scalar_fields.size(); ++i)
    this->addDependentField(scalar_fields[i]);
  
  this->addEvaluatedField(vector_field);
  
  std::string n = "ScalarToVector: " + vector_field.fieldTag().name();
  this->setName(n);
}

//**********************************************************************

template<typename EvalT, typename Traits>				\
ScalarToVector<EvalT,Traits>::
ScalarToVector(const std::vector<PHX::Tag<ScalarT>> & input,
               const PHX::FieldTag & output)
{
  // setup the fields
  vector_field = output;

  scalar_fields.resize(input.size());
  for(std::size_t i=0;i<input.size();i++) 
    scalar_fields[i] = input[i];

  // add dependent/evaluate fields
  this->addEvaluatedField(vector_field);
  
  for (std::size_t i=0; i < scalar_fields.size(); ++i)
    this->addDependentField(scalar_fields[i]);
  
  // name array
  std::string n = "ScalarToVector: " + vector_field.fieldTag().name();
  this->setName(n);
}

//**********************************************************************
template<typename EvalT, typename Traits>
void
ScalarToVector<EvalT, Traits>::
postRegistrationSetup(
  typename Traits::SetupData  /* worksets */,
  PHX::FieldManager<Traits>&  fm)
{
  // Convert std::vector to PHX::View for use on device
  internal_scalar_fields = PHX::View<KokkosScalarFields_t*>("ScalarToVector::internal_scalar_fields", scalar_fields.size());
  for (std::size_t i=0; i < scalar_fields.size(); ++i)
    internal_scalar_fields(i) = scalar_fields[i].get_static_view();
}

//**********************************************************************
template<typename EvalT, typename Traits>
void
ScalarToVector<EvalT, Traits>::
evaluateFields(
  typename Traits::EvalData workset)
{ 

  using Scalar = typename EvalT::ScalarT;

  // Iteration bounds
  const int num_points = vector_field.extent_int(1);
  const int num_vector_scalars = vector_field.extent_int(2);
  const int num_scalars = std::min(internal_scalar_fields.extent_int(0),
                                   num_vector_scalars);

  // Local copies to prevent passing (*this) to device code
  auto vector = vector_field;
  auto scalars = internal_scalar_fields;

  // Loop over cells, points, scalars
  Kokkos::parallel_for (workset.num_cells,KOKKOS_LAMBDA (const int cell){
    for (int pt = 0; pt < num_points; ++pt){

      // Copy over scalars
      for (int sc = 0; sc < num_scalars; ++sc)
        vector(cell,pt,sc) = scalars(sc)(cell,pt);

      // Missing scalars are filled with zero
      for(int sc = num_scalars; sc < num_vector_scalars; ++sc)
        vector(cell,pt,sc) = Scalar(0);
    }
  });

}

//**********************************************************************

}

#endif
