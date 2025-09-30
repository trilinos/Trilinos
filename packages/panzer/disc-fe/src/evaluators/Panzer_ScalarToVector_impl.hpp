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
ScalarToVector(const Teuchos::ParameterList& p)
{
  auto scalar_dl = p.get<Teuchos::RCP<PHX::DataLayout>>("Data Layout Scalar");
  auto vector_dl = p.get< Teuchos::RCP<PHX::DataLayout> >("Data Layout Vector");
  auto scalar_names = p.get<Teuchos::RCP<const std::vector<std::string>>>("Scalar Names");
  auto vector_name = p.get<std::string>("Vector Name");

  scalar_fields_.resize(scalar_names->size());
  for (size_t i=0; i < scalar_names->size(); ++i)
    scalar_fields_[i] = PHX::MDField<const ScalarT,Cell,Point>((*scalar_names)[i], scalar_dl);

  vector_field_ = PHX::MDField<ScalarT,Cell,Point,Dim>(vector_name, vector_dl);

  for (std::size_t i=0; i < scalar_fields_.size(); ++i)
    this->addDependentField(scalar_fields_[i]);

  this->addEvaluatedField(vector_field_);

  std::string n = "ScalarToVector: " + vector_field_.fieldTag().name();
  this->setName(n);
}

//**********************************************************************

template<typename EvalT, typename Traits>				\
ScalarToVector<EvalT,Traits>::
ScalarToVector(const std::vector<PHX::Tag<ScalarT>> & input,
               const PHX::FieldTag & output)
{
  // setup the fields
  vector_field_ = output;

  scalar_fields_.resize(input.size());
  for(std::size_t i=0;i<input.size();i++)
    scalar_fields_[i] = input[i];

  // add dependent/evaluate fields
  this->addEvaluatedField(vector_field_);

  for (std::size_t i=0; i < scalar_fields_.size(); ++i)
    this->addDependentField(scalar_fields_[i]);

  // name array
  std::string n = "ScalarToVector: " + vector_field_.fieldTag().name();
  this->setName(n);
}

//**********************************************************************
template<typename EvalT, typename Traits>
void
ScalarToVector<EvalT, Traits>::
postRegistrationSetup(
  typename Traits::SetupData  /* worksets */,
  PHX::FieldManager<Traits>&  /* fm */ )
{
  scalar_fields_vov_.initialize("Scalar Fields VoV",scalar_fields_.size());
  for (size_t i=0; i < scalar_fields_.size(); ++i)
    scalar_fields_vov_.addView(scalar_fields_[i].get_static_view(),i);

  scalar_fields_vov_.syncHostToDevice();
}

//**********************************************************************
template<typename EvalT, typename Traits>
void
ScalarToVector<EvalT, Traits>::
evaluateFields(
  typename Traits::EvalData workset)
{
  const int num_points = vector_field_.extent_int(1);
  const int num_vector_scalars = vector_field_.extent_int(2);
  auto vec = vector_field_;
  auto scalars = scalar_fields_vov_.getViewDevice();

  // Corner case: user can supply fewer scalar fields than the
  // dimension on vector. If so, fill missing fields with zeroes.
  const int num_scalars = std::min(static_cast<int>(scalars.extent(0)),num_vector_scalars);

  Kokkos::parallel_for (workset.num_cells,KOKKOS_LAMBDA (const int cell){
    for (int pt = 0; pt < num_points; ++pt) {
      for (int sc = 0; sc < num_scalars; ++sc) {
        vec(cell,pt,sc) = scalars(sc)(cell,pt);
      }

      // Missing scalars are filled with zero
      for(int sc = num_scalars; sc < num_vector_scalars; ++sc) {
        vec(cell,pt,sc) = ScalarT(0.0);
      }
    }
  });
}

//**********************************************************************

}

#endif
