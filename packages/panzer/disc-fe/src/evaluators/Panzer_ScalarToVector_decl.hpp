// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef PANZER_EVALUATOR_SCALAR_TO_VECTOR_DECL_HPP
#define PANZER_EVALUATOR_SCALAR_TO_VECTOR_DECL_HPP

#include <vector>
#include "Phalanx_Evaluator_WithBaseImpl.hpp"
#include "Panzer_Evaluator_WithBaseImpl.hpp"
#include "Phalanx_MDField.hpp"
#include "Phalanx_KokkosViewOfViews.hpp"

namespace panzer {

//! Copies a set of scalar fields into a single vector field
template<typename EvalT, typename Traits>
class ScalarToVector : public panzer::EvaluatorWithBaseImpl<Traits>
{
public:
  using ScalarT = typename EvalT::ScalarT;

  ScalarToVector(const Teuchos::ParameterList& p);

  ScalarToVector(const std::vector<PHX::Tag<ScalarT>> & input,
                 const PHX::FieldTag & output);

  void postRegistrationSetup(typename Traits::SetupData d,
                             PHX::FieldManager<Traits>& fm);

  void evaluateFields(typename Traits::EvalData d);

private:

  PHX::MDField<ScalarT,Cell,Point,Dim> vector_field_;
  std::vector<PHX::MDField<const ScalarT,Cell,Point>> scalar_fields_;
  using InnerView = typename PHX::MDField<const ScalarT,Cell,Point>::array_type;
  PHX::ViewOfViews<1,InnerView> scalar_fields_vov_;
};

}

#endif
