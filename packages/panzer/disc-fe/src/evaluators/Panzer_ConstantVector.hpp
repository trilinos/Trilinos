// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef PANZER_EVALUATOR_CONSTANT_VECTOR_HPP
#define PANZER_EVALUATOR_CONSTANT_VECTOR_HPP

#include "PanzerDiscFE_config.hpp"

#include "Phalanx_Evaluator_Macros.hpp"
#include "Phalanx_MDField.hpp"

#include "Panzer_Evaluator_Macros.hpp"

namespace panzer {

/** This sets a vector field to a constant (x, y[, z]) value, the same
  * for every cell and point in the workset.
  * 
  *  NOTE: This evalautor with not work with PHX_ENABLE_SHARED=1 as
  *        the memory is overwritten by other evaluators.
  */
template<typename EvalT, typename Traits>
class ConstantVector
  : public panzer::EvaluatorWithBaseImpl<Traits>,
    public PHX::EvaluatorDerived<EvalT, Traits>
{
public:

  /// \brief Construct from a ParameterList specifying the field name, data layout, and the constant vector components.
  ConstantVector(const Teuchos::ParameterList& p);

  /// \brief Sets the vector field to the constant components given at construction. Called once before the first evaluation.
  void postRegistrationSetup(typename Traits::SetupData d,
                             PHX::FieldManager<Traits>& fm);

  /// \brief No-op; the vector field is already set by postRegistrationSetup().
  void evaluateFields(typename Traits::EvalData d);
  
private:

  using ScalarT = typename EvalT::ScalarT;
  Kokkos::View<double*> vals_;
  PHX::MDField<ScalarT> vec_;
};


}

#endif
