// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef PANZER_COORDINATESEVALUTOR_HPP
#define PANZER_COORDINATESEVALUTOR_HPP

#include "PanzerDiscFE_config.hpp"

#include "Phalanx_Evaluator_Macros.hpp"
#include "Phalanx_MDField.hpp"

#include "Panzer_Evaluator_Macros.hpp"

namespace panzer {

/** This copies one spatial dimension of each cell's node coordinates
  * into a basis-indexed field. It is useful for exposing mesh
  * coordinates as an evaluated field, e.g. for coordinate-dependent
  * closure models.
  */
template<typename EvalT, typename Traits>
class CoordinatesEvaluator
  :
  public panzer::EvaluatorWithBaseImpl<Traits>,
  public PHX::EvaluatorDerived<EvalT, Traits>
{
  public:

    /// \brief Construct from a ParameterList specifying the field name, data layout, and which spatial dimension to extract.
    CoordinatesEvaluator(
      const Teuchos::ParameterList& p);

    /// \brief Sets up the field data handle used by evaluateFields(). Called once before the first evaluation.
    void
    postRegistrationSetup(
      typename Traits::SetupData d,
      PHX::FieldManager<Traits>& fm);

    /// \brief Copies the selected spatial dimension of this workset's cell node coordinates into the field.
    void
    evaluateFields(
      typename Traits::EvalData d);

  private:

    using ScalarT = typename EvalT::ScalarT;
  
  int dimension;
  
  PHX::MDField<ScalarT,Cell,BASIS> coordinate;
  
}; // end of class CoordinatesEvaluator


}

#endif
