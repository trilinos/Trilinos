// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef PANZER_EVALUATOR_CONSTANT_FLUX_DECL_HPP
#define PANZER_EVALUATOR_CONSTANT_FLUX_DECL_HPP

#include "PanzerDiscFE_config.hpp"

#include "Phalanx_Evaluator_Macros.hpp"
#include "Phalanx_MDField.hpp"

#include "Panzer_Evaluator_Macros.hpp"

namespace panzer {

/** This sets a flux field to a constant, per-dimension vector of values,
  * the same for every cell and integration point in the workset.
  * 
  * NOTE: This implementation will not work with PHX_ENABLE_SHARED=1 as
  *       the memory is potentially reused and overwritten by other 
  *       evaluators.
  */
template<typename EvalT, typename Traits>
class ConstantFlux
  :
  public panzer::EvaluatorWithBaseImpl<Traits>,
  public PHX::EvaluatorDerived<EvalT, Traits>
{
  public:

    /// \brief Construct from a ParameterList specifying the flux field name, data layout, and constant per-dimension flux values.
    ConstantFlux(
      const Teuchos::ParameterList& p);

    /// \brief Sets the flux field to the constant values given at construction. Called once before the first evaluation.
    void
    postRegistrationSetup(
      typename Traits::SetupData d,
      PHX::FieldManager<Traits>& fm);

    /// \brief No-op; the flux field is already set by postRegistrationSetup().
    void
    evaluateFields(
      typename Traits::EvalData d);

  private:

    using ScalarT = typename EvalT::ScalarT;
  
  std::vector<ScalarT> values;
  
  PHX::MDField<ScalarT> flux;

}; // end of class ConstantFlux


}

#endif
