// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef PANZER_EVALUATOR_DIRICHLET_RESIDUAL_DECL_HPP
#define PANZER_EVALUATOR_DIRICHLET_RESIDUAL_DECL_HPP

#include "Phalanx_Evaluator_Macros.hpp"
#include "Phalanx_MDField.hpp"

#include "Panzer_Evaluator_Macros.hpp"

namespace panzer {
    
//! Evaluates a Dirichlet BC residual corresponding to a field value
template<typename EvalT, typename Traits>
class DirichletResidual
  :
  public panzer::EvaluatorWithBaseImpl<Traits>,
  public PHX::EvaluatorDerived<EvalT, Traits>
{
  public:

    DirichletResidual(
      const Teuchos::ParameterList& p);

    void
    postRegistrationSetup(
      typename Traits::SetupData d,
      PHX::FieldManager<Traits>& fm);

    void
    evaluateFields(
      typename Traits::EvalData d);

  private:

    using ScalarT = typename EvalT::ScalarT;
  
  PHX::MDField<ScalarT> residual;
  PHX::MDField<const ScalarT> dof;
  PHX::MDField<const ScalarT> value;

  std::size_t cell_data_size;

}; // end of class DirichletResidual


}

#endif
