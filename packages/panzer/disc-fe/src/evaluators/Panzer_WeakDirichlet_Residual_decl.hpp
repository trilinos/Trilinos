// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef PANZER_EVALUATOR_WEAKDIRICHLET_RESIDUAL_DECL_HPP
#define PANZER_EVALUATOR_WEAKDIRICHLET_RESIDUAL_DECL_HPP

#include "Phalanx_Evaluator_Macros.hpp"
#include "Phalanx_MDField.hpp"

#include "Panzer_Evaluator_Macros.hpp"

namespace panzer {
    
  /** \brief Evaluates a Weak Dirichlet BC residual contribution

      computes the surface integral term resulting from integration
      by parts for a particular dof:

      int(n \cdot flux * phi) + int(\sigma (u-g) * phi)
      
  */
template<typename EvalT, typename Traits>
class WeakDirichletResidual
  :
  public panzer::EvaluatorWithBaseImpl<Traits>,
  public PHX::EvaluatorDerived<EvalT, Traits>
{
  public:

    WeakDirichletResidual(
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
  PHX::MDField<ScalarT> normal_dot_flux_plus_pen;
  PHX::MDField<const ScalarT> flux; // i.e., -K \nabla u
  PHX::MDField<const ScalarT> normal;
  PHX::MDField<const ScalarT> sigma;
  PHX::MDField<const ScalarT> dof;
  PHX::MDField<const ScalarT> value;

  std::string basis_name;
  std::size_t basis_index;
  std::size_t num_ip;
  std::size_t num_dim;

}; // end of class WeakDirichletResidual


}

#endif
