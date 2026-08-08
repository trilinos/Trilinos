// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef PANZER_EVALUATOR_TRANSIENT_BASISTIMESSCALAR_DECL_HPP
#define PANZER_EVALUATOR_TRANSIENT_BASISTIMESSCALAR_DECL_HPP

#include <string>
#include "Panzer_Dimension.hpp"
#include "Phalanx_Evaluator_Macros.hpp"
#include "Phalanx_MDField.hpp"
#include "Kokkos_DynRankView.hpp"

#include "Panzer_Evaluator_Macros.hpp"

namespace panzer {

/** This integrates a scalar quantity times a basis function over each
  * cell, contributing to the residual only when the workset indicates
  * transient (time derivative) terms should be evaluated. It is
  * useful for building the mass-matrix-like terms of a time-dependent
  * weak form.
  */
template<typename EvalT, typename Traits>
class Integrator_TransientBasisTimesScalar
  :
  public panzer::EvaluatorWithBaseImpl<Traits>,
  public PHX::EvaluatorDerived<EvalT, Traits>
{
  public:

    /// \brief Construct from a ParameterList specifying the residual field name, integrand name, basis, integration rule, scaling multiplier, and optional field multipliers.
    Integrator_TransientBasisTimesScalar(
      const Teuchos::ParameterList& p);

    /// \brief Looks up the basis index needed by evaluateFields(). Called once before the first evaluation.
    void
    postRegistrationSetup(
      typename Traits::SetupData d,
      PHX::FieldManager<Traits>& fm);

    /// \brief Computes the residual field for this workset per the class description, when the workset indicates transient terms should be evaluated.
    void
    evaluateFields(
      typename Traits::EvalData d);

  private:

    using ScalarT = typename EvalT::ScalarT;
  
  PHX::MDField<ScalarT,Cell,BASIS> residual;
    
  PHX::MDField<const ScalarT,Cell,IP> scalar;

  std::vector<PHX::MDField<const ScalarT,Cell,IP> > field_multipliers;

  std::size_t num_nodes;

  std::size_t num_qp;

  double multiplier;

  std::string basis_name;
  std::size_t basis_index;

  Kokkos::DynRankView<ScalarT,typename PHX::DevLayout<ScalarT>::type,PHX::Device> tmp;

private:
  /// \brief Returns the ParameterList of valid parameters/defaults for this evaluator's constructor.
  Teuchos::RCP<Teuchos::ParameterList> getValidParameters() const;

}; // end of class Integrator_TransientBasisTimesScalar


}

#endif
