// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef PANZER_EVALUATOR_INTERFACE_RESIDUAL_DECL_HPP
#define PANZER_EVALUATOR_INTERFACE_RESIDUAL_DECL_HPP

#include "Phalanx_Evaluator_Macros.hpp"
#include "Phalanx_MDField.hpp"

#include "Panzer_Evaluator_Macros.hpp"

namespace panzer {
    
  /** \brief Evaluates an interface BC residual contribution

      computes the surface integral term resulting from integration
      by parts for a particular dof:

      int(n \cdot (flux * phi) )
  */
template<typename EvalT, typename Traits>
class InterfaceResidual
  :
  public panzer::EvaluatorWithBaseImpl<Traits>,
  public PHX::EvaluatorDerived<EvalT, Traits>
{
  public:

    InterfaceResidual(
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
  PHX::MDField<ScalarT> normal_dot_flux;
  PHX::MDField<const ScalarT> flux;
  PHX::MDField<const ScalarT> normal;

  std::string basis_name;
  std::size_t basis_index;
  std::size_t num_ip;
  std::size_t num_dim;

}; // end of class InterfaceResidual


}

#endif
