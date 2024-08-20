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
    
template<typename EvalT, typename Traits>
class Integrator_TransientBasisTimesScalar
  :
  public panzer::EvaluatorWithBaseImpl<Traits>,
  public PHX::EvaluatorDerived<EvalT, Traits>
{
  public:

    Integrator_TransientBasisTimesScalar(
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
  Teuchos::RCP<Teuchos::ParameterList> getValidParameters() const;

}; // end of class Integrator_TransientBasisTimesScalar


}

#endif
