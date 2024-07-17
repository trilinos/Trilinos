// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef PANZER_EVALUATOR_SCALAR_DECL_HPP
#define PANZER_EVALUATOR_SCALAR_DECL_HPP

#include <string>
#include "Panzer_Dimension.hpp"
#include "Phalanx_Evaluator_Macros.hpp"
#include "Phalanx_MDField.hpp"
#include "Kokkos_DynRankView.hpp"

#include "Panzer_Evaluator_Macros.hpp"

namespace panzer {
    
/** This integrates a scalar quanity over each cell.
  * It is useful for comptuing integral responses.

  \verbatim
    <ParameterList>
      <Parameter name="Integral Name" type="string" value="<Name to give to the integral field>"/>
      <Parameter name="Integrand Name" type="string" value="<Name of integrand>"/>
      <Parameter name="IR" type="RCP<IntegrationRule>" value="<user specified IntegrationRule>"/>
      <Parameter name="Multiplier" type="double" value="<Scaling factor, default=1>"/>
      <Parameter name="Field Multipliers" type="RCP<const vector<string> >" value="<Other scalar multiplier fields>"/>
    </ParameterList>
  \endverbatim
  */
template<typename EvalT, typename Traits>
class Integrator_Scalar
  :
  public panzer::EvaluatorWithBaseImpl<Traits>,
  public PHX::EvaluatorDerived<EvalT, Traits>
{
  public:

    Integrator_Scalar(
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
  
  PHX::MDField<ScalarT> integral;  // result
    
  PHX::MDField<const ScalarT,Cell,IP> scalar; // function to be integrated

  typename PHX::View<PHX::UnmanagedView<const ScalarT**>* >::HostMirror field_multipliers_h;
  PHX::View<PHX::UnmanagedView<const ScalarT**>* > field_multipliers;

  std::size_t num_qp;
  std::size_t quad_index;
  int quad_order;

  double multiplier;

  Kokkos::DynRankView<ScalarT,typename PHX::DevLayout<ScalarT>::type,PHX::Device> tmp;

public:
  // for testing purposes
  const PHX::FieldTag & getFieldTag() const 
  { return integral.fieldTag(); }

private:
  Teuchos::RCP<Teuchos::ParameterList> getValidParameters() const;

}; // end of class Integrator_Scalar


}

#endif
