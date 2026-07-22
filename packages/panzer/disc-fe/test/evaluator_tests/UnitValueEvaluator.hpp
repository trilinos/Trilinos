// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef UNIT_VALUE_EVALUATOR
#define UNIT_VALUE_EVALUATOR

#include <string>
#include "Phalanx_Evaluator_Macros.hpp"
#include "Phalanx_MDField.hpp"
#include "Panzer_IntegrationRule.hpp"

namespace panzer {

template<typename EvalT, typename Traits>
class UnitValueEvaluator
  :
  public PHX::EvaluatorWithBaseImpl<Traits>,
  public PHX::EvaluatorDerived<EvalT, Traits>
{
  public:

    UnitValueEvaluator(
      const Teuchos::ParameterList& p);

    void
    evaluateFields(
      typename Traits::EvalData d);

  private:

    using ScalarT = typename EvalT::ScalarT;

  PHX::MDField<ScalarT,Cell,panzer::IP> unitValue;

}; // end of class UnitValueEvaluator


//**********************************************************************
template<typename EvalT, typename Traits>
UnitValueEvaluator<EvalT, Traits>::
UnitValueEvaluator(
  const Teuchos::ParameterList& p)
{
  // Read from parameters
  const std::string name = p.get<std::string>("Name");
  Teuchos::RCP<panzer::IntegrationRule> ir
     = p.get< Teuchos::RCP<panzer::IntegrationRule> >("IR");

  // grab information from quadrature rule
  unitValue = PHX::MDField<ScalarT,Cell,IP>(name, ir->dl_scalar);

  this->addEvaluatedField(unitValue);
  
  std::string n = "UnitValueEvaluator: " + name;
  this->setName(n);
}

//**********************************************************************
template<typename EvalT, typename Traits>
void
UnitValueEvaluator<EvalT, Traits>::
evaluateFields(
  typename Traits::EvalData  /* workset */)
{ 
  Kokkos::deep_copy(unitValue.get_static_view(), 1.0);
}

//**********************************************************************

}

#endif
