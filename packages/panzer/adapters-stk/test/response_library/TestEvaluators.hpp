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
#include "Panzer_Evaluator_Macros.hpp"

namespace panzer {

template<typename EvalT, typename Traits>
class TestEvaluator
  :
  public panzer::EvaluatorWithBaseImpl<Traits>,
  public PHX::EvaluatorDerived<EvalT, Traits>
{
  public:

    TestEvaluator(
      const Teuchos::ParameterList& p);

    void
    evaluateFields(
      typename Traits::EvalData d);

  private:

    using ScalarT = typename EvalT::ScalarT;

  PHX::MDField<ScalarT,Cell> dogValues;
  PHX::MDField<ScalarT,Cell> hrsValues;

}; // end of class TestEvaluator


//**********************************************************************
template<typename EvalT, typename Traits>
TestEvaluator<EvalT, Traits>::
TestEvaluator(
  const Teuchos::ParameterList& p)
{
  // Read from parameters
  int worksetSize = p.get<int>("Workset Size");

  Teuchos::RCP<PHX::DataLayout> dl = Teuchos::rcp(new PHX::MDALayout<Cell>(worksetSize));
  // grab information from quadrature rule
  dogValues = PHX::MDField<ScalarT,Cell>("Dog", dl);
  hrsValues = PHX::MDField<ScalarT,Cell>("Horse", dl);

  this->addEvaluatedField(dogValues);
  this->addEvaluatedField(hrsValues);
  
  std::string n = "TestEvaluator";
  this->setName(n);
}

//**********************************************************************
template<typename EvalT, typename Traits>
void
TestEvaluator<EvalT, Traits>::
evaluateFields(
  typename Traits::EvalData workset)
{ 
   double extra = 0.0;
   if(this->wda(workset).block_id=="block_1")
      extra = 44.3;

   for(int i=0;i<dogValues.extent_int(0);i++) {
      dogValues(i) = double(i) + 1.0 + extra;
      hrsValues(i) = -double(i) - 5.5 + extra;
   }
}

//**********************************************************************

}

#endif
