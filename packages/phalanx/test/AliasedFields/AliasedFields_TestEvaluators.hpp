// @HEADER
// *****************************************************************************
//        Phalanx: A Partial Differential Equation Field Evaluation 
//       Kernel for Flexible Management of Complex Dependency Chains
//
// Copyright 2008 NTESS and the Phalanx contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef PHX_ALIASED_FIELDS_TEST_EVALUATORS_HPP
#define PHX_ALIASED_FIELDS_TEST_EVALUATORS_HPP

#include "Phalanx_Evaluator_WithBaseImpl.hpp"
#include "Phalanx_Evaluator_Derived.hpp"
#include "Phalanx_MDField.hpp"

template<typename EvalT, typename Traits>
class EvalA : public PHX::EvaluatorWithBaseImpl<Traits>,
              public PHX::EvaluatorDerived<EvalT, Traits>  {

  typedef typename EvalT::ScalarT ScalarT;
  PHX::MDField<ScalarT,CELL,BASIS> a;
  PHX::MDField<const ScalarT,CELL,BASIS> b;

public:
  EvalA(const Teuchos::RCP<PHX::DataLayout> & data_layout)
  {
    this->setName("Eval A");
    a = PHX::MDField<ScalarT,CELL,BASIS>("a",data_layout);
    this->addEvaluatedField(a);
    b = PHX::MDField<const ScalarT,CELL,BASIS>("b",data_layout);
    this->addDependentField(b);
  }

  void postRegistrationSetup(typename Traits::SetupData /* d */,
                             PHX::FieldManager<Traits>& /* vm */) override
  {}
  
  void evaluateFields(typename Traits::EvalData ) override
  {
    a.deep_copy(b);
  }
};


template<typename EvalT, typename Traits>
class EvalC : public PHX::EvaluatorWithBaseImpl<Traits>,
              public PHX::EvaluatorDerived<EvalT, Traits>  {

  typedef typename EvalT::ScalarT ScalarT;
  PHX::MDField<ScalarT,CELL,BASIS> c;

public:
  EvalC(const Teuchos::RCP<PHX::DataLayout> & data_layout)
  {
    this->setName("Eval C");
    c = PHX::MDField<ScalarT,CELL,BASIS>("c",data_layout);
    this->addEvaluatedField(c);
  }

  void postRegistrationSetup(typename Traits::SetupData /* d */,
                             PHX::FieldManager<Traits>& /* vm */) override
  {}
  
  void evaluateFields(typename Traits::EvalData ) override
  {
    c.deep_copy(4.0);
  }
};


#endif
