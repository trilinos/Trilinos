// @HEADER
// *****************************************************************************
//        Phalanx: A Partial Differential Equation Field Evaluation 
//       Kernel for Flexible Management of Complex Dependency Chains
//
// Copyright 2008 NTESS and the Phalanx contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef PHX_SIMPLE_EVALUATOR_HPP
#define PHX_SIMPLE_EVALUATOR_HPP

#include "Phalanx_Evaluator_WithBaseImpl.hpp"
#include "Phalanx_Evaluator_Derived.hpp"
#include "Phalanx_MDField.hpp"

/** \brief Simple test evaluator

    This class is used to test the EvaluatorUnitTester object.
*/
template<typename EvalT, typename Traits>
class SimpleEvaluator : public PHX::EvaluatorWithBaseImpl<Traits>,
                        public PHX::EvaluatorDerived<EvalT, Traits>  {

  using ScalarT = typename EvalT::ScalarT;
  PHX::MDField<ScalarT,CELL,QP> a;
  PHX::MDField<const ScalarT,CELL,QP> b;
  PHX::MDField<const ScalarT,CELL,QP,DIM> c;
  
public:
  SimpleEvaluator(const Teuchos::RCP<PHX::DataLayout>& a_layout,
                  const Teuchos::RCP<PHX::DataLayout>& b_layout,
                  const Teuchos::RCP<PHX::DataLayout>& c_layout);
  void evaluateFields(typename Traits::EvalData workset) override;
  KOKKOS_INLINE_FUNCTION
  void operator () (const Kokkos::TeamPolicy<PHX::exec_space>::member_type& team) const;
};

#include "SimpleEvaluator_Def.hpp"

#endif
