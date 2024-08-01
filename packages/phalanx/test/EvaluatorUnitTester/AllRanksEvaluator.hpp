// @HEADER
// *****************************************************************************
//        Phalanx: A Partial Differential Equation Field Evaluation 
//       Kernel for Flexible Management of Complex Dependency Chains
//
// Copyright 2008 NTESS and the Phalanx contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef PHX_ALL_RANKS_EVALUATOR_HPP
#define PHX_ALL_RANKS_EVALUATOR_HPP

#include "Phalanx_Evaluator_WithBaseImpl.hpp"
#include "Phalanx_Evaluator_Derived.hpp"
#include "Phalanx_MDField.hpp"

/** \brief Evaluator to unit test the EvaluatorUnitTester for field checks for all ranks

    This class is used to test the EvaluatorUnitTester object.
*/
template<typename EvalT, typename Traits>
class AllRanksEvaluator : public PHX::EvaluatorWithBaseImpl<Traits>,
                          public PHX::EvaluatorDerived<EvalT, Traits>  {

  using ScalarT = typename EvalT::ScalarT;
  PHX::MDField<const ScalarT,R> f1;
  PHX::MDField<const ScalarT,R,R> f2;
  PHX::MDField<const ScalarT,R,R,R> f3;
  PHX::MDField<const ScalarT,R,R,R,R> f4;
  PHX::MDField<const ScalarT,R,R,R,R,R> f5;
  PHX::MDField<const ScalarT,R,R,R,R,R,R> f6;
  PHX::MDField<ScalarT,R> x1;
  PHX::MDField<ScalarT,R,R> x2;
  PHX::MDField<ScalarT,R,R,R> x3;
  PHX::MDField<ScalarT,R,R,R,R> x4;
  PHX::MDField<ScalarT,R,R,R,R,R> x5;
  PHX::MDField<ScalarT,R,R,R,R,R,R> x6;
  
public:
  AllRanksEvaluator(const Teuchos::RCP<PHX::DataLayout>& dl1,
                    const Teuchos::RCP<PHX::DataLayout>& dl2,
                    const Teuchos::RCP<PHX::DataLayout>& dl3,
                    const Teuchos::RCP<PHX::DataLayout>& dl4,
                    const Teuchos::RCP<PHX::DataLayout>& dl5,
                    const Teuchos::RCP<PHX::DataLayout>& dl6);
  void evaluateFields(typename Traits::EvalData workset) override;
  KOKKOS_INLINE_FUNCTION
  void operator () (const Kokkos::TeamPolicy<PHX::exec_space>::member_type& team) const;
};

#include "AllRanksEvaluator_Def.hpp"

#endif
