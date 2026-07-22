// @HEADER
// *****************************************************************************
//        Phalanx: A Partial Differential Equation Field Evaluation 
//       Kernel for Flexible Management of Complex Dependency Chains
//
// Copyright 2008 NTESS and the Phalanx contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef PHX_DUPLICATE_FIELD_EVALUATOR_HPP
#define PHX_DUPLICATE_FIELD_EVALUATOR_HPP

#include "Phalanx_Evaluator_WithBaseImpl.hpp"
#include "Phalanx_Evaluator_Derived.hpp"
#include "Phalanx_MDField.hpp"

/** \brief Corner case test evaluator 

    This class is used to test the EvaluatorUnitTester object. It also
    servers to test a corner case for the default implementation of
    the evaluator class. In particular, if a single evaluator
    registers two MDField objects that point to the same underlying
    MDField, both fields need to be bound to the same memory with
    automatic binding. This replicates a bug from a using an
    unordered_map instead of unordered_multimap to store the field
    pointers in the EvaluatorWithBaseImpl class. This is a valid use
    case where some fields in an evaluator have arbitrary names,
    chosen at runtime, and could potentially require the same field as
    another field in the evaluator.
    */
template<typename EvalT, typename Traits>
class DuplicateFieldEvaluator : public PHX::EvaluatorWithBaseImpl<Traits>,
                                public PHX::EvaluatorDerived<EvalT, Traits>  {

  using ScalarT = typename EvalT::ScalarT;
  PHX::MDField<ScalarT,CELL,QP> a;
  PHX::MDField<const ScalarT,CELL,QP> b1;
  PHX::MDField<const ScalarT,CELL,QP> b2; // Points to same memory as b1
  PHX::MDField<const ScalarT,CELL,QP,DIM> c;
  
public:
  DuplicateFieldEvaluator(const Teuchos::RCP<PHX::DataLayout>& a_layout,
                          const Teuchos::RCP<PHX::DataLayout>& b_layout,
                          const Teuchos::RCP<PHX::DataLayout>& c_layout);
  void evaluateFields(typename Traits::EvalData workset) override;
  KOKKOS_INLINE_FUNCTION
  void operator () (const Kokkos::TeamPolicy<PHX::exec_space>::member_type& team) const;
};

#include "DuplicateFieldEvaluator_Def.hpp"

#endif
