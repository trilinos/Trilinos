// @HEADER
// *****************************************************************************
//        Phalanx: A Partial Differential Equation Field Evaluation 
//       Kernel for Flexible Management of Complex Dependency Chains
//
// Copyright 2008 NTESS and the Phalanx contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef PHX_EVALUATOR_UNAMANGED_FIELD_DUMMY_HPP
#define PHX_EVALUATOR_UNAMANGED_FIELD_DUMMY_HPP

#include "Phalanx_Evaluator_WithBaseImpl.hpp"
#include "Phalanx_Evaluator_Derived.hpp"
#include "Phalanx_MDField.hpp"

namespace PHX {

/** \brief Evaluator that performs no computations. Typically used to
    satisfy DAG dependencies for unmanaged fields that are evalatued
    external to the DAG. */
template<typename EvalT, typename Traits, typename FieldT>
class UnmanagedFieldDummy : public PHX::EvaluatorWithBaseImpl<Traits>,
                            public PHX::EvaluatorDerived<EvalT, Traits>  {
public:
  UnmanagedFieldDummy(const FieldT& f)
  {
    this->addEvaluatedField(f.fieldTag());
    this->setName("UnmanagedFieldDummy");
  }
  void evaluateFields(typename Traits::EvalData /* workset */) override {}
};

}

#endif
