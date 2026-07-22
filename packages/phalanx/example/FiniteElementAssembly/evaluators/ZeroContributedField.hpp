// @HEADER
// *****************************************************************************
//        Phalanx: A Partial Differential Equation Field Evaluation 
//       Kernel for Flexible Management of Complex Dependency Chains
//
// Copyright 2008 NTESS and the Phalanx contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef PHX_ZERO_CONTRIBUTED_FIELD_HPP
#define PHX_ZERO_CONTRIBUTED_FIELD_HPP

#include "Phalanx_Evaluator_WithBaseImpl.hpp"
#include "Phalanx_Evaluator_Derived.hpp"
#include "Phalanx_MDField.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Dimension.hpp"

template<typename EvalT, typename Traits>
class ZeroContributedField : public PHX::EvaluatorWithBaseImpl<Traits>,
                             public PHX::EvaluatorDerived<EvalT, Traits>  {

  using ScalarT = typename EvalT::ScalarT;
  PHX::MDField<ScalarT> field;
  
public:
  ZeroContributedField(const std::string& field_name,
                       const Teuchos::RCP<PHX::DataLayout>& layout);
  void evaluateFields(typename Traits::EvalData d) override;
};

#endif
