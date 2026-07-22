// @HEADER
// *****************************************************************************
//        Phalanx: A Partial Differential Equation Field Evaluation 
//       Kernel for Flexible Management of Complex Dependency Chains
//
// Copyright 2008 NTESS and the Phalanx contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef PHX_MOCK_DAG_EVALUATOR_HPP
#define PHX_MOCK_DAG_EVALUATOR_HPP

#include "Phalanx_Evaluator_WithBaseImpl.hpp"
#include "Phalanx_Evaluator_Derived.hpp"

namespace PHX {

  //! Mock evaluator to checking DAG construction
  template<typename EvalT, typename Traits>
  class MockDAG : public PHX::EvaluatorWithBaseImpl<Traits>,
		  public PHX::EvaluatorDerived<EvalT, Traits>  {
  public:
    void postRegistrationSetup(typename Traits::SetupData d,
			       PHX::FieldManager<Traits>& fm);
    void evaluateFields(typename Traits::EvalData d);
    void evaluates(const std::string& field_name, const bool use_dynamic_layout=false);
    void depends(const std::string& field_name, const bool use_dynamic_layout=false);
    void contributes(const std::string& field_name, const bool use_dynamic_layout=false);
    void unshared(const std::string& field_name);
  };

}

#include "Evaluator_MockDAG_Def.hpp"

#endif
