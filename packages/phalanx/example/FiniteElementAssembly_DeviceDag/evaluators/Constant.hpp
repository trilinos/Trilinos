// @HEADER
// *****************************************************************************
//        Phalanx: A Partial Differential Equation Field Evaluation 
//       Kernel for Flexible Management of Complex Dependency Chains
//
// Copyright 2008 NTESS and the Phalanx contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef PHX_EXAMPLE_EVALUATOR_CONSTANT_HPP
#define PHX_EXAMPLE_EVALUATOR_CONSTANT_HPP

#include "Phalanx_Evaluator_WithBaseImpl.hpp"
#include "Phalanx_Evaluator_Derived.hpp"
#include "Phalanx_DeviceEvaluator.hpp"
#include "Phalanx_MDField.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Dimension.hpp"

template<typename EvalT, typename Traits>
class Constant : public PHX::EvaluatorWithBaseImpl<Traits>,
                 public PHX::EvaluatorDerived<EvalT, Traits>  {

  using ScalarT = typename EvalT::ScalarT;
  ScalarT value;
  PHX::MDField<ScalarT,CELL,POINT> constant;
  
public:
  struct MyDevEval : public PHX::DeviceEvaluator<Traits> {
    KOKKOS_FUNCTION
    void evaluate(const typename PHX::DeviceEvaluator<Traits>::member_type& ,
                  const typename Traits::EvalData ) override {}
  };

  Constant(const std::string& field_name,
           const Teuchos::RCP<PHX::DataLayout>& layout,
           const double& value);
  PHX::DeviceEvaluator<Traits>* createDeviceEvaluator() const override;
  void postRegistrationSetup(typename Traits::SetupData d,
                             PHX::FieldManager<Traits>& fm) override;
  void evaluateFields(typename Traits::EvalData d) override;

};

#endif
