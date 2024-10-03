// @HEADER
// *****************************************************************************
//        Phalanx: A Partial Differential Equation Field Evaluation 
//       Kernel for Flexible Management of Complex Dependency Chains
//
// Copyright 2008 NTESS and the Phalanx contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef PHX_EXAMPLE_VP_NONLINEAR_SOURCE_HPP
#define PHX_EXAMPLE_VP_NONLINEAR_SOURCE_HPP

#include "Phalanx_config.hpp"
#ifdef  PHX_ENABLE_KOKKOS_AMT
#include "Phalanx_Evaluator_TaskBase.hpp"
#else
#include "Phalanx_Evaluator_WithBaseImpl.hpp"
#endif
#include "Phalanx_Evaluator_Derived.hpp"
#include "Phalanx_MDField.hpp"

template<typename EvalT, typename Traits>
class NonlinearSource :
#ifdef PHX_ENABLE_KOKKOS_AMT
  public PHX::TaskBase<Traits,NonlinearSource<EvalT,Traits>>,
#else
  public PHX::EvaluatorWithBaseImpl<Traits>,
#endif
  public PHX::EvaluatorDerived<EvalT, Traits> {
  
public:
  
  NonlinearSource(const Teuchos::ParameterList& p);
  
  void postRegistrationSetup(typename Traits::SetupData d,
			     PHX::FieldManager<Traits>& vm);
  
  void evaluateFields(typename Traits::EvalData d);
  
  void preEvaluate(typename Traits::PreEvalData d);
  
  void postEvaluate(typename Traits::PostEvalData d);
  
  KOKKOS_INLINE_FUNCTION
  void operator () (const int i) const;

private:
  
  typedef typename EvalT::ScalarT ScalarT;

  PHX::MDField<ScalarT,Cell,Point> source;
  PHX::MDField<const ScalarT,Cell,Point> density;
  PHX::MDField<const ScalarT,Cell,Point> temp;

  std::size_t cell_data_size;
};

#include "Evaluator_NonlinearSource_Def.hpp"

#endif
