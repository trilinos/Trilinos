// @HEADER
// *****************************************************************************
//        Phalanx: A Partial Differential Equation Field Evaluation 
//       Kernel for Flexible Management of Complex Dependency Chains
//
// Copyright 2008 NTESS and the Phalanx contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef PHX_EXAMPLE_VP_CONSTANT_HPP
#define PHX_EXAMPLE_VP_CONSTANT_HPP

#include "Phalanx_config.hpp"
#ifdef  PHX_ENABLE_KOKKOS_AMT
#include "Phalanx_Evaluator_TaskBase.hpp"
#else
#include "Phalanx_Evaluator_WithBaseImpl.hpp"
#endif
#include "Phalanx_Evaluator_Derived.hpp"
#include "Phalanx_MDField.hpp"

#include "Dimension.hpp"

// Dummy evalauator for testing
template<typename EvalT, typename Traits>
class Constant : 
#ifdef PHX_ENABLE_KOKKOS_AMT
  public PHX::TaskBase<Traits,Constant<EvalT,Traits>>,
#else
  public PHX::EvaluatorWithBaseImpl<Traits>,
#endif
  public PHX::EvaluatorDerived<EvalT, Traits> {
  
public:
  
  Constant(Teuchos::ParameterList& p);
  
  void postRegistrationSetup(typename Traits::SetupData d,
			     PHX::FieldManager<Traits>& vm);
  
  void evaluateFields(typename Traits::EvalData ud);
  
  KOKKOS_INLINE_FUNCTION
  void operator () (const int i) const;  

private:
  
  typedef typename EvalT::ScalarT ScalarT;

  ScalarT value;

  PHX::MDField<ScalarT,Cell,Point> constant;

  int num_points;
  
  //! Not neede for problem, but included for some unit testing
  std::size_t dummy_workset_size;
};

#include "Evaluator_Constant_Def.hpp"

#endif
