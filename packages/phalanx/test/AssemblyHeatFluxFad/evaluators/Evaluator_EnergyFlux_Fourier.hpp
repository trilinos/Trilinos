// @HEADER
// *****************************************************************************
//        Phalanx: A Partial Differential Equation Field Evaluation 
//       Kernel for Flexible Management of Complex Dependency Chains
//
// Copyright 2008 NTESS and the Phalanx contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef PHX_EXAMPLE_VP_FOURIER_HPP
#define PHX_EXAMPLE_VP_FOURIER_HPP

#include "Phalanx_config.hpp"
#ifdef  PHX_ENABLE_KOKKOS_AMT
#include "Phalanx_Evaluator_TaskBase.hpp"
#else
#include "Phalanx_Evaluator_WithBaseImpl.hpp"
#endif
#include "Phalanx_Evaluator_Derived.hpp"
#include "Phalanx_MDField.hpp"

template<typename EvalT, typename Traits>
class Fourier : 
#ifdef PHX_ENABLE_KOKKOS_AMT
  public PHX::TaskBase<Traits,Fourier<EvalT,Traits>>,
#else
  public PHX::EvaluatorWithBaseImpl<Traits>,
#endif
  public PHX::EvaluatorDerived<EvalT, Traits> {
  
public:
  
  Fourier(const Teuchos::ParameterList& p);
  
  void postRegistrationSetup(typename Traits::SetupData d,
			     PHX::FieldManager<Traits>& vm);
  
  void evaluateFields(typename Traits::EvalData d);

  KOKKOS_INLINE_FUNCTION
  void operator () (const int i) const;
  
private:
  
  typedef typename EvalT::ScalarT ScalarT;

  PHX::MDField<ScalarT,Cell,QuadPoint,Dim> flux;
  PHX::MDField<const ScalarT,Cell,QuadPoint> density;
  PHX::MDField<const ScalarT,Cell,QuadPoint> dc;
  PHX::MDField<const ScalarT,Cell,QuadPoint,Dim> grad_temp;
  
  PHX::index_size_type num_qp;
  PHX::index_size_type num_dim;

};

#include "Evaluator_EnergyFlux_Fourier_Def.hpp"

#endif
