// @HEADER
// *****************************************************************************
//        Phalanx: A Partial Differential Equation Field Evaluation 
//       Kernel for Flexible Management of Complex Dependency Chains
//
// Copyright 2008 NTESS and the Phalanx contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef PHX_EXAMPLE_VP_DENSITY_HPP
#define PHX_EXAMPLE_VP_DENSITY_HPP

#include "Phalanx_config.hpp"
#ifdef  PHX_ENABLE_KOKKOS_AMT
#include "Phalanx_Evaluator_TaskBase.hpp"
#else
#include "Phalanx_Evaluator_WithBaseImpl.hpp"
#endif
#include "Phalanx_Evaluator_Derived.hpp"
#include "Phalanx_MDField.hpp"

//struct DensityTag {};

template<typename EvalT, typename Traits>
class Density :
#ifdef PHX_ENABLE_KOKKOS_AMT
  public PHX::TaskBase<Traits,Density<EvalT,Traits>>,
#else
  public PHX::EvaluatorWithBaseImpl<Traits>,
#endif
  public PHX::EvaluatorDerived<EvalT, Traits> {

public:

  Density(const Teuchos::ParameterList& p);

  void postRegistrationSetup(typename Traits::SetupData d,
			     PHX::FieldManager<Traits>& vm);

  void evaluateFields(typename Traits::EvalData ud);

  KOKKOS_INLINE_FUNCTION
    void operator () (const int i) const;

  // KOKKOS_INLINE_FUNCTION
  //   void operator () (const DensityTag, const int i) const;

  // KOKKOS_INLINE_FUNCTION
  //   void operator () (const DensityTag, typename Kokkos::TeamPolicy<>::member_type & team) const;

private:

  typedef typename EvalT::ScalarT ScalarT;

  double constant;

  PHX::MDField<ScalarT,Cell,Point> density;

  // Should normally use a const dependent field. We will use a
  // non-const to test the method addNonConstDependentField.
  // PHX::MDField<const ScalarT,Cell,Point> temp;
  PHX::MDField<ScalarT,Cell,Point> temp;

  std::size_t cell_data_size;

};

#include "Evaluator_Density_Def.hpp"

#endif
