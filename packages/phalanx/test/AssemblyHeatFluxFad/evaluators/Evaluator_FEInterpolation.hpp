// @HEADER
// *****************************************************************************
//        Phalanx: A Partial Differential Equation Field Evaluation 
//       Kernel for Flexible Management of Complex Dependency Chains
//
// Copyright 2008 NTESS and the Phalanx contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef PHX_EXAMPLE_VP_FE_INTERPOLATION_HPP
#define PHX_EXAMPLE_VP_FE_INTERPOLATION_HPP

#include "Phalanx_config.hpp"
#ifdef  PHX_ENABLE_KOKKOS_AMT
#include "Phalanx_Evaluator_TaskBase.hpp"
#else
#include "Phalanx_Evaluator_WithBaseImpl.hpp"
#endif
#include "Phalanx_Evaluator_Derived.hpp"
#include "Phalanx_MDField.hpp"
/** \brief Finite Element Interpolation Evaluator

    This object evaluates a scalar field and it's gradient at the
    quadrature points for a specific variable.

*/
template<typename EvalT, typename Traits>
class FEInterpolation :
#ifdef PHX_ENABLE_KOKKOS_AMT
  public PHX::TaskBase<Traits,FEInterpolation<EvalT,Traits>>,
#else
  public PHX::EvaluatorWithBaseImpl<Traits>,
#endif
  public PHX::EvaluatorDerived<EvalT, Traits> {
  
public:
  
  FEInterpolation(const Teuchos::ParameterList& p);
  
  void postRegistrationSetup(typename Traits::SetupData d,
			     PHX::FieldManager<Traits>& vm);
  
  void evaluateFields(typename Traits::EvalData d);

  KOKKOS_INLINE_FUNCTION
  void operator () (const int i) const;

#ifdef PHX_ENABLE_KOKKOS_AMT
  Kokkos::Future<void,PHX::exec_space>
    createTask(Kokkos::TaskScheduler<PHX::exec_space>& policy,
	       const int& work_size,
               const std::vector<Kokkos::Future<void,PHX::exec_space>>& dependent_futures,
	       typename Traits::EvalData d) override;
#endif
  
private:

  typedef typename EvalT::ScalarT ScalarT;

  //! Values at nodes
  PHX::MDField<const ScalarT,Cell,Node> val_node;

  //! Values at quadrature points
  PHX::MDField<ScalarT,Cell,QuadPoint> val_qp;

  //! Gradient values at quadrature points
  PHX::MDField<ScalarT,Cell,QuadPoint,Dim> val_grad_qp;
   
  PHX::index_size_type num_nodes;
  PHX::index_size_type num_qp;
  PHX::index_size_type num_dim;

  // dummy field for unit testing mixed scalar types (double field in
  // Jacobian evaluation
  PHX::MDField<double,Cell,QuadPoint,Dim> dummy;
 
  Kokkos::View<double**,PHX::Device> phi;
  Kokkos::View<double***,PHX::Device> grad_phi;

};

#include "Evaluator_FEInterpolation_Def.hpp"

#endif
