// @HEADER
// *****************************************************************************
//        Phalanx: A Partial Differential Equation Field Evaluation 
//       Kernel for Flexible Management of Complex Dependency Chains
//
// Copyright 2008 NTESS and the Phalanx contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef PHX_INTEGRATE_SOURCE_TERM_HPP
#define PHX_INTEGRATE_SOURCE_TERM_HPP

#include "Phalanx_Evaluator_WithBaseImpl.hpp"
#include "Phalanx_Evaluator_Derived.hpp"
#include "Phalanx_DeviceEvaluator.hpp"
#include "Phalanx_MDField.hpp"
#include "Dimension.hpp"

template<typename EvalT, typename Traits>
class IntegrateSourceTerm : public PHX::EvaluatorWithBaseImpl<Traits>,
                            public PHX::EvaluatorDerived<EvalT, Traits>  {

  using ScalarT = typename EvalT::ScalarT;
  PHX::MDField<const ScalarT,CELL,QP> source;
  PHX::MDField<ScalarT,CELL,BASIS> residual;

public:

  struct MyDevEvalResidual : public PHX::DeviceEvaluator<Traits> {
    PHX::View<const ScalarT**> source;
    PHX::AtomicView<ScalarT**> residual;
    KOKKOS_FUNCTION MyDevEvalResidual(const PHX::View<const ScalarT**>& in_source,
			      const PHX::View<ScalarT**>& in_residual) :
      source(in_source), residual(in_residual) {}
    KOKKOS_FUNCTION MyDevEvalResidual(const MyDevEvalResidual& src) = default;
    KOKKOS_FUNCTION void evaluate(const typename PHX::DeviceEvaluator<Traits>::member_type& team,
                                  const typename Traits::EvalData data) override;
  };

  struct MyDevEvalJacobian : public PHX::DeviceEvaluator<Traits> {
    PHX::View<const ScalarT**> source;
#ifdef PHX_ENABLE_KOKKOS_AMT
    // Make residual atomic so that AMT mode can sum diffusion and source terms at same time
    Kokkos::View<ScalarT**,typename PHX::DevLayout<ScalarT>::type,PHX::Device,Kokkos::MemoryTraits<Kokkos::Atomic>> residual;
#else
    PHX::View<ScalarT**> residual;
#endif
    KOKKOS_FUNCTION MyDevEvalJacobian(const PHX::View<const ScalarT**>& in_source,
			      const PHX::View<ScalarT**>& in_residual) :
      source(in_source), residual(in_residual) {}
    KOKKOS_FUNCTION MyDevEvalJacobian(const MyDevEvalJacobian& src) = default;
    KOKKOS_FUNCTION void evaluate(const typename PHX::DeviceEvaluator<Traits>::member_type& team,
                                  const typename Traits::EvalData data) override;
  };
  
  IntegrateSourceTerm(const std::string& source_name,
                      const Teuchos::RCP<PHX::DataLayout>& source_layout,
                      const std::string& residual_name,
                      const Teuchos::RCP<PHX::DataLayout>& residual_layout);
  PHX::DeviceEvaluator<Traits>* createDeviceEvaluator() const override;
  void evaluateFields(typename Traits::EvalData workset) override;
};

#endif
