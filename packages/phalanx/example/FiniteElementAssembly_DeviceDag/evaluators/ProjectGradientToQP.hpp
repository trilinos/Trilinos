// @HEADER
// *****************************************************************************
//        Phalanx: A Partial Differential Equation Field Evaluation 
//       Kernel for Flexible Management of Complex Dependency Chains
//
// Copyright 2008 NTESS and the Phalanx contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef PHX_PROJECT_GRADIENT_TO_QP_HPP
#define PHX_PROJECT_GRADIENT_TO_QP_HPP

#include "Phalanx_Evaluator_WithBaseImpl.hpp"
#include "Phalanx_Evaluator_Derived.hpp"
#include "Phalanx_DeviceEvaluator.hpp"
#include "Phalanx_MDField.hpp"
#include "Dimension.hpp"

template<typename EvalT, typename Traits>
class ProjectGradientToQP : public PHX::EvaluatorWithBaseImpl<Traits>,
			    public PHX::EvaluatorDerived<PHX::MyTraits::Residual, Traits>  {

  using ScalarT = typename EvalT::ScalarT;
  PHX::MDField<const ScalarT,CELL,BASIS> field_at_basis;
  PHX::MDField<ScalarT,CELL,QP,DIM> grad_field_at_qp;

public:

  struct MyDevEvalResidual : public PHX::DeviceEvaluator<Traits> {
    PHX::View<const ScalarT**> field_at_basis;
    PHX::AtomicView<ScalarT***> grad_field_at_qp;
    KOKKOS_FUNCTION MyDevEvalResidual(const PHX::View<const ScalarT**>& in_field_at_basis,
				      const PHX::View<ScalarT***>& in_grad_field_at_qp) :
      field_at_basis(in_field_at_basis), grad_field_at_qp(in_grad_field_at_qp) {}
    KOKKOS_FUNCTION MyDevEvalResidual(const MyDevEvalResidual& src) = default;
    KOKKOS_FUNCTION void evaluate(const typename PHX::DeviceEvaluator<Traits>::member_type& team,
                                  typename Traits::EvalData workset) override;
  };

  struct MyDevEvalJacobian : public PHX::DeviceEvaluator<Traits> {
    PHX::View<const ScalarT**> field_at_basis;
    PHX::View<ScalarT***> grad_field_at_qp;
    KOKKOS_FUNCTION MyDevEvalJacobian(const PHX::View<const ScalarT**>& in_field_at_basis,
				      const PHX::View<ScalarT***>& in_grad_field_at_qp) :
      field_at_basis(in_field_at_basis), grad_field_at_qp(in_grad_field_at_qp) {}
    KOKKOS_FUNCTION MyDevEvalJacobian(const MyDevEvalJacobian& src) = default;
    KOKKOS_FUNCTION void evaluate(const typename PHX::DeviceEvaluator<Traits>::member_type& team,
                                  typename Traits::EvalData workset) override;
  };
  
  ProjectGradientToQP(const std::string& field_name,
                      const Teuchos::RCP<PHX::DataLayout>& basis_layout,
                      const Teuchos::RCP<PHX::DataLayout>& grad_qp_layout);
  PHX::DeviceEvaluator<Traits>* createDeviceEvaluator() const override;
  void evaluateFields(typename Traits::EvalData workset) override;
};

#endif
