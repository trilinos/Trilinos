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
#include "Phalanx_MDField.hpp"
#include "Dimension.hpp"

//! Project field values from basis to qp.
template<typename EvalT, typename Traits>
class ProjectGradientToQP : public PHX::EvaluatorWithBaseImpl<Traits>,
                            public PHX::EvaluatorDerived<EvalT, Traits>  {

  using ScalarT = typename EvalT::ScalarT;
  PHX::MDField<const ScalarT,CELL,BASIS> field_at_basis;
  PHX::MDField<ScalarT,CELL,QP,DIM> grad_field_at_qp;
  Kokkos::View<double****,PHX::Device> grad_basis_view;
  
public:
  ProjectGradientToQP(const std::string& field_name,
                      const Teuchos::RCP<PHX::DataLayout>& basis_layout,
                      const Teuchos::RCP<PHX::DataLayout>& grad_qp_layout);
  void evaluateFields(typename Traits::EvalData workset) override;
  KOKKOS_INLINE_FUNCTION
  void operator () (const Kokkos::TeamPolicy<PHX::exec_space>::member_type& team) const;
};

#endif
