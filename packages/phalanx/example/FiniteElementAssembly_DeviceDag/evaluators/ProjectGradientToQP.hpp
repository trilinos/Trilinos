// @HEADER
// ************************************************************************
//
//        Phalanx: A Partial Differential Equation Field Evaluation 
//       Kernel for Flexible Management of Complex Dependency Chains
//                    Copyright 2008 Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Roger Pawlowski (rppawlo@sandia.gov), Sandia
// National Laboratories.
//
// ************************************************************************
// @HEADER

#ifndef PHX_PROJECT_GRADIENT_TO_QP_HPP
#define PHX_PROEJCT_GRADIENT_TO_QP_HPP

#include "Phalanx_Evaluator_WithBaseImpl.hpp"
#include "Phalanx_Evaluator_Derived.hpp"
#include "Phalanx_DeviceEvaluator.hpp"
#include "Phalanx_MDField.hpp"
#include "Dimension.hpp"

//! Project field values from basis to qp.
template<typename EvalT, typename Traits> class ProjectGradientToQP;

// ***********************
// Residual Specialization
// ***********************

template<typename Traits>
class ProjectGradientToQP<PHX::MyTraits::Residual,Traits> : public PHX::EvaluatorWithBaseImpl<Traits>,
							    public PHX::EvaluatorDerived<PHX::MyTraits::Residual, Traits>  {

  using ScalarT = typename PHX::MyTraits::Residual::ScalarT;
  PHX::MDField<const ScalarT,CELL,BASIS> field_at_basis;
  PHX::MDField<ScalarT,CELL,QP,DIM> grad_field_at_qp;

public:
  struct MyDevEval : public PHX::DeviceEvaluator<Traits> {
    PHX::View<const ScalarT**> field_at_basis;
    PHX::View<ScalarT***> grad_field_at_qp;
    KOKKOS_FUNCTION MyDevEval(const PHX::View<const ScalarT**>& in_field_at_basis,
			      const PHX::View<ScalarT***>& in_grad_field_at_qp) :
      field_at_basis(in_field_at_basis), grad_field_at_qp(in_grad_field_at_qp) {}
    KOKKOS_FUNCTION MyDevEval(const MyDevEval& src) = default;
    KOKKOS_FUNCTION void evaluate(const typename PHX::DeviceEvaluator<Traits>::member_type& team,
                                  typename Traits::EvalData workset) override;
  };
  
  ProjectGradientToQP(const std::string& field_name,
                      const Teuchos::RCP<PHX::DataLayout>& basis_layout,
                      const Teuchos::RCP<PHX::DataLayout>& grad_qp_layout);
  PHX::DeviceEvaluator<Traits>* createDeviceEvaluator() const override;
  void evaluateFields(typename Traits::EvalData workset) override;
};

// ***********************
// Jacobian Specialization
// ***********************

template<typename Traits>
class ProjectGradientToQP<PHX::MyTraits::Jacobian,Traits> : public PHX::EvaluatorWithBaseImpl<Traits>,
							    public PHX::EvaluatorDerived<PHX::MyTraits::Jacobian, Traits>  {

  using ScalarT = typename PHX::MyTraits::Jacobian::ScalarT;
  PHX::MDField<const ScalarT,CELL,BASIS> field_at_basis;
  PHX::MDField<ScalarT,CELL,QP,DIM> grad_field_at_qp;

public:
  struct MyDevEval : public PHX::DeviceEvaluator<Traits> {
    PHX::View<const ScalarT**> field_at_basis;
    PHX::View<ScalarT***> grad_field_at_qp;
    KOKKOS_FUNCTION MyDevEval(const PHX::View<const ScalarT**>& in_field_at_basis,
			      const PHX::View<ScalarT***>& in_grad_field_at_qp) :
      field_at_basis(in_field_at_basis), grad_field_at_qp(in_grad_field_at_qp) {}
    KOKKOS_FUNCTION MyDevEval(const MyDevEval& src) = default;
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
