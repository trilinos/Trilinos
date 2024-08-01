// @HEADER
// *****************************************************************************
//        Phalanx: A Partial Differential Equation Field Evaluation 
//       Kernel for Flexible Management of Complex Dependency Chains
//
// Copyright 2008 NTESS and the Phalanx contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Phalanx_CreateDeviceEvaluator.hpp"

//**********************************************************************
// Residual Specialization
//**********************************************************************

//**********************************************************************
template<typename EvalT, typename Traits>
ProjectGradientToQP<EvalT,Traits>::
ProjectGradientToQP(const std::string& field_name,
                    const Teuchos::RCP<PHX::DataLayout>& basis_layout,
                    const Teuchos::RCP<PHX::DataLayout>& grad_qp_layout) :
  field_at_basis(field_name,basis_layout),
  grad_field_at_qp(field_name,grad_qp_layout)
{
  this->addEvaluatedField(grad_field_at_qp);
  this->addDependentField(field_at_basis);
  this->setName("ProjectGradientToQP: "+field_name);
}

//**********************************************************************
template<typename EvalT, typename Traits>
PHX::DeviceEvaluator<Traits>*
ProjectGradientToQP<EvalT,Traits>::createDeviceEvaluator() const
{
  using MyDevEval = typename std::conditional<std::is_same<EvalT,PHX::MyTraits::Residual>::value,MyDevEvalResidual,MyDevEvalJacobian>::type;
  return PHX::createDeviceEvaluator<MyDevEval,Traits,PHX::exec_space,PHX::mem_space>(field_at_basis.get_static_view(),
                                                                                     grad_field_at_qp.get_static_view());
}

//**********************************************************************
template<typename EvalT, typename Traits>
void ProjectGradientToQP<EvalT,Traits>::evaluateFields(typename Traits::EvalData workset)
{
  using MyDevEval = typename std::conditional<std::is_same<EvalT,PHX::MyTraits::Residual>::value,MyDevEvalResidual,MyDevEvalJacobian>::type;
  auto e = PHX::make_dev_eval(MyDevEval(field_at_basis.get_static_view(), grad_field_at_qp.get_static_view()),workset);
  Kokkos::parallel_for(Kokkos::TeamPolicy<PHX::exec_space>(workset.num_cells_,workset.team_size_,workset.vector_size_),e);
}

//**********************************************************************
template<typename EvalT, typename Traits>
KOKKOS_INLINE_FUNCTION
void ProjectGradientToQP<EvalT,Traits>::MyDevEvalResidual::
evaluate(const typename PHX::DeviceEvaluator<Traits>::member_type& team,
         typename Traits::EvalData workset)
{
#if (__GNUC__ == 5) || (__GNUC__ == 6)
  int cell = team.league_rank(); // remove const for gcc 5/6 bug in nested lambdas
#else
  const int cell = team.league_rank();
#endif
  auto grad_basis_view = workset.grad_basis_real_;
  Kokkos::parallel_for(Kokkos::TeamThreadRange(team,0,grad_field_at_qp.extent(1)), [&] (const int& qp) {
    Kokkos::single(Kokkos::PerThread(team), [&] (){ 
      grad_field_at_qp(cell,qp,0) = ScalarT(0.0);
      grad_field_at_qp(cell,qp,1) = ScalarT(0.0);
      grad_field_at_qp(cell,qp,2) = ScalarT(0.0);
    });
    team.team_barrier();
    Kokkos::parallel_for(Kokkos::ThreadVectorRange(team,field_at_basis.extent(1)), [&] (const int& basis) {
      for (int dim = 0; dim < static_cast<int>(grad_field_at_qp.extent(2)); ++dim)
	grad_field_at_qp(cell,qp,dim) += field_at_basis(cell,basis) * grad_basis_view(cell,qp,basis,dim);
    });
  });
  team.team_barrier();
}

//**********************************************************************
template<typename EvalT, typename Traits>
KOKKOS_INLINE_FUNCTION
void ProjectGradientToQP<EvalT,Traits>::MyDevEvalJacobian::
evaluate(const typename PHX::DeviceEvaluator<Traits>::member_type& team,
         typename Traits::EvalData workset)
{
#if (__GNUC__ == 5) || (__GNUC__ == 6)
  int cell = team.league_rank(); // remove const for gcc 5/6 bug in nested lambdas
#else
  const int cell = team.league_rank();
#endif
  auto grad_basis_view = workset.grad_basis_real_;
  Kokkos::parallel_for(Kokkos::TeamThreadRange(team,0,grad_field_at_qp.extent(1)), KOKKOS_LAMBDA (const int& qp) {
    grad_field_at_qp(cell,qp,0) = ScalarT(0.0);
    grad_field_at_qp(cell,qp,1) = ScalarT(0.0);
    grad_field_at_qp(cell,qp,2) = ScalarT(0.0);
    for (int basis = 0; basis < static_cast<int>(field_at_basis.extent(1)); ++basis)
      for (int dim = 0; dim < static_cast<int>(grad_field_at_qp.extent(2)); ++dim)
	grad_field_at_qp(cell,qp,dim) += field_at_basis(cell,basis) * grad_basis_view(cell,qp,basis,dim);
  });
  team.team_barrier(); // needed for device dag
}

//**********************************************************************
