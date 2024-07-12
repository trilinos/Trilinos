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
template<typename EvalT, typename Traits>
IntegrateDiffusionTerm<EvalT,Traits>::
IntegrateDiffusionTerm(const std::string& flux_name,
                       const Teuchos::RCP<PHX::DataLayout>& flux_layout,
                       const std::string& residual_name,
                       const Teuchos::RCP<PHX::DataLayout>& residual_layout) :
  flux(flux_name,flux_layout),
  residual(residual_name,residual_layout)
{
  this->addContributedField(residual);
  this->addDependentField(flux);
  this->setName("IntegrateDiffusionTerm: "+residual_name);
}

//**********************************************************************
template<typename EvalT, typename Traits>
PHX::DeviceEvaluator<Traits>*
IntegrateDiffusionTerm<EvalT,Traits>::createDeviceEvaluator() const
{
  using MyDevEval = typename std::conditional<std::is_same<EvalT,PHX::MyTraits::Residual>::value,MyDevEvalResidual,MyDevEvalJacobian>::type;
  return PHX::createDeviceEvaluator<MyDevEval,Traits,PHX::exec_space,PHX::mem_space>(flux.get_static_view(),
                                                                                     residual.get_static_view());
}

//**********************************************************************
template<typename EvalT, typename Traits>
void IntegrateDiffusionTerm<EvalT,Traits>::evaluateFields(typename Traits::EvalData workset)
{
  using MyDevEval = typename std::conditional<std::is_same<EvalT,PHX::MyTraits::Residual>::value,MyDevEvalResidual,MyDevEvalJacobian>::type;
  auto e = PHX::make_dev_eval(MyDevEval(flux.get_static_view(),residual.get_static_view()),workset);
  Kokkos::parallel_for(Kokkos::TeamPolicy<PHX::exec_space>(workset.num_cells_,workset.team_size_,workset.vector_size_),e);
}

//**********************************************************************
template<typename EvalT, typename Traits>
KOKKOS_FUNCTION
void IntegrateDiffusionTerm<EvalT,Traits>::MyDevEvalResidual::
evaluate(const typename PHX::DeviceEvaluator<Traits>::member_type& team,
         typename Traits::EvalData workset)
{
  const int cell = team.league_rank();
  const auto& grad_basis = workset.grad_basis_real_;
  const auto& weights = workset.weights_;
  const auto& cell_measure = workset.det_jac_;
  const int num_basis = grad_basis.extent(2);
  const int num_qp = grad_basis.extent(1);

  Kokkos::parallel_for(Kokkos::TeamThreadRange(team,0,num_basis), [&] (const int& basis) {
    Kokkos::parallel_for(Kokkos::ThreadVectorRange(team,num_qp), [&] (const int& qp) {
      for (int dim = 0; dim < static_cast<int>(grad_basis.extent(3)); ++dim)
	residual(cell,basis) +=  - grad_basis(cell,qp,basis,dim) * flux(cell,qp,dim) * weights(qp) * cell_measure(cell,qp);
    });
  });
  team.team_barrier();
}

//**********************************************************************
template<typename EvalT, typename Traits>
KOKKOS_FUNCTION
void IntegrateDiffusionTerm<EvalT,Traits>::MyDevEvalJacobian::
evaluate(const typename PHX::DeviceEvaluator<Traits>::member_type& team,
         typename Traits::EvalData workset)
{
  const int cell = team.league_rank();
  const auto& grad_basis = workset.grad_basis_real_;
  const auto& weights = workset.weights_;
  const auto& cell_measure = workset.det_jac_;

  Kokkos::parallel_for(Kokkos::TeamThreadRange(team,0,grad_basis.extent(2)), [&] (const int& basis) {
    for (int qp = 0; qp < static_cast<int>(grad_basis.extent(1)); ++qp)
      for (int dim = 0; dim < static_cast<int>(grad_basis.extent(3)); ++dim)
	residual(cell,basis) +=  - grad_basis(cell,qp,basis,dim) * flux(cell,qp,dim) * weights(qp) * cell_measure(cell,qp);
  });
}

//**********************************************************************
