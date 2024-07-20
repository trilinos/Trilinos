// @HEADER
// *****************************************************************************
//        Phalanx: A Partial Differential Equation Field Evaluation 
//       Kernel for Flexible Management of Complex Dependency Chains
//
// Copyright 2008 NTESS and the Phalanx contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

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
void IntegrateDiffusionTerm<EvalT,Traits>::evaluateFields(typename Traits::EvalData workset)
{
  grad_basis = workset.grad_basis_real_;
  weights = workset.weights_;
  cell_measure = workset.det_jac_;
#ifdef PHX_ENABLE_KOKKOS_AMT
  residual_atomic = residual.get_static_view();
#endif
  Kokkos::parallel_for(Kokkos::TeamPolicy<PHX::exec_space>(workset.num_cells_,workset.team_size_,workset.vector_size_),*this);
}

//**********************************************************************
template<typename EvalT, typename Traits>
KOKKOS_INLINE_FUNCTION
void IntegrateDiffusionTerm<EvalT,Traits>::
operator()(const Kokkos::TeamPolicy<PHX::exec_space>::member_type& team) const
{
  const int cell = team.league_rank();

  Kokkos::parallel_for(Kokkos::TeamThreadRange(team,0,grad_basis.extent(2)), [&] (const int& basis) {
      for (int qp = 0; qp < static_cast<int>(grad_basis.extent(1)); ++qp)
        for (int dim = 0; dim < static_cast<int>(grad_basis.extent(3)); ++dim)
#ifdef PHX_ENABLE_KOKKOS_AMT
          residual_atomic(cell,basis) +=  - grad_basis(cell,qp,basis,dim) * flux(cell,qp,dim) * weights(qp) * cell_measure(cell,qp);
#else
          residual(cell,basis) +=  - grad_basis(cell,qp,basis,dim) * flux(cell,qp,dim) * weights(qp) * cell_measure(cell,qp);
#endif
  });
}

//**********************************************************************
