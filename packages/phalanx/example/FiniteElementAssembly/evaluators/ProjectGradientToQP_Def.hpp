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
void ProjectGradientToQP<EvalT,Traits>::evaluateFields(typename Traits::EvalData workset)
{
  grad_basis_view = workset.grad_basis_real_;
  Kokkos::parallel_for(Kokkos::TeamPolicy<PHX::exec_space>(workset.num_cells_,workset.team_size_,workset.vector_size_),*this);
}

//**********************************************************************
template<typename EvalT, typename Traits>
KOKKOS_INLINE_FUNCTION
void ProjectGradientToQP<EvalT,Traits>::
operator()(const Kokkos::TeamPolicy<PHX::exec_space>::member_type& team) const
{
  const int cell = team.league_rank();
  Kokkos::parallel_for(Kokkos::TeamThreadRange(team,0,grad_field_at_qp.extent(1)), [&] (const int& qp) {
      grad_field_at_qp(cell,qp,0) = ScalarT(0.0);
      grad_field_at_qp(cell,qp,1) = ScalarT(0.0);
      grad_field_at_qp(cell,qp,2) = ScalarT(0.0);
      for (int basis = 0; basis < static_cast<int>(field_at_basis.extent(1)); ++basis)
        for (int dim = 0; dim < static_cast<int>(grad_field_at_qp.extent(2)); ++dim)
          grad_field_at_qp(cell,qp,dim) += field_at_basis(cell,basis) * grad_basis_view(cell,qp,basis,dim);
  });
}

//**********************************************************************
