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
ProjectValueToQP<EvalT,Traits>::
ProjectValueToQP(const std::string& field_name,
                 const Teuchos::RCP<PHX::DataLayout>& basis_layout,
                 const Teuchos::RCP<PHX::DataLayout>& qp_layout) :
  field_at_basis(field_name,basis_layout),
  field_at_qp(field_name,qp_layout)
{
  this->addEvaluatedField(field_at_qp);
  this->addDependentField(field_at_basis);
  this->setName("ProjectValueToQP: "+field_name);
}

//**********************************************************************
template<typename EvalT, typename Traits>
void ProjectValueToQP<EvalT,Traits>::evaluateFields(typename Traits::EvalData workset)
{
  basis_view = workset.basis_;
  Kokkos::parallel_for(Kokkos::TeamPolicy<PHX::exec_space>(workset.num_cells_,workset.team_size_,workset.vector_size_),*this);
}

//**********************************************************************
template<typename EvalT, typename Traits>
KOKKOS_INLINE_FUNCTION
void ProjectValueToQP<EvalT,Traits>::
operator()(const Kokkos::TeamPolicy<PHX::exec_space>::member_type& team) const
{
  const int cell = team.league_rank();
  Kokkos::parallel_for(Kokkos::TeamThreadRange(team,0,field_at_qp.extent(1)), [&] (const int& qp) {
      field_at_qp(cell,qp) = ScalarT(0.0);
      for (int basis = 0; basis < static_cast<int>(field_at_basis.extent(1)); ++basis)
        field_at_qp(cell,qp) += field_at_basis(cell,basis) * basis_view(qp,basis);
  });
}

//**********************************************************************
