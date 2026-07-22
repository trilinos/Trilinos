// @HEADER
// *****************************************************************************
//        Phalanx: A Partial Differential Equation Field Evaluation 
//       Kernel for Flexible Management of Complex Dependency Chains
//
// Copyright 2008 NTESS and the Phalanx contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Phalanx_FieldTag_Tag.hpp"

//**********************************************************************
template<typename EvalT, typename Traits>
IntegrateSourceTerm<EvalT,Traits>::
IntegrateSourceTerm(const std::string& source_name,
                    const Teuchos::RCP<PHX::DataLayout>& source_layout,
                    const std::string& residual_name,
                    const Teuchos::RCP<PHX::DataLayout>& residual_layout)
{
  PHX::Tag<ScalarT> ft_source(source_name,source_layout);
  PHX::Tag<ScalarT> ft_residual(residual_name,residual_layout);

  this->addContributedField(ft_residual,residual);
  this->addDependentField(ft_source,source);
  this->setName("IntegrateSourceTerm: "+residual_name);
}

//**********************************************************************
template<typename EvalT, typename Traits>
void IntegrateSourceTerm<EvalT,Traits>::evaluateFields(typename Traits::EvalData workset)
{
  basis_view = workset.basis_;
  weights = workset.weights_;
  cell_measure = workset.det_jac_;
  Kokkos::parallel_for(Kokkos::TeamPolicy<PHX::exec_space>(workset.num_cells_,workset.team_size_,workset.vector_size_),*this);
}

//**********************************************************************
template<typename EvalT, typename Traits>
KOKKOS_INLINE_FUNCTION
void IntegrateSourceTerm<EvalT,Traits>::
operator()(const Kokkos::TeamPolicy<PHX::exec_space>::member_type& team) const
{
  const int cell = team.league_rank();

  Kokkos::parallel_for(Kokkos::TeamThreadRange(team,0,basis_view.extent(1)), [&] (const int& basis) {
      for (int qp = 0; qp < static_cast<int>(basis_view.extent(0)); ++qp)
        residual(cell,basis) +=  basis_view(qp,basis) * source(cell,qp) * weights(qp) * cell_measure(cell,qp);
  });
}

//**********************************************************************
