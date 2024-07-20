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
SimpleEvaluator<EvalT,Traits>::
SimpleEvaluator(const Teuchos::RCP<PHX::DataLayout>& a_layout,
                const Teuchos::RCP<PHX::DataLayout>& b_layout,
                const Teuchos::RCP<PHX::DataLayout>& c_layout) :
  a("a",a_layout),
  b("b",b_layout),
  c("c",c_layout)
{
  this->addEvaluatedField(a);
  this->addDependentField(b);
  this->addDependentField(c);
  this->setName("SimpleEvaluator");
}

//**********************************************************************
template<typename EvalT, typename Traits>
void SimpleEvaluator<EvalT,Traits>::evaluateFields(typename Traits::EvalData workset_size)
{
  Kokkos::parallel_for(Kokkos::TeamPolicy<PHX::exec_space>(workset_size,Kokkos::AUTO()),*this);
}

//**********************************************************************
template<typename EvalT, typename Traits>
KOKKOS_INLINE_FUNCTION
void SimpleEvaluator<EvalT,Traits>::
operator()(const Kokkos::TeamPolicy<PHX::exec_space>::member_type& team) const
{
  const int cell = team.league_rank();
  const int num_qp = static_cast<int>(a.extent(1));

  Kokkos::parallel_for(Kokkos::TeamThreadRange(team,0,num_qp), [=] (const int& qp) {
    const int num_dim = static_cast<int>(c.extent(2));
    a(cell,qp) = 0.0;
    for (int i = 0; i < num_dim; ++i)
      a(cell,qp) += b(cell,qp) * b(cell,qp) * c(cell,qp,i);
  });
}

//**********************************************************************
