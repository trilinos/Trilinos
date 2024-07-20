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
DuplicateFieldEvaluator<EvalT,Traits>::
DuplicateFieldEvaluator(const Teuchos::RCP<PHX::DataLayout>& a_layout,
                        const Teuchos::RCP<PHX::DataLayout>& b_layout,
                        const Teuchos::RCP<PHX::DataLayout>& c_layout) :
  a("a",a_layout),
  b1("b",b_layout),
  b2("b",b_layout), // purposely duplicate b1 for corner case
  c("c",c_layout)
{
  this->addEvaluatedField(a);
  this->addDependentField(b1);
  this->addDependentField(b2);
  this->addDependentField(c);
  this->setName("DuplicateFieldEvaluator");
}

//**********************************************************************
template<typename EvalT, typename Traits>
void DuplicateFieldEvaluator<EvalT,Traits>::evaluateFields(typename Traits::EvalData workset_size)
{
  Kokkos::parallel_for(Kokkos::TeamPolicy<PHX::exec_space>(workset_size,Kokkos::AUTO()),*this);
}

//**********************************************************************
template<typename EvalT, typename Traits>
KOKKOS_INLINE_FUNCTION
void DuplicateFieldEvaluator<EvalT,Traits>::
operator()(const Kokkos::TeamPolicy<PHX::exec_space>::member_type& team) const
{
  const int cell = team.league_rank();
  const int num_qp = static_cast<int>(a.extent(1));

  Kokkos::parallel_for(Kokkos::TeamThreadRange(team,0,num_qp), [=] (const int& qp) {
    const int num_dim = static_cast<int>(c.extent(2));
    a(cell,qp) = 0.0;
    for (int i = 0; i < num_dim; ++i)
      a(cell,qp) += b1(cell,qp) * b2(cell,qp) * c(cell,qp,i);
  });
}

//**********************************************************************
