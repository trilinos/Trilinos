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
ZeroContributedField<EvalT,Traits>::
ZeroContributedField(const std::string& field_name,
                     const Teuchos::RCP<PHX::DataLayout>& layout) :
  field(field_name,layout)
{
  // "Evalauted" is always called before "Contributed" for the same field
  this->addEvaluatedField(field);
  this->setName("ZeroContributedField: " + field.fieldTag().identifier());
}

//**********************************************************************
template<typename EvalT, typename Traits>
PHX::DeviceEvaluator<Traits>*
ZeroContributedField<EvalT,Traits>::createDeviceEvaluator() const
{
  using MyDevEval = typename std::conditional<std::is_same<EvalT,PHX::MyTraits::Residual>::value,MyDevEvalResidual,MyDevEvalJacobian>::type;
  return PHX::createDeviceEvaluator<MyDevEval,Traits,PHX::exec_space,PHX::mem_space>(field.get_view());
}

//**********************************************************************
template<typename EvalT, typename Traits>
KOKKOS_INLINE_FUNCTION
void
ZeroContributedField<EvalT,Traits>::MyDevEvalResidual::
evaluate(const typename PHX::DeviceEvaluator<Traits>::member_type& team,
         typename Traits::EvalData )
{
  const int cell = team.league_rank();
  const int num_basis = field_.extent(1);
  if (team.team_rank() == 0) {
    Kokkos::parallel_for(Kokkos::ThreadVectorRange(team,num_basis), [&] (const int& basis) {
      field_(cell,basis) = ScalarT(0.);
    });
  }
}

//**********************************************************************
template<typename EvalT, typename Traits>
KOKKOS_INLINE_FUNCTION
void
ZeroContributedField<EvalT,Traits>::MyDevEvalJacobian::
evaluate(const typename PHX::DeviceEvaluator<Traits>::member_type& team,
         typename Traits::EvalData )
{
  const int cell = team.league_rank();
  const int num_basis = field_.extent(1);
  Kokkos::parallel_for(Kokkos::TeamThreadRange(team,0,num_basis), [&] (const int& basis) {
    field_(cell,basis) = ScalarT(0.);
  });
}

//**********************************************************************
template<typename EvalT, typename Traits>
void
ZeroContributedField<EvalT,Traits>::evaluateFields(typename Traits::EvalData )
{
  field.deep_copy(ScalarT(0.0));
}

//**********************************************************************
