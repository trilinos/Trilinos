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

#include "Phalanx_CreateDeviceEvaluator.hpp"

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
PHX::DeviceEvaluator<Traits>*
ProjectValueToQP<EvalT,Traits>::createDeviceEvaluator() const
{
  using MyDevEval = typename std::conditional<std::is_same<EvalT,PHX::MyTraits::Residual>::value,MyDevEvalResidual,MyDevEvalJacobian>::type;
  return PHX::createDeviceEvaluator<MyDevEval,Traits,PHX::exec_space,PHX::mem_space>(field_at_basis.get_static_view(),
                                                                                     field_at_qp.get_static_view());
}

//**********************************************************************
template<typename EvalT, typename Traits>
void ProjectValueToQP<EvalT,Traits>::evaluateFields(typename Traits::EvalData workset)
{
  using MyDevEval = typename std::conditional<std::is_same<EvalT,PHX::MyTraits::Residual>::value,MyDevEvalResidual,MyDevEvalJacobian>::type;
  auto e = PHX::make_dev_eval(MyDevEval(field_at_basis.get_static_view(),field_at_qp.get_static_view()),workset);
  Kokkos::parallel_for(Kokkos::TeamPolicy<PHX::exec_space>(workset.num_cells_,workset.team_size_,workset.vector_size_),e);
}

//**********************************************************************
template<typename EvalT, typename Traits>
KOKKOS_INLINE_FUNCTION
void ProjectValueToQP<EvalT,Traits>::MyDevEvalResidual::
evaluate(const typename PHX::DeviceEvaluator<Traits>::member_type& team,
         typename Traits::EvalData workset)
{
#if (__GNUC__ == 5) || (__GNUC__ == 6)
  int cell = team.league_rank(); // remove const for gcc 5/6 bug in nested lambdas
#else
  const int cell = team.league_rank();
#endif
  auto basis_view = workset.basis_;
  Kokkos::parallel_for(Kokkos::TeamThreadRange(team,0,field_at_qp.extent(1)), [&] (const int& qp) {
    field_at_qp(cell,qp) = ScalarT(0.0);
    Kokkos::parallel_for(Kokkos::ThreadVectorRange(team,field_at_basis.extent(1)), [&] (const int& basis) {
	field_at_qp(cell,qp) += field_at_basis(cell,basis) * basis_view(qp,basis);
    });
  });
}

//**********************************************************************
template<typename EvalT, typename Traits>
KOKKOS_INLINE_FUNCTION
void ProjectValueToQP<EvalT,Traits>::MyDevEvalJacobian::
evaluate(const typename PHX::DeviceEvaluator<Traits>::member_type& team,
         typename Traits::EvalData workset)
{
  const int cell = team.league_rank();
  auto basis_view = workset.basis_;
  Kokkos::parallel_for(Kokkos::TeamThreadRange(team,0,field_at_qp.extent(1)), [&] (const int& qp) {
    field_at_qp(cell,qp) = ScalarT(0.0);
    for (int basis = 0; basis < static_cast<int>(field_at_basis.extent(1)); ++basis)
      field_at_qp(cell,qp) += field_at_basis(cell,basis) * basis_view(qp,basis);
  });
}

//**********************************************************************
