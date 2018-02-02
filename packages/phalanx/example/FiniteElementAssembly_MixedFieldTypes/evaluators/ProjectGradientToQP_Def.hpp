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
  Kokkos::parallel_for(Kokkos::TeamThreadRange(team,0,grad_field_at_qp.extent(1)), KOKKOS_LAMBDA (const int& qp) {
      grad_field_at_qp(cell,qp,0) = ScalarT(0.0);
      grad_field_at_qp(cell,qp,1) = ScalarT(0.0);
      grad_field_at_qp(cell,qp,2) = ScalarT(0.0);
      for (int basis = 0; basis < static_cast<int>(field_at_basis.extent(1)); ++basis)
        for (int dim = 0; dim < static_cast<int>(grad_field_at_qp.extent(2)); ++dim)
          grad_field_at_qp(cell,qp,dim) += field_at_basis(cell,basis) * grad_basis_view(cell,qp,basis,dim);
  });
}

//**********************************************************************
