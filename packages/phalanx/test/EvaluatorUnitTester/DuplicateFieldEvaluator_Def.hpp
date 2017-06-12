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
DuplicateFieldEvaluator<EvalT,Traits>::
DuplicateFieldEvaluator(const Teuchos::RCP<const PHX::DataLayout>& a_layout,
                        const Teuchos::RCP<const PHX::DataLayout>& b_layout,
                        const Teuchos::RCP<const PHX::DataLayout>& c_layout) :
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
