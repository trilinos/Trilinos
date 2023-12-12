// @HEADER
// ***********************************************************************
//
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//                 Copyright (2011) Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
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
// Questions? Contact Roger P. Pawlowski (rppawlo@sandia.gov) and
// Eric C. Cyr (eccyr@sandia.gov)
// ***********************************************************************
// @HEADER

#ifndef PANZER_CONSTANT_VECTOR_IMPL_HPP
#define PANZER_CONSTANT_VECTOR_IMPL_HPP

namespace panzer {

//**********************************************************************
template<typename EvalT, typename Traits>
ConstantVector<EvalT, Traits>::
ConstantVector(
  const Teuchos::ParameterList& p) :
  vec_(p.get<std::string>("Name"), 
       p.get< Teuchos::RCP<PHX::DataLayout> >("Data Layout") )
{
  this->addEvaluatedField(vec_);

  // Make this unshared so that it is not overwritten
  this->addUnsharedField(vec_.fieldTag().clone());

  const int dim = vec_.fieldTag().dataLayout().extent(2);

  vals_ = Kokkos::View<double*>("ConstantVector::vals",dim);
  auto vals_host = Kokkos::create_mirror_view(Kokkos::HostSpace(),vals_);

  vals_host(0) = p.get<double>("Value X");
  if(dim>1)
    vals_host(1) = p.get<double>("Value Y");
  if(dim>2)
    vals_host(2) = p.get<double>("Value Z");

  Kokkos::deep_copy(vals_,vals_host);

  std::string n = "ConstantVector: " + vec_.fieldTag().name();
  this->setName(n);
}

//**********************************************************************
template<typename EvalT, typename Traits>
void
ConstantVector<EvalT, Traits>::
postRegistrationSetup(typename Traits::SetupData  /* worksets */,
                      PHX::FieldManager<Traits>&  fm)
{
  auto vals = this->vals_;
  auto vec = this->vec_;
  Kokkos::MDRangePolicy<PHX::Device,Kokkos::Rank<3>> policy({0,0,0},{static_cast<int64_t>(vec.extent(0)),
        static_cast<int64_t>(vec.extent(1)),static_cast<int64_t>(vec.extent(2))});
  Kokkos::parallel_for("panzer::ConstantVector",policy,KOKKOS_LAMBDA(const int c, const int p, const int d){
    vec(c,p,d) = vals(d);
  });
}

//**********************************************************************
template<typename EvalT, typename Traits>
void
ConstantVector<EvalT, Traits>::
evaluateFields(typename Traits::EvalData  /* d */)
{}

//**********************************************************************

}

#endif
