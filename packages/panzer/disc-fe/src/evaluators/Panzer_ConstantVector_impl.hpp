// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
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
