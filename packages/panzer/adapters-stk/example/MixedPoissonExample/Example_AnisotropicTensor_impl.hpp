// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef __Example_AnisotropicTensor_impl_hpp__
#define __Example_AnisotropicTensor_impl_hpp__

#include "Panzer_HierarchicParallelism.hpp"

namespace Example {

//**********************************************************************
template <typename EvalT,typename Traits>
AnisotropicTensor<EvalT,Traits>::AnisotropicTensor(const std::string & name,
                                                   const panzer::IntegrationRule & ir)
{
  tensor = PHX::MDField<ScalarT,Cell,Point,Dim,Dim>(name, ir.dl_tensor);

  this->addEvaluatedField(tensor);

  std::string n = "Anisotropic Tensor: " + name;
  this->setName(n);
}

//**********************************************************************
template <typename EvalT,typename Traits>
void AnisotropicTensor<EvalT,Traits>::evaluateFields(typename Traits::EvalData workset)
{
  auto t = tensor.get_static_view();
  const AnisotropicTensorEntries K;

  const int num_points = t.extent_int(1);

  auto policy = panzer::HP::inst().teamPolicy<ScalarT,PHX::ExecutionSpace>(workset.num_cells);
  Kokkos::parallel_for("MixedPoisson AnisotropicTensor",policy,KOKKOS_LAMBDA (const Kokkos::TeamPolicy<PHX::ExecutionSpace>::member_type& team) {
    const int cell = team.league_rank();
    Kokkos::parallel_for(Kokkos::TeamThreadRange(team,0,num_points), [&] (const int pt) {
      t(cell,pt,0,0) = K.xx; t(cell,pt,0,1) = K.xy; t(cell,pt,0,2) = K.xz;
      t(cell,pt,1,0) = K.xy; t(cell,pt,1,1) = K.yy; t(cell,pt,1,2) = K.yz;
      t(cell,pt,2,0) = K.xz; t(cell,pt,2,1) = K.yz; t(cell,pt,2,2) = K.zz;
    });
  });
}

//**********************************************************************
}

#endif
