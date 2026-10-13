// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef __Example_LinearCoefficient_impl_hpp__
#define __Example_LinearCoefficient_impl_hpp__

#include "Panzer_Workset_Utilities.hpp"
#include "Panzer_HierarchicParallelism.hpp"

namespace Example {

//**********************************************************************
template <typename EvalT,typename Traits>
LinearCoefficient<EvalT,Traits>::LinearCoefficient(const std::string & name,
                                                   const panzer::IntegrationRule & ir,
                                                   const int in_direction)
  : direction(in_direction)
{
  ir_degree = ir.cubature_degree;

  coeff = PHX::MDField<ScalarT,Cell,Point>(name, ir.dl_scalar);

  this->addEvaluatedField(coeff);

  std::string n = "Linear Coefficient: " + name;
  this->setName(n);
}

//**********************************************************************
template <typename EvalT,typename Traits>
void LinearCoefficient<EvalT,Traits>::postRegistrationSetup(typename Traits::SetupData sd,
                                                            PHX::FieldManager<Traits>& /* fm */)
{
  ir_index = panzer::getIntegrationRuleIndex(ir_degree,(*sd.worksets_)[0], this->wda);
}

//**********************************************************************
template <typename EvalT,typename Traits>
void LinearCoefficient<EvalT,Traits>::evaluateFields(typename Traits::EvalData workset)
{
  auto c = coeff.get_static_view();
  auto ip_coords = this->wda(workset).int_rules[ir_index]->ip_coordinates.get_static_view();
  const int d = direction;

  const int num_points = c.extent_int(1);

  auto policy = panzer::HP::inst().teamPolicy<ScalarT,PHX::ExecutionSpace>(workset.num_cells);
  Kokkos::parallel_for("MixedPoisson LinearCoefficient",policy,KOKKOS_LAMBDA (const Kokkos::TeamPolicy<PHX::ExecutionSpace>::member_type& team) {
    const int cell = team.league_rank();
    Kokkos::parallel_for(Kokkos::TeamThreadRange(team,0,num_points), [&] (const int pt) {
      c(cell,pt) = 1.0 + ip_coords(cell,pt,d);
    });
  });
}

//**********************************************************************
}

#endif
