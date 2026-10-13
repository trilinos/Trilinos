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

  Kokkos::MDRangePolicy<PHX::exec_space,Kokkos::Rank<2>> policy({0,0},{static_cast<int>(workset.num_cells),c.extent_int(1)});
  Kokkos::parallel_for("MixedPoisson LinearCoefficient",policy,KOKKOS_LAMBDA (const int cell,const int pt) {
    c(cell,pt) = 1.0 + ip_coords(cell,pt,d);
  });
}

//**********************************************************************
}

#endif
