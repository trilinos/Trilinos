// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef __Example_SineSolution_impl_hpp__
#define __Example_SineSolution_impl_hpp__

#include <cmath>

#include "Panzer_BasisIRLayout.hpp"
#include "Panzer_Workset.hpp"
#include "Panzer_Workset_Utilities.hpp"

namespace Example {

//**********************************************************************
template <typename EvalT,typename Traits>
SineSolution<EvalT,Traits>::SineSolution(const std::string & name,
                                         const panzer::IntegrationRule & ir)
{
  using Teuchos::RCP;

  Teuchos::RCP<PHX::DataLayout> data_layout = ir.dl_scalar;
  ir_degree = ir.cubature_degree;

  solution = PHX::MDField<ScalarT,Cell,Point>(name, data_layout);

  this->addEvaluatedField(solution);
  
  std::string n = "Sine Solution";
  this->setName(n);
}

//**********************************************************************
template <typename EvalT,typename Traits>
void SineSolution<EvalT,Traits>::postRegistrationSetup(typename Traits::SetupData sd,           
                                                       PHX::FieldManager<Traits>& /* fm */)
{
  ir_index = panzer::getIntegrationRuleIndex(ir_degree,(*sd.worksets_)[0], this->wda);
}

//**********************************************************************
template <typename EvalT,typename Traits>
void SineSolution<EvalT,Traits>::evaluateFields(typename Traits::EvalData workset)
{ 
  using panzer::index_t;
  auto ip_coordinates = this->wda(workset).int_rules[ir_index]->ip_coordinates.get_static_view();
  auto solution_v = solution.get_static_view();

  Kokkos::parallel_for (workset.num_cells, KOKKOS_LAMBDA (const index_t cell) {
    for (int point = 0; point < solution_v.extent_int(1); ++point) {

      const double & x = ip_coordinates(cell,point,0);
      const double & y = ip_coordinates(cell,point,1);
      const double & z = ip_coordinates(cell,point,2);

      solution_v(cell,point) = std::sin(2.0*M_PI*x)*std::sin(2*M_PI*y)*std::sin(2.0*M_PI*z);
    }
  });
}

//**********************************************************************
}

#endif
