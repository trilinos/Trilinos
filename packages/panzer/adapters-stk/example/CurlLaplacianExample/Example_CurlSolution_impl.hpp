// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef __Example_CurlSolution_impl_hpp__
#define __Example_CurlSolution_impl_hpp__

#include <cmath>

#include "Panzer_BasisIRLayout.hpp"
#include "Panzer_Workset.hpp"
#include "Panzer_Workset_Utilities.hpp"

namespace Example {

//**********************************************************************
template <typename EvalT,typename Traits>
CurlSolution<EvalT,Traits>::CurlSolution(const std::string & name,
                                         const panzer::IntegrationRule & ir)
{
  using Teuchos::RCP;

  Teuchos::RCP<PHX::DataLayout> dl_vector = ir.dl_vector;
  Teuchos::RCP<PHX::DataLayout> dl_scalar = ir.dl_scalar;
  ir_degree = ir.cubature_degree;

  solution      = PHX::MDField<ScalarT,Cell,Point,Dim>(name, dl_vector);
  solution_curl = PHX::MDField<ScalarT,Cell,Point>("CURL_"+name, dl_scalar);

  // only thing that works
  TEUCHOS_ASSERT(ir.spatial_dimension==2);

  this->addEvaluatedField(solution);
  this->addEvaluatedField(solution_curl);
  
  std::string n = "Curl Solution";
  this->setName(n);
}

//**********************************************************************
template <typename EvalT,typename Traits>
void CurlSolution<EvalT,Traits>::postRegistrationSetup(typename Traits::SetupData sd,           
                                                       PHX::FieldManager<Traits>& /* fm */)
{
  ir_index = panzer::getIntegrationRuleIndex(ir_degree,(*sd.worksets_)[0], this->wda);
}

//**********************************************************************
template <typename EvalT,typename Traits>
void CurlSolution<EvalT,Traits>::evaluateFields(typename Traits::EvalData workset)
{ 
  using panzer::index_t;
  auto ip_coordinates = this->wda(workset).int_rules[ir_index]->ip_coordinates.get_static_view();
  auto solution_v = solution.get_static_view();
  auto solution_curl_v = solution_curl.get_static_view();

  Kokkos::parallel_for (workset.num_cells, KOKKOS_LAMBDA (const index_t cell) {
    for (int point = 0; point < solution_v.extent_int(1); ++point) {

      const double & x = ip_coordinates(cell,point,0);
      const double & y = ip_coordinates(cell,point,1);

      solution_v(cell,point,0) = -(y-1.0)*y + cos(2.0*M_PI*x)*sin(2.0*M_PI*y);
      solution_v(cell,point,1) = -(x-1.0)*x + sin(2.0*M_PI*x)*cos(2.0*M_PI*y);

      solution_curl_v(cell,point) = -2.0*x+2.0*y;
    }
  });
}

//**********************************************************************
}

#endif
