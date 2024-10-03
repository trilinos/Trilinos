// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef EXAMPLE_SIMPLE_SOURCE_IMPL_HPP
#define EXAMPLE_SIMPLE_SOURCE_IMPL_HPP

#include <cmath>

#include "Panzer_BasisIRLayout.hpp"
#include "Panzer_Workset.hpp"
#include "Panzer_Workset_Utilities.hpp"

namespace Example {

//**********************************************************************
template <typename EvalT,typename Traits>
SimpleSource<EvalT,Traits>::SimpleSource(const std::string & name,
                                         const panzer::IntegrationRule & ir)
{
  using Teuchos::RCP;

  Teuchos::RCP<PHX::DataLayout> data_layout = ir.dl_vector;
  ir_degree = ir.cubature_degree;

  source = PHX::MDField<ScalarT,Cell,Point,Dim>(name, data_layout);

  this->addEvaluatedField(source);
  
  std::string n = "Simple Source";
  this->setName(n);
}

//**********************************************************************
template <typename EvalT,typename Traits>
void SimpleSource<EvalT,Traits>::postRegistrationSetup(typename Traits::SetupData sd,           
                                                       PHX::FieldManager<Traits>& /* fm */)
{
  ir_index = panzer::getIntegrationRuleIndex(ir_degree,(*sd.worksets_)[0], this->wda);
}

//**********************************************************************
template <typename EvalT,typename Traits>
void SimpleSource<EvalT,Traits>::evaluateFields(typename Traits::EvalData workset)
{ 
  using panzer::index_t;
  auto ip_coordinates = workset.int_rules[ir_index]->ip_coordinates.get_static_view();
  auto source_v = source.get_static_view();
 
  Kokkos::parallel_for ("SimpleSource", workset.num_cells, KOKKOS_LAMBDA (const index_t cell) {
    for (int point = 0; point < source_v.extent_int(1); ++point) {

      const double & x = ip_coordinates(cell,point,0);
      const double & y = ip_coordinates(cell,point,1);

      source_v(cell,point,0) = 2.0+y-y*y + cos(2.0*M_PI*x)*sin(2.0*M_PI*y);
      source_v(cell,point,1) = 2.0+x-x*x + sin(2.0*M_PI*x)*cos(2.0*M_PI*y);

      // if three d
      if(source_v.extent(2)==3)
        source_v(cell,point,2) = 0.0;
    }
  });
}

//**********************************************************************
}

#endif
