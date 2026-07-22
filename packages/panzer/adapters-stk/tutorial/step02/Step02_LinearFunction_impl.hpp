// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef __Step02_LinearFunction_impl_hpp__
#define __Step02_LinearFunction_impl_hpp__

#include <cmath>

#include "Panzer_BasisIRLayout.hpp"
#include "Panzer_Workset.hpp"
#include "Panzer_Workset_Utilities.hpp"

namespace user_app {

//**********************************************************************
template <typename EvalT,typename Traits>
LinearFunction<EvalT,Traits>::LinearFunction(const std::string & name,
                                             double acoeff,double bcoeff,
                                             const panzer::IntegrationRule & ir)
  : acoeff_(acoeff) 
  , bcoeff_(bcoeff) 
  , ir_degree_(ir.cubature_degree)
{
  using Teuchos::RCP;

  Teuchos::RCP<PHX::DataLayout> data_layout = ir.dl_scalar;

  result = PHX::MDField<ScalarT,panzer::Cell,panzer::Point>(name, data_layout);
  this->addEvaluatedField(result);

  this->setName("Linear Function("+name+")");
}

//**********************************************************************
template <typename EvalT,typename Traits>
void LinearFunction<EvalT,Traits>::postRegistrationSetup(typename Traits::SetupData sd,           
                                                          PHX::FieldManager<Traits>& /* fm */)
{
  ir_index_ = panzer::getIntegrationRuleIndex(ir_degree_,(*sd.worksets_)[0]);
}

//**********************************************************************
template <typename EvalT,typename Traits>
void LinearFunction<EvalT,Traits>::evaluateFields(typename Traits::EvalData workset)
{ 
  for (panzer::index_t cell = 0; cell < workset.num_cells; ++cell) {
    for (int point = 0; point < result.extent_int(1); ++point) {

      const double& x = workset.int_rules[ir_index_]->ip_coordinates(cell,point,0);
      const double& y = workset.int_rules[ir_index_]->ip_coordinates(cell,point,1);

      result(cell,point) = acoeff_*x + acoeff_*y + bcoeff_;
    }
  }
}

//**********************************************************************
}

#endif
