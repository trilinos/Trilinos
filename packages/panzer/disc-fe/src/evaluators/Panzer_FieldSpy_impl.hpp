// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef __Panzer_FieldSpy_impl_hpp__
#define __Panzer_FieldSpy_impl_hpp__

#include <cmath>

#include "Phalanx_DataLayout_MDALayout.hpp"

#include "Panzer_Workset.hpp"

namespace panzer {

//**********************************************************************
template <typename EvalT,typename Traits>
FieldSpy<EvalT,Traits>::FieldSpy(const std::string & name,
                                 const Teuchos::RCP<PHX::DataLayout> & data_layout)
{
  using Teuchos::RCP;
  using Teuchos::rcp;

  dummyField = rcp(new PHX::Tag<ScalarT>("Field Spy: " + name,rcp(new PHX::MDALayout<panzer::Dummy>(0))));
  this->addEvaluatedField(*dummyField);

  source = PHX::MDField<const ScalarT,Cell,Point>(name, data_layout);
  this->addDependentField(source);
  
  std::string n = "Field Spy";
  this->setName(n);
}

//**********************************************************************
template <typename EvalT,typename Traits>
void FieldSpy<EvalT,Traits>::evaluateFields(typename Traits::EvalData workset)
{ 
  std::cout << "SPY: Name = \"" << source.fieldTag().identifier() << "\" at t = " << workset.time << "\n";
  for (index_t cell=0;cell<workset.num_cells;++cell) {
    std::cout << "SPY: ";
    for (int point = 0; point < source.extent_int(1); ++point) {
      std::cout << Sacado::scalarValue(source(cell,point)) << " ";
    }
    std::cout << std::endl;
  }
}

//**********************************************************************
}

#endif
