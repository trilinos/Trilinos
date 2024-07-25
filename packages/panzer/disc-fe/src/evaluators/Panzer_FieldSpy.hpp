// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef __Panzer_FieldSpy_hpp__
#define __Panzer_FieldSpy_hpp__

#include "PanzerDiscFE_config.hpp"
#include "Phalanx_Evaluator_WithBaseImpl.hpp"
#include "Phalanx_Evaluator_Derived.hpp"
#include "Phalanx_FieldManager.hpp"
#include "Panzer_Dimension.hpp"
#include "Panzer_FieldLibrary.hpp"
#include <string>

namespace panzer {
    
/** A Output evaluator for writing out fields.
  */
template<typename EvalT, typename Traits>
class FieldSpy : public PHX::EvaluatorWithBaseImpl<Traits>,
                 public PHX::EvaluatorDerived<EvalT, Traits>  {

public:
    FieldSpy(const std::string & name,
             const Teuchos::RCP<PHX::DataLayout> & data_layout);
                                                                        
    void evaluateFields(typename Traits::EvalData d);               

    const PHX::FieldTag & getRequiredFieldTag() const 
    { return *dummyField; }

private:
  typedef typename EvalT::ScalarT ScalarT;

  Teuchos::RCP<PHX::FieldTag> dummyField;
  PHX::MDField<const ScalarT,panzer::Cell,panzer::Point> source;
};

}

#endif
