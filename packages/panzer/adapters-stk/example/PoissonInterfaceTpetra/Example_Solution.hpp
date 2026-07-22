// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef __Example_Solution_hpp__
#define __Example_Solution_hpp__

#include "PanzerAdaptersSTK_config.hpp"

#include "Phalanx_Evaluator_WithBaseImpl.hpp"
#include "Phalanx_Evaluator_Derived.hpp"
#include "Phalanx_FieldManager.hpp"

#include "Panzer_Dimension.hpp"
#include "Panzer_FieldLibrary.hpp"

#include <string>

#include "Panzer_Evaluator_WithBaseImpl.hpp"

namespace Example {
    
using panzer::Cell;
using panzer::Point;
using panzer::Dim;

/** The analytic solution to the Poisson equation for the sin*cos source with
  * DBC on left and right.
  */
template<typename EvalT, typename Traits>
class Solution : public panzer::EvaluatorWithBaseImpl<Traits>,
                 public PHX::EvaluatorDerived<EvalT, Traits>  {

public:
  Solution(const std::string & name,
           const panzer::IntegrationRule & ir,
           const bool linear_Robin = true);
                                                                        
  void postRegistrationSetup(typename Traits::SetupData d,           
                             PHX::FieldManager<Traits>& fm);        
                                                                     
  void evaluateFields(typename Traits::EvalData d);               


private:
  typedef typename EvalT::ScalarT ScalarT;

  PHX::MDField<ScalarT,Cell,Point> solution;
  int ir_degree, ir_index;
  bool linear_Robin;
};

}

#include "Example_Solution_impl.hpp"

#endif
