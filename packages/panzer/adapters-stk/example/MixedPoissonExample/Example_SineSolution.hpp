// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef __Example_SineSolution_hpp__
#define __Example_SineSolution_hpp__

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

/** The analytic solution to the mixed poisson equation for the sine source.
  */
template<typename EvalT, typename Traits>
class SineSolution : public panzer::EvaluatorWithBaseImpl<Traits>,
                        public PHX::EvaluatorDerived<EvalT, Traits>  {

public:
    SineSolution(const std::string & name,
                       const panzer::IntegrationRule & ir);
                                                                        
    void postRegistrationSetup(typename Traits::SetupData d,           
                               PHX::FieldManager<Traits>& fm);        
                                                                     
    void evaluateFields(typename Traits::EvalData d);               


private:
  typedef typename EvalT::ScalarT ScalarT;

  // Simulation solution
  PHX::MDField<ScalarT,Cell,Point> solution;
  int ir_degree, ir_index;
};

}

#include "Example_SineSolution_impl.hpp"

#endif
