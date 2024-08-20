// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef __Step02_LinearFunction_hpp__
#define __Step02_LinearFunction_hpp__

#include "Phalanx_Evaluator_WithBaseImpl.hpp"
#include "Phalanx_Evaluator_Derived.hpp"
#include "Phalanx_FieldManager.hpp"

#include "Panzer_Dimension.hpp"
#include "Panzer_FieldLibrary.hpp"

#include <string>

namespace user_app {
    
/** A source for the poisson equation that results in the solution
  * (with homogeneous dirichlet conditions) \f$sin(2\pi x)sin(2\pi y)\f$.
  */
template<typename EvalT, typename Traits>
class LinearFunction : public PHX::EvaluatorWithBaseImpl<Traits>,
                       public PHX::EvaluatorDerived<EvalT, Traits>  {

public:
    LinearFunction(const std::string & name,
                   double acoeff,double bcoeff,
                   const panzer::IntegrationRule & ir);
                                                                        
    void postRegistrationSetup(typename Traits::SetupData d,           
                               PHX::FieldManager<Traits>& fm);        
                                                                     
    void evaluateFields(typename Traits::EvalData d);               


private:
  typedef typename EvalT::ScalarT ScalarT;

  // Simulation source
  PHX::MDField<ScalarT,panzer::Cell,panzer::Point> result;

  double acoeff_;
  double bcoeff_;
  int ir_degree_;
  int ir_index_;
};

}

#endif
