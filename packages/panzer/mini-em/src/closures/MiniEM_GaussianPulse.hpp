// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef MINIEM_GAUSSIAN_PULSE_DECL_HPP
#define MINIEM_GAUSSIAN_PULSE_DECL_HPP

#include "PanzerAdaptersSTK_config.hpp"

#include "Phalanx_config.hpp"
#include "Phalanx_Evaluator_WithBaseImpl.hpp"
#include "Phalanx_Evaluator_Derived.hpp"
#include "Phalanx_FieldManager.hpp"

#include "Panzer_Dimension.hpp"
#include "Panzer_FieldLibrary.hpp"

#include <string>

#include "Panzer_Evaluator_WithBaseImpl.hpp"

namespace mini_em {
    
  using panzer::Cell;
  using panzer::Point;
  using panzer::Dim;

/** Gaussian pulse current source
  * given dt, time width is defined to be resolved by dt
  * Jz = exp(-r^2/alpha^2)*exp(-(t-2*beta)^2/beta^2)
  * dt = beta/5
  */
template<typename EvalT, typename Traits>
class GaussianPulse : public panzer::EvaluatorWithBaseImpl<Traits>,
                      public PHX::EvaluatorDerived<EvalT, Traits>  {

public:
    GaussianPulse(const std::string & name,
                  const panzer::IntegrationRule & ir,
                  const panzer::FieldLayoutLibrary & fl,
                  const double & dt);

    void postRegistrationSetup(typename Traits::SetupData d,
                               PHX::FieldManager<Traits>& fm);
                                                                        
    void evaluateFields(typename Traits::EvalData d);               


private:
  typedef typename EvalT::ScalarT ScalarT;

  // Simulation source
  PHX::MDField<ScalarT,Cell,Point,Dim> current;
  int ir_degree, ir_index, ir_dim;
  double alpha, beta;
};

}

#include "MiniEM_GaussianPulse_impl.hpp"

#endif
