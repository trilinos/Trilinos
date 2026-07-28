// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef MINIEM_MAXWELLANALYTICSOLUTION_HPP
#define MINIEM_MAXWELLANALYTICSOLUTION_HPP

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

  /** Computes the manufactured analytic solution
    * E(x,y[,z],t) = sin(timeScale*t) * grad[sin(pi*x)*sin(pi*y)[*sin(pi*z)]]
    * used for verifying the mini-em Maxwell solver.
   */
  template<typename EvalT, typename Traits>
  class MaxwellAnalyticSolution : public panzer::EvaluatorWithBaseImpl<Traits>,
                                  public PHX::EvaluatorDerived<EvalT, Traits>  {

  public:
    /// \brief Construct from the output field name, integration rule, field layout library, time-oscillation scale, and the name of the DOF basis whose integration points to evaluate at.
    MaxwellAnalyticSolution(const std::string & name,
                            const panzer::IntegrationRule & ir,
                            const panzer::FieldLayoutLibrary & fl,
                            const double timeScale,
                            const std::string& basisName="E_edge");

    /// \brief Looks up the integration rule index needed by evaluateFields(). Called once before the first evaluation.
    void postRegistrationSetup(typename Traits::SetupData d,
                               PHX::FieldManager<Traits>& fm);

    /// \brief Computes the analytic solution field at this workset's integration points and time, per the class description.
    void evaluateFields(typename Traits::EvalData d);


  private:
    typedef typename EvalT::ScalarT ScalarT;

    // Simulation source
    PHX::MDField<ScalarT,Cell,Point,Dim> source;
    int ir_degree, ir_index, ir_dim;

    double timeScale_;

    using device_type = PHX::Device;
  };

}

#include "MiniEM_MaxwellAnalyticSolution_impl.hpp"

#endif
