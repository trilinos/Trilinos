// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef MINIEM_DARCYANALYTICSOLUTION_HPP
#define MINIEM_DARCYANALYTICSOLUTION_HPP

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
  * u(x,y[,z],t) = sin(pi*t)*sin(pi*x)*sin(pi*y)[*sin(pi*z)] used for
  * verifying the mini-em Darcy solver.
  */
template<typename EvalT, typename Traits>
class DarcyAnalyticSolution : public panzer::EvaluatorWithBaseImpl<Traits>,
                      public PHX::EvaluatorDerived<EvalT, Traits>  {

public:
    /// \brief Construct from the output field name, integration rule, field layout library, and the name of the DOF basis whose integration points to evaluate at.
    DarcyAnalyticSolution(const std::string & name,
                  const panzer::IntegrationRule & ir,
                  const panzer::FieldLayoutLibrary & fl,
                  const std::string& basisName="u");

    /// \brief Looks up the integration rule index needed by evaluateFields(). Called once before the first evaluation.
    void postRegistrationSetup(typename Traits::SetupData d,
                               PHX::FieldManager<Traits>& fm);

    /// \brief Computes the analytic solution field at this workset's integration points and time, per the class description.
    void evaluateFields(typename Traits::EvalData d);


private:
  typedef typename EvalT::ScalarT ScalarT;

  // Simulation source
  PHX::MDField<ScalarT,Cell,Point> source;
  int ir_degree, ir_index, ir_dim;

  using device_type = PHX::Device;
};

}

#include "MiniEM_DarcyAnalyticSolution_impl.hpp"

#endif
