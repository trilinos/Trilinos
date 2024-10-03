// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef MINIEM_MAXWELLANALYTICFORCING_HPP
#define MINIEM_MAXWELLANALYTICFORCING_HPP

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

  /** Maxwell source with know analytic solution
   */
  template<typename EvalT, typename Traits>
  class MaxwellAnalyticForcing : public panzer::EvaluatorWithBaseImpl<Traits>,
                                 public PHX::EvaluatorDerived<EvalT, Traits>  {

  public:
    MaxwellAnalyticForcing(const std::string & name,
                           const panzer::IntegrationRule & ir,
                           const panzer::FieldLayoutLibrary & fl,
                           const double epsilon,
                           const double timeScale,
                           const std::string& basisName="E_edge");

    void postRegistrationSetup(typename Traits::SetupData d,
                               PHX::FieldManager<Traits>& fm);

    void evaluateFields(typename Traits::EvalData d);


  private:
    typedef typename EvalT::ScalarT ScalarT;

    // Simulation source
    PHX::MDField<ScalarT,Cell,Point,Dim> source;
    int ir_degree, ir_index, ir_dim;
    double epsilon_;
    double timeScale_;

    using device_type = PHX::Device;
  };

}

#include "MiniEM_MaxwellAnalyticForcing_impl.hpp"

#endif
