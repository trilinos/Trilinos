// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef MINIEM_RTC_CONDUCTIVITY_HPP
#define MINIEM_RTC_CONDUCTIVITY_HPP

#include "PanzerAdaptersSTK_config.hpp"

#include "Phalanx_config.hpp"
#include "Phalanx_Evaluator_WithBaseImpl.hpp"
#include "Phalanx_Evaluator_Derived.hpp"
#include "Phalanx_FieldManager.hpp"

#include "Panzer_FieldLibrary.hpp"

#include <string>

#include "Panzer_Evaluator_WithBaseImpl.hpp"

#include "RTC_FunctionRTC.hh"

namespace mini_em {

  using panzer::Cell;
  using panzer::Point;

  /** \brief Evaluates a scalar field at each point using a
    * user-provided function body, compiled and executed at run time
    * (via the RTC "RunTime Compiler" library) as a function of the
    * point's (x,y[,z]) coordinates.
    */
  template<typename EvalT, typename Traits>
  class RTC : public panzer::EvaluatorWithBaseImpl<Traits>,
              public PHX::EvaluatorDerived<EvalT, Traits>  {

  public:
    /// \brief Constructor.
    RTC(const std::string & name,
        const panzer::IntegrationRule & ir,
        const panzer::FieldLayoutLibrary & fl,
        const std::string &funBody,
        const std::string& DoF_);

    /// \brief Evaluates the run-time-compiled function at this workset's points per the class description.
    void evaluateFields(typename Traits::EvalData d);


  private:
    using ScalarT = typename EvalT::ScalarT;

    PHX::MDField<ScalarT,Cell,Point> values;
    PHX::MDField<const ScalarT,Cell,Point,Dim> coords;
    int ir_degree, ir_dim;
    PG_RuntimeCompiler::Function fun_;
  };

}

#include "MiniEM_RTC_impl.hpp"

#endif
