// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef MINIEM_PIECEWISECONSTANT_DECL_HPP
#define MINIEM_PIECEWISECONSTANT_DECL_HPP

#include "PanzerAdaptersSTK_config.hpp"

#include "Phalanx_config.hpp"
#include "Phalanx_Evaluator_WithBaseImpl.hpp"
#include "Phalanx_Evaluator_Derived.hpp"
#include "Phalanx_FieldManager.hpp"

#include "Panzer_FieldLibrary.hpp"

#include <string>

#include "Panzer_Evaluator_WithBaseImpl.hpp"

namespace mini_em {

  using panzer::Cell;
  using panzer::Point;

  template<typename EvalT, typename Traits>
  class PiecewiseConstant : public panzer::EvaluatorWithBaseImpl<Traits>,
                       public PHX::EvaluatorDerived<EvalT, Traits>  {

  public:
    PiecewiseConstant(const std::string & name,
                      const panzer::IntegrationRule & ir,
                      const panzer::FieldLayoutLibrary & fl,
                      const double value0,
                      const double value1,
                      const double xl,
                      const double xr,
                      const double yl,
                      const double yr,
                      const double zl,
                      const double zr,
                      const std::string& DoF_);

    void evaluateFields(typename Traits::EvalData d);


  private:
    typedef typename EvalT::ScalarT ScalarT;

    PHX::MDField<ScalarT,Cell,Point> values;
    PHX::MDField<const ScalarT,Cell,Point,Dim> coords;
    int ir_degree, ir_dim;
    double value0_, value1_;
    double xl_, xr_, yl_, yr_, zl_, zr_;
  };

}

#include "MiniEM_PiecewiseConstant_impl.hpp"

#endif
