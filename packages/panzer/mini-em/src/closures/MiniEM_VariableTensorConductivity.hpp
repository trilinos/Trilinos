// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef MINIEM_VARIABLETENSORCONDUCTIVITY_DECL_HPP
#define MINIEM_VARIABLETENSORCONDUCTIVITY_DECL_HPP

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
  class VariableTensorConductivity : public panzer::EvaluatorWithBaseImpl<Traits>,
                       public PHX::EvaluatorDerived<EvalT, Traits>  {

  public:
    VariableTensorConductivity(const std::string & name,
                               const panzer::IntegrationRule & ir,
                               const panzer::FieldLayoutLibrary & fl,
                               const double & sigma0_,
                               const double & sigma1_,
                               const double & sigma2_,
                               const double & betax0_,
                               const double & betay0_,
                               const double & betaz0_,
                               const double & betax1_,
                               const double & betay1_,
                               const double & betaz1_,
                               const double & betax2_,
                               const double & betay2_,
                               const double & betaz2_,
                               const std::string& DoF_);

    void postRegistrationSetup(typename Traits::SetupData d,
                               PHX::FieldManager<Traits>& fm);

    void evaluateFields(typename Traits::EvalData d);


  private:
    typedef typename EvalT::ScalarT ScalarT;

    PHX::MDField<ScalarT,Cell,Point,Dim,Dim> conductivity;
    PHX::MDField<const ScalarT,Cell,Point,Dim> coords;
    int ir_degree, ir_dim, ir_index;
    double sigma0, betax0, betay0, betaz0;
    double sigma1, betax1, betay1, betaz1;
    double sigma2, betax2, betay2, betaz2;
  };

}

#include "MiniEM_VariableTensorConductivity_impl.hpp"

#endif
