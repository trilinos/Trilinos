// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef MINIEM_TENSORCONDUCTIVITY_DECL_HPP
#define MINIEM_TENSORCONDUCTIVITY_DECL_HPP

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
  class TensorConductivity : public panzer::EvaluatorWithBaseImpl<Traits>,
                       public PHX::EvaluatorDerived<EvalT, Traits>  {

  public:
    TensorConductivity(const std::string & name,
                 const panzer::IntegrationRule & ir,
                 const panzer::FieldLayoutLibrary & fl,
                 const double & sigma_,
                 const double & betax_,
                 const double & betay_,
                 const double & betaz_,
                 const std::string& DoF_);

    void evaluateFields(typename Traits::EvalData d);


  private:
    typedef typename EvalT::ScalarT ScalarT;

    PHX::MDField<ScalarT,Cell,Point,Dim,Dim> conductivity;
    PHX::MDField<const ScalarT,Cell,Point,Dim> coords;
    int ir_degree, ir_dim;
    double sigma, betax, betay, betaz;
  };

}

#include "MiniEM_TensorConductivity_impl.hpp"

#endif
