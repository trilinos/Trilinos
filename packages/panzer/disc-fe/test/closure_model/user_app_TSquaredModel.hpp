// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef USER_APP_TSQUARED_MODEL_HPP
#define USER_APP_TSQUARED_MODEL_HPP

#include "Phalanx_Evaluator_Macros.hpp"
#include "Phalanx_MDField.hpp"

namespace user_app {

  /// Example evaluator that computes a source term k*T^2 at the integration points.
  template<typename EvalT, typename Traits>
  class TSquaredModel : public PHX::EvaluatorWithBaseImpl<Traits>,
                        public PHX::EvaluatorDerived<EvalT, Traits>
  {
  public:
    TSquaredModel(const Teuchos::ParameterList& p);
    void evaluateFields(typename Traits::EvalData d);

  private:
    using ScalarT = typename EvalT::ScalarT;
    double k_;
    PHX::MDField<const ScalarT,panzer::Cell,panzer::IP> t_;
    PHX::MDField<ScalarT,panzer::Cell,panzer::IP> k_t_squared_;
  };

} // namespace user_app

#include "user_app_TSquaredModel_impl.hpp"

#endif
