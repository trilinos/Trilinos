// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef USER_APP_CONVECTION_HPP
#define USER_APP_CONVECTION_HPP

#include "Panzer_Dimension.hpp"
#include "Phalanx_Evaluator_Macros.hpp"
#include "Phalanx_MDField.hpp"

#include "Panzer_Evaluator_Macros.hpp"

using panzer::Cell;
using panzer::Point;
using panzer::Dim;

namespace user_app {
    
  //! Convection: conv = multiplier * a \cdot \nabla x
template<typename EvalT, typename Traits>
class Convection
  :
  public panzer::EvaluatorWithBaseImpl<Traits>,
  public PHX::EvaluatorDerived<EvalT, Traits>
{
  public:

    Convection(
      const Teuchos::ParameterList& p);

    void
    evaluateFields(
      typename Traits::EvalData d);

  private:

    using ScalarT = typename EvalT::ScalarT;

  PHX::MDField<ScalarT,Cell,Point> conv;
    
  PHX::MDField<const ScalarT,Cell,Point,Dim> a;
    
  PHX::MDField<const ScalarT,Cell,Point,Dim> grad_x;
    
  double multiplier;

}; // end of class Convection


}

#include "user_app_Convection_impl.hpp"

#endif
