#ifndef USER_APP_CONVECTION_HPP
#define USER_APP_CONVECTION_HPP

#include "Panzer_Dimension.hpp"
#include "Phalanx_Evaluator_Macros.hpp"
#include "Phalanx_MDField.hpp"

using panzer::Cell;
using panzer::Point;
using panzer::Dim;

namespace user_app {
    
  //! Convection: conv = multiplier * a \cdot \nabla x
PHX_EVALUATOR_CLASS(Convection)

  PHX::MDField<ScalarT,Cell,Point> conv;
    
  PHX::MDField<ScalarT,Cell,Point,Dim> a;
    
  PHX::MDField<ScalarT,Cell,Point,Dim> grad_x;
    
  double multiplier;

PHX_EVALUATOR_CLASS_END

}

#include "user_app_ConvectionT.hpp"

#endif
