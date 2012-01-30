#ifndef PANZER_EVALUATOR_DOF_GRADIENT_DECL_HPP
#define PANZER_EVALUATOR_DOF_GRADIENT_DECL_HPP

#include "Phalanx_Evaluator_Macros.hpp"
#include "Phalanx_Field.hpp"

namespace panzer {
    
//! Interpolates basis DOF values to IP DOF Gradient values
PHX_EVALUATOR_CLASS(DOFGradient)
  
  PHX::MDField<ScalarT,Cell,Point> dof_value;
  PHX::MDField<ScalarT,Cell,Point,Dim> dof_gradient;

  std::string basis_name;
  std::size_t basis_index;

PHX_EVALUATOR_CLASS_END

}

#endif
