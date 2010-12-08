#ifndef PANZER_EVALUATOR_DOF_HPP
#define PANZER_EVALUATOR_DOF_HPP

#include <string>
#include "Phalanx_Evaluator_Macros.hpp"
#include "Phalanx_Field.hpp"

namespace panzer {
    
//! Interpolates basis DOF values to IP DOF values
PHX_EVALUATOR_CLASS(DOF)
  
  PHX::MDField<ScalarT,Cell,Point> dof_basis;
  PHX::MDField<ScalarT,Cell,Point> dof_ip;

  std::size_t cell_data_size;

  std::string basis_name;
  std::size_t basis_index;

PHX_EVALUATOR_CLASS_END

}

#include "Panzer_DOFT.hpp"

#endif
