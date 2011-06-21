#ifndef PANZER_EVALUATOR_DIRICHLET_RESIDUAL_HPP
#define PANZER_EVALUATOR_DIRICHLET_RESIDUAL_HPP

#include "Phalanx_Evaluator_Macros.hpp"
#include "Phalanx_MDField.hpp"

namespace panzer {
    
//! Evaluates a Dirichlet BC residual corresponding to a field value
PHX_EVALUATOR_CLASS(DirichletResidual)
  
  PHX::MDField<ScalarT> residual;
  PHX::MDField<ScalarT> dof;
  PHX::MDField<ScalarT> value;

  std::size_t cell_data_size;

PHX_EVALUATOR_CLASS_END

}

#include "Panzer_Dirichlet_ResidualT.hpp"

#endif
