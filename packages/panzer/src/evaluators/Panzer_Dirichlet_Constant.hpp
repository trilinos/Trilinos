#ifndef PANZER_EVALUATOR_DIRICHLET_CONSTANT_HPP
#define PANZER_EVALUATOR_DIRICHLET_CONSTANT_HPP

#include "Phalanx_Evaluator_Macros.hpp"
#include "Phalanx_MDField.hpp"

namespace panzer {
    
//! Sums entries on a single data layout
PHX_EVALUATOR_CLASS(DirichletConstant)
  
  PHX::MDField<ScalarT> residual;
  PHX::MDField<ScalarT> dof;
  double value;
  std::size_t cell_data_size;

PHX_EVALUATOR_CLASS_END

}

#include "Panzer_Dirichlet_ConstantT.hpp"

#endif
