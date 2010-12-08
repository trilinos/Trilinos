#ifndef PANZER_EVALUATOR_SUM_HPP
#define PANZER_EVALUATOR_SUM_HPP

#include "Phalanx_Evaluator_Macros.hpp"
#include "Phalanx_MDField.hpp"

namespace panzer {
    
//! Sums entries on a single data layout
PHX_EVALUATOR_CLASS(Sum)
  
  PHX::MDField<ScalarT> sum;
  std::vector< PHX::MDField<ScalarT> > values;

  std::size_t cell_data_size;

PHX_EVALUATOR_CLASS_END

}

#include "Panzer_SumT.hpp"

#endif
