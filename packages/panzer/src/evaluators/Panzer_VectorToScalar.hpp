#ifndef PANZER_EVALUATOR_VECTOR_TO_SCALAR_HPP
#define PANZER_EVALUATOR_VECTOR_TO_SCALAR_HPP

#include <vector>
#include "Phalanx_Evaluator_Macros.hpp"
#include "Phalanx_Field.hpp"

namespace panzer {
    
PHX_EVALUATOR_CLASS(VectorToScalar)
  
  std::vector< PHX::MDField<ScalarT,Cell,Point> > scalar_fields;
  PHX::MDField<ScalarT,Cell,Point,Dim> vector_field;

PHX_EVALUATOR_CLASS_END

}

#include "Panzer_VectorToScalarT.hpp"

#endif
