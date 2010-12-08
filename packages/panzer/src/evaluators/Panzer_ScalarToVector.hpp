#ifndef PANZER_EVALUATOR_SCALAR_TO_VECTOR_HPP
#define PANZER_EVALUATOR_SCALAR_TO_VECTOR_HPP

#include <vector>
#include "Phalanx_Evaluator_Macros.hpp"
#include "Phalanx_Field.hpp"

namespace panzer {
    
//! Interpolates basis DOF values to IP DOF values
PHX_EVALUATOR_CLASS(ScalarToVector)
  
  std::vector< PHX::MDField<ScalarT,Cell,Point> > scalar_fields;
  PHX::MDField<ScalarT,Cell,Point,Dim> vector_field;

PHX_EVALUATOR_CLASS_END

}

#include "Panzer_ScalarToVectorT.hpp"

#endif
