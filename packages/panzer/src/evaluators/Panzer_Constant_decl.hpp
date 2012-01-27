#ifndef PANZER_EVALUATOR_CONSTANT_DECL_HPP
#define PANZER_EVALUATOR_CONSTANT_DECL_HPP

#include "Panzer_config.hpp"

#include "Phalanx_Evaluator_Macros.hpp"
#include "Phalanx_Field.hpp"

namespace panzer {
    
PHX_EVALUATOR_CLASS(Constant)
  
  ScalarT value;
  
  PHX::Field<ScalarT> constant;
  
PHX_EVALUATOR_CLASS_END

}

#endif
