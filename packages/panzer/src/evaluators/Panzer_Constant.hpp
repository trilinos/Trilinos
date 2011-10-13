#ifndef PANZER_EVALUATOR_CONSTANT_HPP
#define PANZER_EVALUATOR_CONSTANT_HPP

#include "Panzer_config.hpp"

#include "Phalanx_Evaluator_Macros.hpp"
#include "Phalanx_Field.hpp"

namespace panzer {
    
PHX_EVALUATOR_CLASS(Constant)
  
  ScalarT value;
  
  PHX::Field<ScalarT> constant;
  
PHX_EVALUATOR_CLASS_END

}

#ifndef PANZER_EXPLICIT_TEMPLATE_INSTANTIATION
#include "Panzer_ConstantT.hpp"
#endif

#endif
