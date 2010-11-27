#ifndef USER_APP_CONSTANT_MODEL_HPP
#define USER_APP_CONSTANT_MODEL_HPP

#include "Phalanx_Evaluator_Macros.hpp"
#include "Phalanx_Field.hpp"

namespace user_app {
    
PHX_EVALUATOR_CLASS(ConstantModel)
  
  ScalarT value;
  
  PHX::Field<ScalarT> constant;
  
PHX_EVALUATOR_CLASS_END

}

#include "user_app_ConstantModelT.hpp"

#endif
