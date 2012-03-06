#ifndef PANZER_EVALUATOR_DOF_CURL_DECL_HPP
#define PANZER_EVALUATOR_DOF_CURL_DECL_HPP

#include "Phalanx_Evaluator_Macros.hpp"
#include "Phalanx_Field.hpp"

namespace panzer {
    
//! Interpolates basis DOF values to IP DOF Curl values
PHX_EVALUATOR_CLASS(DOFCurl)
  
  PHX::MDField<ScalarT,Cell,Point> dof_value;
  PHX::MDField<ScalarT> dof_curl;

  std::string basis_name;
  std::size_t basis_index;

  PHX::MDField<ScalarT,Cell,BASIS> dof_orientation;

PHX_EVALUATOR_CLASS_END

}

#endif
