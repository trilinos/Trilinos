#ifndef PANZER_EXPLICIT_TEMPLATE_INSTANTIATION_HPP
#define PANZER_EXPLICIT_TEMPLATE_INSTANTIATION_HPP

#include "Panzer_Traits.hpp"

// ONE template argument 
#define PANZER_INSTANTIATE_TEMPLATE_CLASS_RESIDUAL_ONE_T(name) \
  template class name<panzer::Traits::Residual>; 

#define PANZER_INSTANTIATE_TEMPLATE_CLASS_JACOBIAN_ONE_T(name) \
  template class name<panzer::Traits::Jacobian>; 

#define PANZER_INSTANTIATE_TEMPLATE_CLASS_ONE_T(name) \
  PANZER_INSTANTIATE_TEMPLATE_CLASS_RESIDUAL_ONE_T(name) \
  PANZER_INSTANTIATE_TEMPLATE_CLASS_JACOBIAN_ONE_T(name)

// TWO template arguments
#define PANZER_INSTANTIATE_TEMPLATE_CLASS_RESIDUAL_TWO_T(name) \
  template class name<panzer::Traits::Residual, panzer::Traits>; 

#define PANZER_INSTANTIATE_TEMPLATE_CLASS_JACOBIAN_TWO_T(name) \
  template class name<panzer::Traits::Jacobian, panzer::Traits>; 

#define PANZER_INSTANTIATE_TEMPLATE_CLASS_TWO_T(name) \
  PANZER_INSTANTIATE_TEMPLATE_CLASS_RESIDUAL_TWO_T(name) \
  PANZER_INSTANTIATE_TEMPLATE_CLASS_JACOBIAN_TWO_T(name)

#endif
