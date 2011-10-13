#ifndef PANZER_EXPLICIT_TEMPLATE_INSTANTIATION_HPP
#define PANZER_EXPLICIT_TEMPLATE_INSTANTIATION_HPP

#include "Panzer_Traits.hpp"

// ONE template argument 
#define PANZER_INSTANTIATE_TEMPLATE_CLASS_RESIDUAL_ONE_T(name) \
  template class name<panzer::Traits::Residual>; 

#define PANZER_INSTANTIATE_TEMPLATE_CLASS_JACOBIAN_ONE_T(name) \
  template class name<panzer::Traits::Jacobian>; 

// stochastic galerkin objects
#ifdef HAVE_STOKHOS
   #define PANZER_INSTANTIATE_TEMPLATE_CLASS_SGRESIDUAL_ONE_T(name) \
     template class name<panzer::Traits::SGResidual>; 

   #define PANZER_INSTANTIATE_TEMPLATE_CLASS_SGJACOBIAN_ONE_T(name) \
     template class name<panzer::Traits::SGJacobian>; 
#else
   #define PANZER_INSTANTIATE_TEMPLATE_CLASS_SGRESIDUAL_ONE_T(name)
   #define PANZER_INSTANTIATE_TEMPLATE_CLASS_SGJACOBIAN_ONE_T(name)
#endif

#define PANZER_INSTANTIATE_TEMPLATE_CLASS_ONE_T(name) \
  PANZER_INSTANTIATE_TEMPLATE_CLASS_RESIDUAL_ONE_T(name) \
  PANZER_INSTANTIATE_TEMPLATE_CLASS_JACOBIAN_ONE_T(name) \
  PANZER_INSTANTIATE_TEMPLATE_CLASS_SGRESIDUAL_ONE_T(name) \
  PANZER_INSTANTIATE_TEMPLATE_CLASS_SGJACOBIAN_ONE_T(name)

// TWO template arguments
#define PANZER_INSTANTIATE_TEMPLATE_CLASS_RESIDUAL_TWO_T(name) \
  template class name<panzer::Traits::Residual, panzer::Traits>; 

#define PANZER_INSTANTIATE_TEMPLATE_CLASS_JACOBIAN_TWO_T(name) \
  template class name<panzer::Traits::Jacobian, panzer::Traits>; 

// stochastic galerkin objects
#ifdef HAVE_STOKHOS
   #define PANZER_INSTANTIATE_TEMPLATE_CLASS_SGRESIDUAL_TWO_T(name) \
     template class name<panzer::Traits::SGResidual, panzer::Traits>; 
   #define PANZER_INSTANTIATE_TEMPLATE_CLASS_SGJACOBIAN_TWO_T(name) \
     template class name<panzer::Traits::SGJacobian, panzer::Traits>; 
#else
   #define PANZER_INSTANTIATE_TEMPLATE_CLASS_SGRESIDUAL_TWO_T(name) 
   #define PANZER_INSTANTIATE_TEMPLATE_CLASS_SGJACOBIAN_TWO_T(name) 
#endif

#define PANZER_INSTANTIATE_TEMPLATE_CLASS_TWO_T(name) \
  PANZER_INSTANTIATE_TEMPLATE_CLASS_RESIDUAL_TWO_T(name) \
  PANZER_INSTANTIATE_TEMPLATE_CLASS_JACOBIAN_TWO_T(name) \
  PANZER_INSTANTIATE_TEMPLATE_CLASS_SGRESIDUAL_TWO_T(name) \
  PANZER_INSTANTIATE_TEMPLATE_CLASS_SGJACOBIAN_TWO_T(name)


#endif
