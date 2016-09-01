#ifndef TEMPUS_EXPLICIT_TEMPLATE_INSTANTIATION_HPP
#define TEMPUS_EXPLICIT_TEMPLATE_INSTANTIATION_HPP

#include "Tempus_config.hpp"

// Always instantiate on double
#define TEMPUS_INSTANTIATE_TEMPLATE_CLASS_ON_DOUBLE(name) \
  template class name<double>; 

// Complex not yet supported. Just need to add define to cmake logic
#ifdef TEMPUS_BUILD_COMPLEX_SUPPORT 
#define TEMPUS_INSTANTIATE_TEMPLATE_CLASS_ON_COMPLEX_DOUBLE(name) \
  template class name<std::complex<double>>;
#else 
#define TEMPUS_INSTANTIATE_TEMPLATE_CLASS_ON_COMPLEX_DOUBLE(name)
#endif

#define TEMPUS_INSTANTIATE_TEMPLATE_CLASS(name) \
  TEMPUS_INSTANTIATE_TEMPLATE_CLASS_ON_DOUBLE(name) \
  TEMPUS_INSTANTIATE_TEMPLATE_CLASS_ON_COMPLEX_DOUBLE(name)

#endif
