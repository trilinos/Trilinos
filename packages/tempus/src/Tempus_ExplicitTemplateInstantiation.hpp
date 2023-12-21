// @HEADER
// ****************************************************************************
//                Tempus: Copyright (2017) Sandia Corporation
//
// Distributed under BSD 3-clause license (See accompanying file Copyright.txt)
// ****************************************************************************
// @HEADER

#ifndef Tempus_ExplicitTemplateInstantiation_hpp
#define Tempus_ExplicitTemplateInstantiation_hpp

#include "Tempus_config.hpp"

#define TEMPUS_INSTANTIATE_TEMPLATE_CLASS_TPETRA(name, SC, LO, GO, Node) \
  template class name<SC, LO, GO, Node>;

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
