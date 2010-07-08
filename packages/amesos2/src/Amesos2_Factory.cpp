#ifndef AMESOS2_FACTORY_DECL_HPP
#define AMESOS2_FACTORY_DECL_HPP

#include "Amesos2_Factory_decl.hpp"

#ifdef HAVE_AMESOS2_EXPLICIT_INSTANTIATION
#  include "Amesos2_Factory_def.hpp"
#  include "Teuchos_ExplicitInstantiationHelpers.hpp"
// I don't think this makes sense to do for a class with only static methods...
// namespace Amesos {
// TEUCHOS_CLASS_TEMPLATE_INSTANT_SCALAR_TYPES(Amesos::Factory)
// }
#endif  // HAVE_AMESOS2_EXPLICIT_INSTANTIATION

#endif  // AMESOS2_FACTORY_DECL_HPP
