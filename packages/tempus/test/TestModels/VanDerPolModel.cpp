#include "Tempus_ExplicitTemplateInstantiation.hpp"

#ifdef HAVE_TEMPUS_EXPLICIT_INSTANTIATION
#include "VanDerPolModel.hpp"
#include "VanDerPolModel_impl.hpp"

namespace Tempus_Test {
  TEMPUS_INSTANTIATE_TEMPLATE_CLASS(VanDerPolModel)
}

#endif
