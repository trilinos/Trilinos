#include "Tempus_ExplicitTemplateInstantiation.hpp"

#ifdef HAVE_TEMPUS_EXPLICIT_INSTANTIATION
#include "SinCosModel.hpp"
#include "SinCosModel_impl.hpp"

namespace Tempus_Test {
  TEMPUS_INSTANTIATE_TEMPLATE_CLASS(SinCosModel)
}

#endif
