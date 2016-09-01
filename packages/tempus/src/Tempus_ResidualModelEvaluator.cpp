#include "Tempus_ExplicitTemplateInstantiation.hpp"

#ifdef HAVE_TEMPUS_EXPLICIT_INSTANTIATION
#include "Tempus_ResidualModelEvaluator.hpp"
#include "Tempus_ResidualModelEvaluator_impl.hpp"

namespace Tempus {
  TEMPUS_INSTANTIATE_TEMPLATE_CLASS(ResidualModelEvaluator)
}

#endif
