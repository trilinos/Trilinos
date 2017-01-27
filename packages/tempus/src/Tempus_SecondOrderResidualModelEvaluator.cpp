#include "Tempus_ExplicitTemplateInstantiation.hpp"

#ifdef HAVE_TEMPUS_EXPLICIT_INSTANTIATION
#include "Tempus_SecondOrderResidualModelEvaluator.hpp"
#include "Tempus_SecondOrderResidualModelEvaluator_impl.hpp"

namespace Tempus {
  TEMPUS_INSTANTIATE_TEMPLATE_CLASS(SecondOrderResidualModelEvaluator)
}

#endif
