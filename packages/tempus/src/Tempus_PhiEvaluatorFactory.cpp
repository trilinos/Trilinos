//@HEADER
// *****************************************************************************
// TODO
// *****************************************************************************
//@HEADER

#include "Tempus_ExplicitTemplateInstantiation.hpp"

#ifdef HAVE_TEMPUS_EXPLICIT_INSTANTIATION
#include "Tempus_PhiEvaluatorFactory.hpp"
#include "Tempus_PhiEvaluatorFactory_impl.hpp"

namespace Tempus {
TEMPUS_INSTANTIATE_TEMPLATE_CLASS(PhiEvaluatorFactory)
}

#endif
