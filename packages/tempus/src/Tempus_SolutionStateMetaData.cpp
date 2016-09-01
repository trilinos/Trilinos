#include "Tempus_ExplicitTemplateInstantiation.hpp"

#ifdef HAVE_TEMPUS_EXPLICIT_INSTANTIATION
#include "Tempus_SolutionStateMetaData.hpp"
#include "Tempus_SolutionStateMetaData_impl.hpp"

namespace Tempus {
  TEMPUS_INSTANTIATE_TEMPLATE_CLASS(SolutionStateMetaData)
}

#endif
