#include "Tempus_ExplicitTemplateInstantiation.hpp"

#ifdef HAVE_TEMPUS_EXPLICIT_INSTANTIATION
#include "CDR_Model.hpp"
#include "CDR_Model_impl.hpp"

namespace Tempus_Test {
  TEMPUS_INSTANTIATE_TEMPLATE_CLASS(CDR_Model)
}

#endif
