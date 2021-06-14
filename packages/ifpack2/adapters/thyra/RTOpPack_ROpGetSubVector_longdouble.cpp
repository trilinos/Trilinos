#include "Teuchos_ConfigDefs.hpp"
#include "RTOpPack_ROpGetSubVector_decl.hpp"

#if defined(HAVE_RTOP_EXPLICIT_INSTANTIATION) && defined(HAVE_TEUCHOS_LONG_DOUBLE)

#include "RTOpPack_ROpGetSubVector_def.hpp"
#include "Teuchos_ExplicitInstantiationHelpers.hpp"

namespace RTOpPack {

TEUCHOS_CLASS_TEMPLATE_INSTANT_LONG_DOUBLE(ROpGetSubVector)

} // namespace RTOpPack

#endif // HAVE_RTOP_EXCPLICIT_INSTANTIATION && HAVE_TEUCHOS_LONG_DOUBLE
