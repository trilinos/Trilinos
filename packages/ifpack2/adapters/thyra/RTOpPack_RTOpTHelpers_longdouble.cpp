#include "Teuchos_ConfigDefs.hpp"
#include "RTOpPack_RTOpTHelpers_decl.hpp"

#if defined(HAVE_RTOP_EXPLICIT_INSTANTIATION) && defined(HAVE_TEUCHOS_LONG_DOUBLE)

#include "RTOpPack_RTOpTHelpers_def.hpp"
#include "Teuchos_ExplicitInstantiationHelpers.hpp"

namespace RTOpPack {

  TEUCHOS_MACRO_TEMPLATE_INSTANT_LONG_DOUBLE(RTOPPACK_RTOPT_HELPERS_INSTANT_SCALAR)

} // namespace RTOpPack

#endif // HAVE_RTOP_EXPLICIT_INSTANTIATION && HAVE_TEUCHOS_LONG_DOUBLE
