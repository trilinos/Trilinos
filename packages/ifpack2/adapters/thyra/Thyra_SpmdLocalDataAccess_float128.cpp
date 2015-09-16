#include "Teuchos_ConfigDefs.hpp"
#include "Thyra_ConfigDefs.hpp"

#if defined(HAVE_THYRA_EXPLICIT_INSTANTIATION) && defined(HAVE_TEUCHOSCORE_QUADMATH)

#include "Thyra_SpmdLocalDataAccess_decl.hpp"
#include "Thyra_SpmdLocalDataAccess_def.hpp"
#include "Teuchos_ExplicitInstantiationHelpers.hpp"

namespace Thyra {

  TEUCHOS_MACRO_TEMPLATE_INSTANT_FLOAT128(THYRA_SPMD_LOCAL_DATA_ACCESS_INSTANT)

} // namespace Thyra

#endif // HAVE_THYRA_EXPLICIT_INSTANTIATION && HAVE_TEUCHOSCORE_QUADMATH
