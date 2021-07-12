#include "Teuchos_ConfigDefs.hpp"
#include "Thyra_ConfigDefs.hpp"

#if defined(HAVE_THYRA_EXPLICIT_INSTANTIATION) && defined(HAVE_TEUCHOS_LONG_DOUBLE)

#include "Thyra_MultiVectorDefaultBase_decl.hpp"
#include "Thyra_MultiVectorDefaultBase_def.hpp"
#include "Teuchos_ExplicitInstantiationHelpers.hpp"

namespace Thyra {

TEUCHOS_CLASS_TEMPLATE_INSTANT_LONG_DOUBLE(MultiVectorDefaultBase)

} // namespace Thyra

#endif // HAVE_THYRA_EXPLICIT_INSTANTIATION && HAVE_TEUCHOS_LONG_DOUBLE
