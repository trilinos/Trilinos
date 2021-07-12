#include "Teuchos_ConfigDefs.hpp"
#include "Thyra_ConfigDefs.hpp"

#if defined(HAVE_THYRA_EXPLICIT_INSTANTIATION) && defined(HAVE_TEUCHOS_LONG_DOUBLE)

#include "Thyra_SpmdVectorSpaceDefaultBase_decl.hpp"
#include "Thyra_SpmdVectorSpaceDefaultBase_def.hpp"
#include "Teuchos_ExplicitInstantiationHelpers.hpp"

namespace Thyra {

TEUCHOS_CLASS_TEMPLATE_INSTANT_LONG_DOUBLE(SpmdVectorSpaceDefaultBase)

} // namespace Thyra

#endif // HAVE_THYRA_EXPLICIT_INSTANTIATION && HAVE_TEUCHOS_LONG_DOUBLE
