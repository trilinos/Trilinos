#include "Thyra_DefaultSpmdVector_decl.hpp"

#ifdef HAVE_THYRA_EXPLICIT_INSTANTIATION

#include "Thyra_DefaultSpmdVector_def.hpp"
#include "Teuchos_ExplicitInstantiationHelpers.hpp"

namespace Thyra {

TEUCHOS_CLASS_TEMPLATE_INSTANT_SCALAR_TYPES(DefaultSpmdVector)

} // namespace Thyra

#endif // HAVE_THYRA_EXPLICIT_INSTANTIATION
