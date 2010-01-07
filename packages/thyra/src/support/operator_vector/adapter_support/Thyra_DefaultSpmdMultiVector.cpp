#include "Thyra_DefaultSpmdMultiVector_decl.hpp"

#ifdef HAVE_THYRA_EXPLICIT_INSTANTIATION

#include "Thyra_DefaultSpmdMultiVector_def.hpp"
#include "Teuchos_ExplicitInstantiationHelpers.hpp"

namespace Thyra {

TEUCHOS_CLASS_TEMPLATE_INSTANT_SCALAR_TYPES(DefaultSpmdMultiVector)

} // namespace Thyra

#endif // HAVE_THYRA_EXPLICIT_INSTANTIATION
