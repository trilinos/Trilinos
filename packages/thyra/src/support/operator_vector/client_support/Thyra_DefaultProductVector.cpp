#include "Thyra_DefaultProductVector_decl.hpp"

#if (defined(HAVE_THYRA_EXPLICIT_INSTANTIATION) || defined(THYRA_DEFAULT_PRODUCT_VECTOR_EXPLICIT_INSTANTIATION))

#include "Thyra_DefaultProductVector_def.hpp"
#include "Teuchos_ExplicitInstantiationHelpers.hpp"

namespace Thyra {

TEUCHOS_CLASS_TEMPLATE_INSTANT_SCALAR_TYPES(DefaultProductVector)

} // namespace Thyra

#endif // HAVE_THYRA_EXPLICIT_INSTANTIATION
