#include "Thyra_DefaultProductMultiVector_decl.hpp"

#if defined(HAVE_THYRA_EXPLICIT_INSTANTIATION) || defined(THYRA_DEFAULT_PRODUCT_MULTI_VECTOR_EXPLICIT_INSTANTIATION)

#include "Thyra_DefaultProductMultiVector_def.hpp"
#include "Teuchos_ExplicitInstantiationHelpers.hpp"

namespace Thyra {

TEUCHOS_MACRO_TEMPLATE_INSTANT_SCALAR_TYPES(THYRA_DEFAULT_PRODUCT_MULTI_VECTOR_INSTANT)

} // namespace Thyra

#endif // HAVE_THYRA_EXPLICIT_INSTANTIATION
