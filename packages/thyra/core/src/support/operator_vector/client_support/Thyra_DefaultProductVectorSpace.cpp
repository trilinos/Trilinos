#include "Thyra_DefaultProductVectorSpace_decl.hpp"

#if (defined(HAVE_THYRA_EXPLICIT_INSTANTIATION) || defined(THYRA_DEFAULT_PRODUCT_VECTOR_SPACE_EXPLICIT_INSTANTIATION))

#include "Thyra_DefaultProductVectorSpace_def.hpp"
#include "Teuchos_ExplicitInstantiationHelpers.hpp"

namespace Thyra {

TEUCHOS_CLASS_TEMPLATE_INSTANT_SCALAR_TYPES(DefaultProductVectorSpace)

} // namespace Thyra

#endif // HAVE_THYRA_EXPLICIT_INSTANTIATION
