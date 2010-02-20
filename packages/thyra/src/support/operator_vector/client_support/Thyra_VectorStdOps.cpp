#include "Thyra_VectorStdOps_decl.hpp"

#ifdef HAVE_THYRA_EXPLICIT_INSTANTIATION

#include "Thyra_VectorStdOps_def.hpp"
#include "Teuchos_ExplicitInstantiationHelpers.hpp"

namespace Thyra {

TEUCHOS_MACRO_TEMPLATE_INSTANT_SCALAR_TYPES(THYRA_VECTOR_STD_OPS_INSTANT)

TEUCHOS_MACRO_TEMPLATE_INSTANT_REAL_SCALAR_TYPES(THYRA_VECTOR_STD_OPS_REAL_INSTANT)

} // namespace Thyra

#endif // HAVE_THYRA_EXPLICIT_INSTANTIATION
