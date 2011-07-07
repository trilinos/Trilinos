#include "Thyra_describeLinearOp_decl.hpp"

#ifdef HAVE_THYRA_EXPLICIT_INSTANTIATION

#include "Thyra_describeLinearOp_def.hpp"
#include "Teuchos_ExplicitInstantiationHelpers.hpp"

namespace Thyra {

TEUCHOS_MACRO_TEMPLATE_INSTANT_SCALAR_TYPES(THYRA_DESCRIBE_LINEAR_INSTANT)

} // namespace Thyra

#endif // HAVE_THYRA_EXPLICIT_INSTANTIATION
