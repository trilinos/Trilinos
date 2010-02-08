#include "Thyra_DefaultAddedLinearOp_decl.hpp"

#ifdef HAVE_THYRA_EXPLICIT_INSTANTIATION

#include "Thyra_DefaultAddedLinearOp_def.hpp"
#include "Teuchos_ExplicitInstantiationHelpers.hpp"

namespace Thyra {

TEUCHOS_MACRO_TEMPLATE_INSTANT_SCALAR_TYPES(THYRA_DEFAULT_ADDED_LINEAR_OP_INSTANT)
//TEUCHOS_CLASS_TEMPLATE_INSTANT_SCALAR_TYPES(DefaultAddedLinearOp)

} // namespace Thyra

#endif // HAVE_THYRA_EXPLICIT_INSTANTIATION
