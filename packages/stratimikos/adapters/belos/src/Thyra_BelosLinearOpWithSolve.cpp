#include "Thyra_BelosLinearOpWithSolve_decl.hpp"

#ifdef HAVE_THYRA_EXPLICIT_INSTANTIATION

#include "Thyra_BelosLinearOpWithSolve_def.hpp"
#include "Teuchos_ExplicitInstantiationHelpers.hpp"

namespace Thyra {

TEUCHOS_CLASS_TEMPLATE_INSTANT_REAL_SCALAR_TYPES(BelosLinearOpWithSolve)

} // namespace Thyra

#endif // HAVE_THYRA_EXPLICIT_INSTANTIATION
