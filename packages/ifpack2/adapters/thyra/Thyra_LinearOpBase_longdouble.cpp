#include "Teuchos_ConfigDefs.hpp"
#include "Thyra_ConfigDefs.hpp"

#if defined(HAVE_THYRA_EXPLICIT_INSTANTIATION) && defined(HAVE_TEUCHOS_LONG_DOUBLE)

#include "Thyra_LinearOpBase_decl.hpp"
#include "Thyra_LinearOpBase_def.hpp"
#include "Teuchos_ExplicitInstantiationHelpers.hpp"

namespace Thyra {

// This instantiates both LinearOpBase<__float128> and
// apply(const LinearOpBase<__float128>&, ...).
TEUCHOS_MACRO_TEMPLATE_INSTANT_LONG_DOUBLE(THYRA_LINEAR_OP_BASE_INSTANT)

} // namespace Thyra

#endif // HAVE_THYRA_EXPLICIT_INSTANTIATION && HAVE_TEUCHOS_LONG_DOUBLE
