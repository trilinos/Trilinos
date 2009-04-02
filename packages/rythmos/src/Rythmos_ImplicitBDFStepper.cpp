#include "Rythmos_ImplicitBDFStepper_decl.hpp"

#ifdef HAVE_RYTHMOS_EXPLICIT_INSTANTIATION

#include "Rythmos_ImplicitBDFStepper_def.hpp"
#include "Teuchos_ExplicitInstantiationHelpers.hpp"

namespace Rythmos {

TEUCHOS_MACRO_TEMPLATE_INSTANT_DOUBLE(RYTHMOS_IMPLICITBDF_STEPPER_INSTANT) 

} // namespace Rythmos

#endif // HAVE_RYTHMOS_EXPLICIT_INSTANTIATION

