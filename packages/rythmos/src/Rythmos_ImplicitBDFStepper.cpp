#include "Rythmos_ImplicitBDFStepper_decl.hpp"

#ifdef HAVE_RYTHMOS_EXPLICIT_INSTANTIATION

#include "Rythmos_ImplicitBDFStepper_def.hpp"
#include "Rythmos_ExplicitInstantiationHelpers.hpp"

namespace Rythmos {

RYTHMOS_MACRO_TEMPLATE_INSTANT_SCALAR_TYPES(RYTHMOS_IMPLICITBDF_STEPPER_INSTANT) 

} // namespace Rythmos

#endif // HAVE_RYTHMOS_EXPLICIT_INSTANTIATION

