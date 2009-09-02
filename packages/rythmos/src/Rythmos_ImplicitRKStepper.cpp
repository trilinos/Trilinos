#include "Rythmos_ImplicitRKStepper_decl.hpp"

#ifdef HAVE_RYTHMOS_EXPLICIT_INSTANTIATION

#include "Rythmos_ImplicitRKStepper_def.hpp"
#include "Rythmos_ExplicitInstantiationHelpers.hpp"

namespace Rythmos {

RYTHMOS_MACRO_TEMPLATE_INSTANT_SCALAR_TYPES(RYTHMOS_IMPLICIT_RK_STEPPER_INSTANT) 

} // namespace Rythmos

#endif // HAVE_RYTHMOS_EXPLICIT_INSTANTIATION

