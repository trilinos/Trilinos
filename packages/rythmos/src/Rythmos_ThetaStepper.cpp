#include "Rythmos_ThetaStepper_decl.hpp"

#ifdef HAVE_RYTHMOS_EXPLICIT_INSTANTIATION

#ifdef HAVE_RYTHMOS_EXPERIMENTAL

#include "Rythmos_ThetaStepper_def.hpp"
#include "Rythmos_ExplicitInstantiationHelpers.hpp"

namespace Rythmos {

RYTHMOS_MACRO_TEMPLATE_INSTANT_SCALAR_TYPES(RYTHMOS_THETA_STEPPER_INSTANT) 

} // namespace Rythmos

#endif // HAVE_RYTHMOS_EXPERIMENTAL

#endif // HAVE_RYTHMOS_EXPLICIT_INSTANTIATION

