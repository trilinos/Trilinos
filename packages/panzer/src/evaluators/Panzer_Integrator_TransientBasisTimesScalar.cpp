#include "Panzer_config.hpp"

#ifdef HAVE_PANZER_EXPLICIT_INSTANTIATION

#include "Panzer_ExplicitTemplateInstantiation.hpp"

#include "Panzer_Integrator_TransientBasisTimesScalar_decl.hpp"
#include "Panzer_Integrator_TransientBasisTimesScalar_impl.hpp"

PANZER_INSTANTIATE_TEMPLATE_CLASS_TWO_T(panzer::Integrator_TransientBasisTimesScalar)

#endif
