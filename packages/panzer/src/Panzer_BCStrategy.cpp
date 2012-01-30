#include "Panzer_config.hpp"

#ifdef HAVE_PANZER_EXPLICIT_INSTANTIATION

#include "Panzer_ExplicitTemplateInstantiation.hpp"

#include "Panzer_BCStrategy_decl.hpp"
#include "Panzer_BCStrategy_impl.hpp"

PANZER_INSTANTIATE_TEMPLATE_CLASS_ONE_T(panzer::BCStrategy)

#endif
