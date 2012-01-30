#include "Panzer_config.hpp"

#ifdef HAVE_PANZER_EXPLICIT_INSTANTIATION

#include "Panzer_ExplicitTemplateInstantiation.hpp"

#include "Panzer_GatherOrientation_decl.hpp"
#include "Panzer_GatherOrientation_impl.hpp"

PANZER_INSTANTIATE_TEMPLATE_CLASS_FOUR_T(panzer::GatherOrientation,int,int)
PANZER_INSTANTIATE_TEMPLATE_CLASS_FOUR_T(panzer::GatherOrientation,short,int)
PANZER_INSTANTIATE_TEMPLATE_CLASS_FOUR_T(panzer::GatherOrientation,char,int)

#endif
