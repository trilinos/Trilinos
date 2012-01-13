#include "Panzer_config.hpp"

#ifdef PANZER_EXPLICIT_TEMPLATE_INSTANTIATION

#include "Panzer_ExplicitTemplateInstantiation.hpp"

#include "Panzer_GatherIntegrationCoordinates.hpp"
#include "Panzer_GatherIntegrationCoordinatesT.hpp"

PANZER_INSTANTIATE_TEMPLATE_CLASS_TWO_T(panzer::GatherIntegrationCoordinates)

#endif
