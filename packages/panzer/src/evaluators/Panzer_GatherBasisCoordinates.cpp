#include "Panzer_config.hpp"

#ifdef PANZER_EXPLICIT_TEMPLATE_INSTANTIATION

#include "Panzer_ExplicitTemplateInstantiation.hpp"

#include "Panzer_GatherBasisCoordinates.hpp"
#include "Panzer_GatherBasisCoordinatesT.hpp"

PANZER_INSTANTIATE_TEMPLATE_CLASS_TWO_T(panzer::GatherBasisCoordinates)

#endif
