#include "Panzer_config.hpp"

#ifdef PANZER_EXPLICIT_TEMPLATE_INSTANTIATION

#include "Panzer_ExplicitTemplateInstantiation.hpp"

#include "Panzer_STK_ScatterCellAvgQuantity.hpp"
#include "Panzer_STK_ScatterCellAvgQuantityT.hpp"

PANZER_INSTANTIATE_TEMPLATE_CLASS_TWO_T(panzer_stk::ScatterCellAvgQuantity)

#endif
