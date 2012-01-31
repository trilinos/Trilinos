#include "Panzer_config.hpp"

#ifdef HAVE_PANZER_EXPLICIT_INSTANTIATION

#include "Panzer_ExplicitTemplateInstantiation.hpp"

#include "Panzer_STK_ScatterCellAvgQuantity_decl.hpp"
#include "Panzer_STK_ScatterCellAvgQuantity_impl.hpp"

PANZER_INSTANTIATE_TEMPLATE_CLASS_TWO_T(panzer_stk::ScatterCellAvgQuantity)

#endif
