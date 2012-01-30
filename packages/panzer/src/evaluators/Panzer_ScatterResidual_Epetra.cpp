#include "Panzer_config.hpp"

#ifdef HAVE_PANZER_EXPLICIT_INSTANTIATION

#include "Panzer_ExplicitTemplateInstantiation.hpp"

#include "Panzer_ScatterResidual_Epetra_decl.hpp"
#include "Panzer_ScatterResidual_Epetra_impl.hpp"

#include "Panzer_ScatterResidual_EpetraSG_decl.hpp"
#include "Panzer_ScatterResidual_EpetraSG_impl.hpp"

PANZER_INSTANTIATE_TEMPLATE_CLASS_FOUR_T(panzer::ScatterResidual_Epetra,int,int)
PANZER_INSTANTIATE_TEMPLATE_CLASS_FOUR_T(panzer::ScatterResidual_Epetra,short,int)
PANZER_INSTANTIATE_TEMPLATE_CLASS_FOUR_T(panzer::ScatterResidual_Epetra,char,int)

#endif
