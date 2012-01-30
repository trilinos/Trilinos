#include "Panzer_config.hpp"

#ifdef HAVE_PANZER_EXPLICIT_INSTANTIATION

#include "Panzer_ExplicitTemplateInstantiation.hpp"

#include "Panzer_ScatterDirichletResidual_Epetra_decl.hpp"
#include "Panzer_ScatterDirichletResidual_Epetra_impl.hpp"

#include "Panzer_ScatterDirichletResidual_EpetraSG_decl.hpp"
#include "Panzer_ScatterDirichletResidual_EpetraSG_impl.hpp"

PANZER_INSTANTIATE_TEMPLATE_CLASS_FOUR_T(panzer::ScatterDirichletResidual_Epetra,int,int)
PANZER_INSTANTIATE_TEMPLATE_CLASS_FOUR_T(panzer::ScatterDirichletResidual_Epetra,short,int)
PANZER_INSTANTIATE_TEMPLATE_CLASS_FOUR_T(panzer::ScatterDirichletResidual_Epetra,char,int)

#endif
