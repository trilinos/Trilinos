#include "Panzer_config.hpp"

#ifdef HAVE_PANZER_EXPLICIT_INSTANTIATION

#include "Panzer_ExplicitTemplateInstantiation.hpp"

#include "Panzer_GatherSolution_Epetra_decl.hpp"
#include "Panzer_GatherSolution_Epetra_impl.hpp"

#include "Panzer_GatherSolution_EpetraSG_decl.hpp"
#include "Panzer_GatherSolution_EpetraSG_impl.hpp"

PANZER_INSTANTIATE_TEMPLATE_CLASS_FOUR_T(panzer::GatherSolution_Epetra,int,int)
PANZER_INSTANTIATE_TEMPLATE_CLASS_FOUR_T(panzer::GatherSolution_Epetra,short,int)
PANZER_INSTANTIATE_TEMPLATE_CLASS_FOUR_T(panzer::GatherSolution_Epetra,char,int)

#endif
