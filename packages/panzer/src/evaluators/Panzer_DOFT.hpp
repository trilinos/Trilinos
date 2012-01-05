
#include <algorithm>
#include "Panzer_IntegrationRule.hpp"
#include "Panzer_Basis.hpp"
#include "Panzer_Workset_Utilities.hpp"
#include "Intrepid_FunctionSpaceTools.hpp"

namespace panzer {

//**********************************************************************
PHX_EVALUATOR_CTOR(DOF,p) :
  dof_basis( p.get<std::string>("Name"), 
	     p.get< Teuchos::RCP<panzer::BasisIRLayout> >("Basis")->functional),
  dof_ip( p.get<std::string>("Name"), 
	  p.get< Teuchos::RCP<panzer::IntegrationRule> >("IR")->dl_scalar),
  basis_name(p.get< Teuchos::RCP<panzer::BasisIRLayout> >("Basis")->name())
{
  this->addEvaluatedField(dof_ip);
  this->addDependentField(dof_basis);
  
  std::string n = "DOF: " + dof_basis.fieldTag().name();
  this->setName(n);
}

//**********************************************************************
PHX_POST_REGISTRATION_SETUP(DOF,sd,fm)
{
  this->utils.setFieldData(dof_basis,fm);
  this->utils.setFieldData(dof_ip,fm);

  cell_data_size = dof_basis.size() / 
    dof_basis.fieldTag().dataLayout().dimension(0);

  basis_index = panzer::getBasisIndex(basis_name, (*sd.worksets_)[0]);
}

//**********************************************************************
PHX_EVALUATE_FIELDS(DOF,workset)
{ 
  // Zero out arrays (probably not needed anymore intrepid does replace)
  std::size_t length = workset.num_cells * cell_data_size;
  for (std::size_t i = 0; i < length; ++i)
    dof_ip[i] = 0.0;

  if(workset.num_cells>0)
    Intrepid::FunctionSpaceTools::
      evaluate<ScalarT>(dof_ip,dof_basis,(workset.bases[basis_index])->basis);
}

//**********************************************************************

}
