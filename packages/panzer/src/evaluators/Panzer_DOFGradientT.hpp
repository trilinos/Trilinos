#include "Panzer_IntegrationRule.hpp"
#include "Panzer_Basis.hpp"
#include "Panzer_Workset_Utilities.hpp"
#include "Intrepid_FunctionSpaceTools.hpp"

namespace panzer {

//**********************************************************************
PHX_EVALUATOR_CTOR(DOFGradient,p) :
  dof_value( p.get<std::string>("Name"), 
	     p.get< Teuchos::RCP<panzer::BasisIRLayout> >("Basis")->functional),
  dof_gradient( p.get<std::string>("Gradient Name"), 
		p.get< Teuchos::RCP<panzer::IntegrationRule> >("IR")->dl_vector ),
  basis_name(p.get< Teuchos::RCP<panzer::BasisIRLayout> >("Basis")->name())
{
  this->addEvaluatedField(dof_gradient);
  this->addDependentField(dof_value);
  
  std::string n = "DOFGradient: " + dof_gradient.fieldTag().name();
  this->setName(n);
}

//**********************************************************************
PHX_POST_REGISTRATION_SETUP(DOFGradient,sd,fm)
{
  this->utils.setFieldData(dof_value,fm);
  this->utils.setFieldData(dof_gradient,fm);

  cell_data_size = dof_gradient.size() / 
    dof_gradient.fieldTag().dataLayout().dimension(0);

  basis_index = panzer::getBasisIndex(basis_name, (*sd.worksets_)[0]);
}

//**********************************************************************
PHX_EVALUATE_FIELDS(DOFGradient,workset)
{ 
  // Zero out arrays (probably don't need this anymore)
  std::size_t length = workset.num_cells * cell_data_size;
  for (std::size_t i = 0; i < length; ++i)
    dof_gradient[i] = 0.0;

  if(workset.num_cells>0)
     Intrepid::FunctionSpaceTools::evaluate<ScalarT>(dof_gradient,dof_value,(workset.bases[basis_index])->grad_basis);
}

//**********************************************************************

}
