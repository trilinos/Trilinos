
#include <algorithm>
#include "Panzer_IntegrationRule.hpp"
#include "Panzer_BasisIRLayout.hpp"
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

  basis_index = panzer::getBasisIndex(basis_name, (*sd.worksets_)[0]);
}

//**********************************************************************
PHX_EVALUATE_FIELDS(DOF,workset)
{ 
/*
  if (typeid(EvalT) == typeid(panzer::Traits::Jacobian)) {
     std::cout << dof_basis.fieldTag().name() << " BASIS =\n";
     for (int i=0; i < dof_basis.dimension(0); ++i) {
        std::cout << "   Cell " << workset.cell_local_ids[i] << std::endl;
        for (int j=0; j < dof_basis.dimension(1); ++j)
            std::cout << "      " << dof_basis(i,j) << std::endl;
     }

     for (int i=0; i < (workset.bases[basis_index])->basis.dimension(0);i++) {
        std::cout << "   Cell " << workset.cell_local_ids[i] << std::endl;
        for (int j=0; j < (workset.bases[basis_index])->basis.dimension(1);j++) {
           std::cout << "      basis " << j << " = ";
           for (int ip=0; ip<(workset.bases[basis_index])->basis.dimension(2);ip++)
              std::cout << (workset.bases[basis_index])->basis(i,j,ip) << " ";
           std::cout << std::endl;
        }
     }
  }
*/

  // Zero out arrays (intrepid does a sum! 1/17/2012)
  for (int i = 0; i < dof_ip.size(); ++i)
    dof_ip[i] = 0.0;

/*
  if (typeid(EvalT) == typeid(panzer::Traits::Jacobian)) {
     std::cout << "BEFORE " << dof_ip.fieldTag().name() << " IP =\n";
     for (int i=0; i < dof_ip.dimension(0); ++i) {
        std::cout << "   Cell " << workset.cell_local_ids[i] << std::endl;
        for (int j=0; j < dof_ip.dimension(1); ++j)
            std::cout << "      " << dof_ip(i,j) << std::endl;
     }
  }
*/

  if(workset.num_cells>0)
    Intrepid::FunctionSpaceTools::
      evaluate<ScalarT>(dof_ip,dof_basis,(workset.bases[basis_index])->basis);

/*
  if (typeid(EvalT) == typeid(panzer::Traits::Jacobian)) {
     std::cout << "AFTER " << dof_ip.fieldTag().name() << " IP =\n";
     for (int i=0; i < dof_ip.dimension(0); ++i) {
        std::cout << "   Cell " << workset.cell_local_ids[i] << std::endl;
        for (int j=0; j < dof_ip.dimension(1); ++j)
            std::cout << "      " << dof_ip(i,j) << std::endl;
     }
  }
*/
}

//**********************************************************************

}
