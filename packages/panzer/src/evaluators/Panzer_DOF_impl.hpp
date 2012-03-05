#ifndef PANZER_DOF_IMPL_HPP
#define PANZER_DOF_IMPL_HPP

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
  basis_name(p.get< Teuchos::RCP<panzer::BasisIRLayout> >("Basis")->name())
{
  Teuchos::RCP<const PureBasis> basis 
     = p.get< Teuchos::RCP<BasisIRLayout> >("Basis")->getBasis();
  requires_orientation = basis->requiresOrientations();

  // swap between scalar basis value, or vector basis value
  if(basis->isScalarBasis())
     dof_ip = PHX::MDField<ScalarT>(
                p.get<std::string>("Name"), 
     	        p.get< Teuchos::RCP<panzer::IntegrationRule> >("IR")->dl_scalar);
  else if(basis->isVectorBasis())
     dof_ip = PHX::MDField<ScalarT>(
                p.get<std::string>("Name"), 
     	        p.get< Teuchos::RCP<panzer::IntegrationRule> >("IR")->dl_vector);
  else
  { TEUCHOS_ASSERT(false); }

  this->addEvaluatedField(dof_ip);
  this->addDependentField(dof_basis);

  if(requires_orientation) {
     dof_orientation = PHX::MDField<ScalarT,Cell,BASIS>(p.get<std::string>("Name")+" Orientation", 
	                                                basis->functional);
     
     this->addDependentField(dof_orientation);
  }
  
  std::string n = "DOF: " + dof_basis.fieldTag().name();
  this->setName(n);
}

//**********************************************************************
PHX_POST_REGISTRATION_SETUP(DOF,sd,fm)
{
  this->utils.setFieldData(dof_basis,fm);
  this->utils.setFieldData(dof_ip,fm);

  if(requires_orientation)
     this->utils.setFieldData(dof_orientation,fm);

  basis_index = panzer::getBasisIndex(basis_name, (*sd.worksets_)[0]);
}

//**********************************************************************
PHX_EVALUATE_FIELDS(DOF,workset)
{ 
  // Zero out arrays (intrepid does a sum! 1/17/2012)
  for (int i = 0; i < dof_ip.size(); ++i)
    dof_ip[i] = 0.0;

  if(workset.num_cells>0) {
    if(requires_orientation) {
       Intrepid::FieldContainer<double> bases = (workset.bases[basis_index])->basis;

       // assign ScalarT "dof_orientation" to double "orientation"
       Intrepid::FieldContainer<double> orientation(dof_orientation.dimension(0),
                                                    dof_orientation.dimension(1));
       for(int i=0;i<dof_orientation.dimension(0);i++)
          for(int j=0;j<dof_orientation.dimension(1);j++)
             orientation(i,j) = Sacado::ScalarValue<ScalarT>::eval(dof_orientation(i,j));


       // make sure things are orientated correctly
       Intrepid::FunctionSpaceTools::
          applyFieldSigns<ScalarT>(bases,orientation);

       // evaluate at quadrature points
       Intrepid::FunctionSpaceTools::
         evaluate<ScalarT>(dof_ip,dof_basis,bases);
    }
    else // no orientation needed
       Intrepid::FunctionSpaceTools::
         evaluate<ScalarT>(dof_ip,dof_basis,(workset.bases[basis_index])->basis);
  }
}

//**********************************************************************

}

#endif
