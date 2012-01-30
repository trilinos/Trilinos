#ifndef PANZER_EVALUATOR_GRADBASISDOTVECTOR_IMPL_HPP
#define PANZER_EVALUATOR_GRADBASISDOTVECTOR_IMPL_HPP

#include "Intrepid_FunctionSpaceTools.hpp"
#include "Panzer_IntegrationRule.hpp"
#include "Panzer_BasisIRLayout.hpp"
#include "Panzer_Workset_Utilities.hpp"

namespace panzer {

//**********************************************************************
PHX_EVALUATOR_CTOR(Integrator_GradBasisDotVector,p) :
  residual( p.get<std::string>("Residual Name"), 
	    p.get< Teuchos::RCP<panzer::BasisIRLayout> >("Basis")->functional),
  flux( p.get<std::string>("Flux Name"), 
	p.get< Teuchos::RCP<panzer::IntegrationRule> >("IR")->dl_vector ),
  basis_name(p.get< Teuchos::RCP<panzer::BasisIRLayout> >("Basis")->name())
{
  this->addEvaluatedField(residual);
  this->addDependentField(flux);
  
  multiplier = p.get<double>("Multiplier");

  if (p.isType<Teuchos::RCP<const std::vector<std::string> > >("Field Multipliers")) 
  {
    const std::vector<std::string>& field_multiplier_names = 
      *(p.get<Teuchos::RCP<const std::vector<std::string> > >("Field Multipliers"));

    for (std::vector<std::string>::const_iterator name = field_multiplier_names.begin(); 
      name != field_multiplier_names.end(); ++name) 
    {
      PHX::MDField<ScalarT,Cell,IP> tmp_field(*name, p.get< Teuchos::RCP<panzer::IntegrationRule> >("IR")->dl_scalar);
      field_multipliers.push_back(tmp_field);
    }
  }

  for (typename std::vector<PHX::MDField<ScalarT,Cell,IP> >::iterator field = field_multipliers.begin();
       field != field_multipliers.end(); ++field)
    this->addDependentField(*field);

  std::string n = 
    "Integrator_GradBasisDotVector: " + residual.fieldTag().name();

  this->setName(n);
}

//**********************************************************************
PHX_POST_REGISTRATION_SETUP(Integrator_GradBasisDotVector,sd,fm)
{
  this->utils.setFieldData(residual,fm);
  this->utils.setFieldData(flux,fm);

  for (typename std::vector<PHX::MDField<ScalarT,Cell,IP> >::iterator field = field_multipliers.begin();
       field != field_multipliers.end(); ++field)
    this->utils.setFieldData(*field,fm);

  num_nodes = residual.dimension(1);
  num_qp = flux.dimension(1);
  num_dim = flux.dimension(2);

  basis_index = panzer::getBasisIndex(basis_name, (*sd.worksets_)[0]);

  tmp = Intrepid::FieldContainer<ScalarT>(flux.dimension(0), num_qp, num_dim); 
}

//**********************************************************************
PHX_EVALUATE_FIELDS(Integrator_GradBasisDotVector,workset)
{ 
  for (int i=0; i < residual.size(); ++i)
    residual[i] = 0.0;
  
  for (std::size_t cell = 0; cell < workset.num_cells; ++cell)
  {
    for (std::size_t qp = 0; qp < num_qp; ++qp)
    {
      ScalarT tmpVar = 1.0;
      for (typename std::vector<PHX::MDField<ScalarT,Cell,IP> >::iterator field = field_multipliers.begin();
           field != field_multipliers.end(); ++field)
        tmpVar = tmpVar * (*field)(cell,qp);  

      for (std::size_t dim = 0; dim < num_dim; ++dim)
        tmp(cell,qp,dim) = multiplier * tmpVar * flux(cell,qp,dim);
    }
  }
  
  if(workset.num_cells>0)
     Intrepid::FunctionSpaceTools::
       integrate<ScalarT>(residual, tmp, 
   		       (workset.bases[basis_index])->weighted_grad_basis, 
		       Intrepid::COMP_BLAS);
}

//**********************************************************************

}

#endif
