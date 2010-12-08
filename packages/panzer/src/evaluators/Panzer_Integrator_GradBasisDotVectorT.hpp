#ifndef PANZER_EVALUATOR_GRADBASISDOTVECTOR_T_HPP
#define PANZER_EVALUATOR_GRADBASISDOTVECTOR_T_HPP

#include "Intrepid_FunctionSpaceTools.hpp"
#include "Panzer_IntegrationRule.hpp"
#include "Panzer_Basis.hpp"

namespace panzer {

//**********************************************************************
PHX_EVALUATOR_CTOR(Integrator_GradBasisDotVector,p) :
  residual( p.get<std::string>("Residual Name"), 
	    p.get< Teuchos::RCP<panzer::Basis> >("Basis")->functional),
  flux( p.get<std::string>("Flux Name"), 
	p.get< Teuchos::RCP<panzer::IntegrationRule> >("IR")->dl_vector ),
  basis_name(p.get< Teuchos::RCP<panzer::Basis> >("Basis")->name())
{
  this->addEvaluatedField(residual);
  this->addDependentField(flux);
  
  multiplier = p.get<double>("Multiplier");

  std::string n = 
    "Integrator_GradBasisDotVector: " + residual.fieldTag().name();

  this->setName(n);
}

//**********************************************************************
PHX_POST_REGISTRATION_SETUP(Integrator_GradBasisDotVector,sd,fm)
{
  this->utils.setFieldData(residual,fm);
  this->utils.setFieldData(flux,fm);

  num_nodes = residual.dimension(1);
  num_qp = flux.dimension(1);
  num_dim = flux.dimension(2);

  basis_index = 
    std::distance((*sd.worksets_)[0].basis_names->begin(),
		  std::find((*sd.worksets_)[0].basis_names->begin(),
			    (*sd.worksets_)[0].basis_names->end(),
			    basis_name));

  tmp = Intrepid::FieldContainer<ScalarT>(flux.dimension(0), num_qp, num_dim); 
}

//**********************************************************************
PHX_EVALUATE_FIELDS(Integrator_GradBasisDotVector,workset)
{ 
  for (int i=0; i < residual.size(); ++i)
    residual[i] = 0.0;
  
  for (std::size_t cell = 0; cell < workset.num_cells; ++cell)
    for (std::size_t qp = 0; qp < num_qp; ++qp)
      for (std::size_t dim = 0; dim < num_dim; ++dim)
	tmp(cell,qp,dim) = multiplier * flux(cell,qp,dim);
  
  Intrepid::FunctionSpaceTools::
    integrate<ScalarT>(residual, tmp, 
		       (workset.bases[basis_index])->weighted_grad_basis, 
		       Intrepid::COMP_BLAS);
}

//**********************************************************************

}

#endif
