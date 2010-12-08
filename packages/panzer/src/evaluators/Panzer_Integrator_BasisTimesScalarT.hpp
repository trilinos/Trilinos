#ifndef PANZER_EVALUATOR_BASISTIMESSCALAR_T_HPP
#define PANZER_EVALUATOR_BASISTIMESSCALAR_T_HPP

#include "Intrepid_FunctionSpaceTools.hpp"
#include "Panzer_IntegrationRule.hpp"
#include "Panzer_Basis.hpp"

namespace panzer {

//**********************************************************************
PHX_EVALUATOR_CTOR(Integrator_BasisTimesScalar,p) :
  residual( p.get<std::string>("Residual Name"), 
	    p.get< Teuchos::RCP<panzer::Basis> >("Basis")->functional),
  scalar( p.get<std::string>("Value Name"), 
	  p.get< Teuchos::RCP<panzer::IntegrationRule> >("IR")->dl_scalar),
  basis_name(p.get< Teuchos::RCP<panzer::Basis> >("Basis")->name())
{
  this->addEvaluatedField(residual);
  this->addDependentField(scalar);
    
  multiplier = p.get<double>("Multiplier");

  std::string n = "Integrator_BasisTimesScalar: " + residual.fieldTag().name();
  this->setName(n);
}

//**********************************************************************
PHX_POST_REGISTRATION_SETUP(Integrator_BasisTimesScalar,sd,fm)
{
  this->utils.setFieldData(residual,fm);
  this->utils.setFieldData(scalar,fm);

  num_nodes = residual.dimension(1);
  num_qp = scalar.dimension(1);

  basis_index = 
    std::distance((*sd.worksets_)[0].basis_names->begin(),
		  std::find((*sd.worksets_)[0].basis_names->begin(),
			    (*sd.worksets_)[0].basis_names->end(),
			    basis_name));

  tmp = Intrepid::FieldContainer<ScalarT>(scalar.dimension(0), num_qp); 
}

//**********************************************************************
PHX_EVALUATE_FIELDS(Integrator_BasisTimesScalar,workset)
{ 
  for (int i=0; i < residual.size(); ++i)
    residual[i] = 0.0;

  for (std::size_t cell = 0; cell < workset.num_cells; ++cell)
    for (std::size_t qp = 0; qp < num_qp; ++qp)
      tmp(cell,qp) = multiplier * scalar(cell,qp);
  
  Intrepid::FunctionSpaceTools::
    integrate<ScalarT>(residual, tmp, 
		       (workset.bases[basis_index])->weighted_basis, 
		       Intrepid::COMP_BLAS);
}

//**********************************************************************

}

#endif

