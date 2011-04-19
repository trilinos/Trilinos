#ifndef PANZER_EVALUATOR_TRANSIENT_BASISTIMESSCALAR_T_HPP
#define PANZER_EVALUATOR_TRANSIENT_BASISTIMESSCALAR_T_HPP

#include "Intrepid_FunctionSpaceTools.hpp"
#include "Panzer_IntegrationRule.hpp"
#include "Panzer_Basis.hpp"

namespace panzer {

//**********************************************************************
PHX_EVALUATOR_CTOR(Integrator_TransientBasisTimesScalar,p) :
  residual( p.get<std::string>("Residual Name"), 
	    p.get< Teuchos::RCP<panzer::Basis> >("Basis")->functional),
  scalar( p.get<std::string>("Value Name"), 
	  p.get< Teuchos::RCP<panzer::IntegrationRule> >("IR")->dl_scalar),
  basis_name(p.get< Teuchos::RCP<panzer::Basis> >("Basis")->name())
{
  Teuchos::RCP<Teuchos::ParameterList> valid_params = this->getValidParameters();
  p.validateParameters(*valid_params);

  this->addEvaluatedField(residual);
  this->addDependentField(scalar);
    
  multiplier = p.get<double>("Multiplier");


  if (p.isType<Teuchos::RCP<const std::vector<std::string> > >("Field Multipliers")) {
    const std::vector<std::string>& field_multiplier_names = 
      *(p.get<Teuchos::RCP<const std::vector<std::string> > >("Field Multipliers"));

    for (std::vector<std::string>::const_iterator name = 
	   field_multiplier_names.begin(); 
	 name != field_multiplier_names.end(); ++name) {
      PHX::MDField<ScalarT,Cell,IP> tmp_field(*name, p.get< Teuchos::RCP<panzer::IntegrationRule> >("IR")->dl_scalar);
      field_multipliers.push_back(tmp_field);
    }
  }

  for (typename std::vector<PHX::MDField<ScalarT,Cell,IP> >::iterator field = field_multipliers.begin();
       field != field_multipliers.end(); ++field)
    this->addDependentField(*field);

  std::string n = "Integrator_TransientBasisTimesScalar: " + residual.fieldTag().name();
  this->setName(n);
}

//**********************************************************************
PHX_POST_REGISTRATION_SETUP(Integrator_TransientBasisTimesScalar,sd,fm)
{
  this->utils.setFieldData(residual,fm);
  this->utils.setFieldData(scalar,fm);
  
  for (typename std::vector<PHX::MDField<ScalarT,Cell,IP> >::iterator field = field_multipliers.begin();
       field != field_multipliers.end(); ++field)
    this->utils.setFieldData(*field,fm);

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
PHX_EVALUATE_FIELDS(Integrator_TransientBasisTimesScalar,workset)
{ 
  if (workset.evaluate_transient_terms) {
    
    for (int i=0; i < residual.size(); ++i)
      residual[i] = 0.0;
    
    for (std::size_t cell = 0; cell < workset.num_cells; ++cell) {
      for (std::size_t qp = 0; qp < num_qp; ++qp) {
	tmp(cell,qp) = multiplier * scalar(cell,qp);
	for (typename std::vector<PHX::MDField<ScalarT,Cell,IP> >::iterator field = field_multipliers.begin();
	     field != field_multipliers.end(); ++field)
	  tmp(cell,qp) = tmp(cell,qp) * (*field)(cell,qp);  
      }
    }
    Intrepid::FunctionSpaceTools::
      integrate<ScalarT>(residual, tmp, 
			 (workset.bases[basis_index])->weighted_basis, 
			 Intrepid::COMP_BLAS);
  }
}

//**********************************************************************
template<typename EvalT, typename Traits>
Teuchos::RCP<Teuchos::ParameterList> 
Integrator_TransientBasisTimesScalar<EvalT, Traits>::getValidParameters() const
{
  Teuchos::RCP<Teuchos::ParameterList> p = Teuchos::rcp(new Teuchos::ParameterList);
  p->set<std::string>("Residual Name", "?");
  p->set<std::string>("Value Name", "?");
  Teuchos::RCP<panzer::Basis> basis;
  p->set("Basis", basis);
  Teuchos::RCP<panzer::IntegrationRule> ir;
  p->set("IR", ir);
  p->set<double>("Multiplier", 1.0);
  Teuchos::RCP<const std::vector<std::string> > fms;
  p->set("Field Multipliers", fms);
  return p;
}

//**********************************************************************

}

#endif

