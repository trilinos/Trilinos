#ifndef PANZER_EQUATIONSET_DEFAULT_IMPL_T_HPP
#define PANZER_EQUATIONSET_DEFAULT_IMPL_T_HPP

#include "Panzer_GatherOrientation.hpp"
#include "Panzer_DOF.hpp"
#include "Panzer_DOFGradient.hpp"

#include "Phalanx_MDField.hpp"
#include "Phalanx_DataLayout.hpp"
#include "Phalanx_DataLayout_MDALayout.hpp"
#include "Teuchos_RCP.hpp"

// ***********************************************************************
template <typename EvalT>
panzer::EquationSet_DefaultImpl<EvalT>::
EquationSet_DefaultImpl(const panzer::InputEquationSet& ies,
			const panzer::CellData& cell_data,
			const bool build_transient_support) :
  m_input_eq_set(ies),
  m_cell_data(cell_data),
  m_build_transient_support(build_transient_support),
  m_eqset_prefix("")
{ 
  m_dof_names = Teuchos::rcp(new std::vector<std::string>);
  m_dof_gradient_names = Teuchos::rcp(new std::vector<std::string>);
  m_dof_time_derivative_names = Teuchos::rcp(new std::vector<std::string>);
  m_residual_names = Teuchos::rcp(new std::vector<std::string>);
  m_eval_plist = Teuchos::rcp(new Teuchos::ParameterList);
}

// ***********************************************************************
template <typename EvalT>
void panzer::EquationSet_DefaultImpl<EvalT>::
setupDOFs(int equation_dimension)
{
  this->m_int_rule = 
    Teuchos::rcp(new panzer::IntegrationRule(m_input_eq_set.integration_order,
					     m_cell_data));
  
  this->m_basis = Teuchos::rcp(new panzer::Basis(m_input_eq_set.basis,
						 *(this->m_int_rule)));
  
  this->m_provided_dofs.clear();
  int index = 0;
  for (std::vector<std::string>::const_iterator i = this->m_dof_names->begin();
       i != this->m_dof_names->end(); ++i, ++index)
    this->m_provided_dofs.push_back(std::make_pair((*(this->m_dof_names))[index], this->m_basis));

  this->m_eval_plist->set("IR", this->m_int_rule);
  this->m_eval_plist->set("Basis", this->m_basis);
  this->m_eval_plist->set("Equation Dimension", equation_dimension);
  this->m_eval_plist->set("DOF Names", this->m_dof_names);  
}

// ***********************************************************************
template <typename EvalT>
void panzer::EquationSet_DefaultImpl<EvalT>::
buildAndRegisterGatherScatterEvaluators(PHX::FieldManager<panzer::Traits>& fm,
					const std::vector<std::pair<std::string,Teuchos::RCP<panzer::Basis> > > & dofs,
                                        const LinearObjFactory<panzer::Traits> & lof,
					const Teuchos::ParameterList& user_data) const
{
  using Teuchos::ParameterList;
  using Teuchos::RCP;
  using Teuchos::rcp;

  // this turns off the scatter contribution, and does
  // only the gather
  bool ignoreScatter = false;
  if(user_data.isParameter("Ignore Scatter")) 
     ignoreScatter = user_data.get<bool>("Ignore Scatter");

  // ********************
  // DOFs (unknowns)
  // ********************

  // Gather, includes construction of orientation gathers
  {
    ParameterList p("Gather");
    p.set("Basis", m_basis);
    p.set("DOF Names", m_dof_names);
    p.set("Indexer Names", m_dof_names);
    
    RCP< PHX::Evaluator<panzer::Traits> > op = lof.buildGather<EvalT>(p);
    
    fm.template registerEvaluator<EvalT>(op);
  }

  if(!ignoreScatter) {
     // Scatter
     RCP<std::map<std::string,std::string> > names_map =
       rcp(new std::map<std::string,std::string>);
   
     TEUCHOS_ASSERT(m_dof_names->size() == m_residual_names->size());
   
     for (std::size_t i = 0; i < m_dof_names->size(); ++i)
       names_map->insert(std::make_pair((*m_residual_names)[i],(*m_dof_names)[i]));
   
     {
       ParameterList p("Scatter");
       p.set("Scatter Name", this->m_scatter_name);
       p.set("Basis", this->m_basis);
       p.set("Dependent Names", this->m_residual_names);
       p.set("Dependent Map", names_map);
   
       RCP< PHX::Evaluator<panzer::Traits> > op = lof.buildScatter<EvalT>(p);
         // rcp(new panzer::ScatterResidual_Epetra<EvalT,panzer::Traits>(p));
       
       fm.template registerEvaluator<EvalT>(op);
     }
     
     // Require variables
     {
       PHX::Tag<typename EvalT::ScalarT> tag(this->m_scatter_name, 
   					  Teuchos::rcp(new PHX::MDALayout<Dummy>(0)));
       fm.template requireField<EvalT>(tag);
     }
  }

  // DOFs: Scalar value @ basis --> Scalar value @ IP 
  for (std::vector<std::string>::const_iterator dof_name = this->m_dof_names->begin();
       dof_name != this->m_dof_names->end(); ++dof_name) {
    ParameterList p;
    p.set("Name", *dof_name);
    p.set("Basis", this->m_basis); 
    p.set("IR", this->m_int_rule);
    
    RCP< PHX::Evaluator<panzer::Traits> > op = 
      rcp(new panzer::DOF<EvalT,panzer::Traits>(p));
    
    fm.template registerEvaluator<EvalT>(op);
  }

  TEUCHOS_ASSERT(this->m_dof_names->size() == this->m_dof_gradient_names->size());

  // Gradients of DOFs: Scalar value @ basis --> Vector value @ IP
  {
    std::vector<std::string>::const_iterator dof_name = this->m_dof_names->begin();
    std::vector<std::string>::const_iterator dof_grad_name = this->m_dof_gradient_names->begin();
    for (;dof_grad_name != this->m_dof_gradient_names->end(); ++dof_grad_name,++dof_name) {
      ParameterList p;
      p.set("Name", *dof_name);
      p.set("Gradient Name", *dof_grad_name);
      p.set("Basis", this->m_basis); 
      p.set("IR", this->m_int_rule);
      
      RCP< PHX::Evaluator<panzer::Traits> > op = 
	rcp(new panzer::DOFGradient<EvalT,panzer::Traits>(p));

      fm.template registerEvaluator<EvalT>(op);
    }
  }

  // **************************
  // Time derivative terms
  // **************************

  // Gather of time derivative terms
  {
    ParameterList p("Gather");
    p.set("Basis", m_basis);
    p.set("DOF Names", m_dof_time_derivative_names);
    p.set("Indexer Names", m_dof_names);
    p.set("Use Time Derivative Solution Vector", true);
    
    RCP< PHX::Evaluator<panzer::Traits> > op = lof.buildGather<EvalT>(p);
    //   rcp(new panzer::GatherSolution_Epetra<EvalT,panzer::Traits>(p));
    
    fm.template registerEvaluator<EvalT>(op);
  }
  
  // Time derivative of DOFs: Scalar value @ basis --> Scalar value @ IP 
  for (std::vector<std::string>::const_iterator dof_name = this->m_dof_time_derivative_names->begin();
       dof_name != this->m_dof_time_derivative_names->end(); ++dof_name) {
    ParameterList p;
    p.set("Name", *dof_name);
    p.set("Basis", this->m_basis); 
    p.set("IR", this->m_int_rule);
    
    RCP< PHX::Evaluator<panzer::Traits> > op = 
      rcp(new panzer::DOF<EvalT,panzer::Traits>(p));
    
    fm.template registerEvaluator<EvalT>(op);
  }

}

// ***********************************************************************
template <typename EvalT>
void panzer::EquationSet_DefaultImpl<EvalT>::
buildAndRegisterClosureModelEvaluators(PHX::FieldManager<panzer::Traits>& fm,
				       const std::vector<std::pair<std::string,Teuchos::RCP<panzer::Basis> > > & dofs,
				       const panzer::ClosureModelFactory_TemplateManager<panzer::Traits>& factory,
				       const Teuchos::ParameterList& models,
				       const Teuchos::ParameterList& user_data) const
{
/*
  Teuchos::RCP< std::vector< Teuchos::RCP<PHX::Evaluator<panzer::Traits> > > > evaluators = 
    factory.getAsObject<EvalT>()->buildClosureModels(this->m_input_eq_set.model_id, this->m_input_eq_set, models,
						     *(this->m_eval_plist), user_data, fm);
  
  for (std::vector< Teuchos::RCP<PHX::Evaluator<panzer::Traits> > >::size_type i=0; i < evaluators->size(); ++i)
    fm.template registerEvaluator<EvalT>((*evaluators)[i]);
*/
  buildAndRegisterClosureModelEvaluators(fm,dofs,factory,this->m_input_eq_set.model_id,models,user_data);
}

// ***********************************************************************

template <typename EvalT>
void panzer::EquationSet_DefaultImpl<EvalT>::
buildAndRegisterClosureModelEvaluators(PHX::FieldManager<panzer::Traits>& fm,
				   const std::vector<std::pair<std::string,Teuchos::RCP<panzer::Basis> > > & dofs,
				   const panzer::ClosureModelFactory_TemplateManager<panzer::Traits>& factory,
				   const std::string& model_name,
				   const Teuchos::ParameterList& models,
				   const Teuchos::ParameterList& user_data) const
{
  Teuchos::RCP< std::vector< Teuchos::RCP<PHX::Evaluator<panzer::Traits> > > > evaluators = 
    factory.getAsObject<EvalT>()->buildClosureModels(model_name, this->m_input_eq_set, models, 
                                                     *(this->m_eval_plist), user_data, fm);
    
  for (std::vector< Teuchos::RCP<PHX::Evaluator<panzer::Traits> > >::size_type i=0; i < evaluators->size(); ++i)
    fm.template registerEvaluator<EvalT>((*evaluators)[i]);
}

// ***********************************************************************
template <typename EvalT>
void panzer::EquationSet_DefaultImpl<EvalT>::
buildAndRegisterInitialConditionEvaluators(PHX::FieldManager<panzer::Traits>& fm,
					   const std::vector<std::pair<std::string,Teuchos::RCP<panzer::Basis> > > & dofs,
					   const panzer::ClosureModelFactory_TemplateManager<panzer::Traits>& factory,
					   const std::string& model_name,
					   const Teuchos::ParameterList& models,
					   const LinearObjFactory<panzer::Traits> & lof,
					   const Teuchos::ParameterList& user_data) const
{

  {
    Teuchos::ParameterList p("Scatter");
    p.set("Scatter Name", this->m_scatter_name);
    p.set("Basis", this->m_basis);
    p.set("Dependent Names", this->m_dof_names);

    Teuchos::RCP< PHX::Evaluator<panzer::Traits> > op = lof.buildScatterInitialCondition<EvalT>(p);
      // rcp(new panzer::ScatterResidual_Epetra<EvalT,panzer::Traits>(p));
    
    fm.template registerEvaluator<EvalT>(op);
  }

  // Require variables
  {
    PHX::Tag<typename EvalT::ScalarT> tag(this->m_scatter_name, 
					  Teuchos::rcp(new PHX::MDALayout<Dummy>(0)));
    fm.template requireField<EvalT>(tag);
  }  

  // Add in closure models
  {
    Teuchos::RCP< std::vector< Teuchos::RCP<PHX::Evaluator<panzer::Traits> > > > evaluators = 
      factory.getAsObject<EvalT>()->buildClosureModels(model_name, this->m_input_eq_set, models, *(this->m_eval_plist), user_data, fm);
    
    for (std::vector< Teuchos::RCP<PHX::Evaluator<panzer::Traits> > >::size_type i=0; i < evaluators->size(); ++i)
      fm.template registerEvaluator<EvalT>((*evaluators)[i]);
  }

}

// ***********************************************************************
template <typename EvalT>
const Teuchos::RCP<Teuchos::ParameterList>
panzer::EquationSet_DefaultImpl<EvalT>::getEvaluatorParameterList() const
{
  return m_eval_plist;
}

// ***********************************************************************
template <typename EvalT>
const std::vector<std::string>&
panzer::EquationSet_DefaultImpl<EvalT>::getDOFNames() const
{
  return *m_dof_names;
}
// ***********************************************************************
template <typename EvalT>
const std::vector<std::pair<std::string,Teuchos::RCP<panzer::Basis> > >&
panzer::EquationSet_DefaultImpl<EvalT>::getProvidedDOFs() const
{
  return m_provided_dofs;
}

// ***********************************************************************

#endif
