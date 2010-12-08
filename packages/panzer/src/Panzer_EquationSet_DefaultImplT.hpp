#ifndef PANZER_EQUATIONSET_DEFAULT_IMPL_T_HPP
#define PANZER_EQUATIONSET_DEFAULT_IMPL_T_HPP

#include "Panzer_GatherSolution_Epetra.hpp"
#include "Panzer_ScatterResidual_Epetra.hpp"
#include "Panzer_DOF.hpp"
#include "Panzer_DOFGradient.hpp"

// ***********************************************************************
template <typename EvalT>
panzer::EquationSet_DefaultImpl<EvalT>::
EquationSet_DefaultImpl(const panzer::InputEquationSet& ies,
			const panzer::CellData& cell_data) :
  m_input_eq_set(ies),
  m_cell_data(cell_data),
  m_eqset_prefix("")
{ }

// ***********************************************************************
template <typename EvalT>
void panzer::EquationSet_DefaultImpl<EvalT>::
buildAndRegisterGatherScatterEvaluators(PHX::FieldManager<panzer::Traits>& fm,
					const std::vector<std::pair<std::string,Teuchos::RCP<panzer::Basis> > > & dofs) const
{
  using Teuchos::ParameterList;
  using Teuchos::RCP;
  using Teuchos::rcp;

  // ********************
  // DOFs (unknowns)
  // ********************

  // Gather
  {
    ParameterList p("Gather");
    p.set("Basis", m_basis);
    p.set("DOF Names", m_dof_names);
    
    RCP< PHX::Evaluator<panzer::Traits> > op = 
      rcp(new panzer::GatherSolution_Epetra<EvalT,panzer::Traits>(p));
    
    fm.template registerEvaluator<EvalT>(op);
  }
  
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

    RCP< PHX::Evaluator<panzer::Traits> > op = 
      rcp(new panzer::ScatterResidual_Epetra<EvalT,panzer::Traits>(p));
    
    fm.template registerEvaluator<EvalT>(op);
  }
  
  // Require variables
  {
    std::string reqField = m_eqset_prefix+"Scatter_EnergyRawEvaluators";

    PHX::Tag<typename EvalT::ScalarT> tag(this->m_scatter_name, 
					  Teuchos::rcp(new PHX::MDALayout<Dummy>(0)));
    fm.template requireField<EvalT>(tag);
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

}

// ***********************************************************************
template <typename EvalT>
void panzer::EquationSet_DefaultImpl<EvalT>::
buildAndRegisterModelEvaluators(PHX::FieldManager<panzer::Traits>& fm,
				const std::vector<std::pair<std::string,Teuchos::RCP<panzer::Basis> > > & dofs,
				const std::map<std::string,Teuchos::RCP<panzer::ModelFactory_TemplateManager<panzer::Traits> > >& factories,
				const std::vector<Teuchos::ParameterList>& models) const
{
  Teuchos::RCP< std::vector< Teuchos::RCP<PHX::Evaluator<panzer::Traits> > > > evaluators;
  
  Teuchos::RCP<Teuchos::ParameterList> default_params = this->getEvaluatorParameterList();

  evaluators = factories.find(m_input_eq_set.model_factory)->second->getAsObject<EvalT>()->buildModels(m_input_eq_set, models, *default_params);
  
  // Loop over evaluators and register them with field manager
  for (std::size_t i=0; i < evaluators->size(); ++i)
    fm.template registerEvaluator<EvalT>((*evaluators)[i]);
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
