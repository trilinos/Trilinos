#ifndef PANZER_EQUATIONSET_DEFAULT_IMPL_T_HPP
#define PANZER_EQUATIONSET_DEFAULT_IMPL_T_HPP

// ***********************************************************************
template <typename EvalT>
panzer::EquationSet_DefaultImpl<EvalT>::
EquationSet_DefaultImpl(const panzer::InputEquationSet& ies,
			const panzer::CellData& cell_data) :
  m_input_eq_set(ies),
  m_cell_data(cell_data)
{ }

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

#endif
