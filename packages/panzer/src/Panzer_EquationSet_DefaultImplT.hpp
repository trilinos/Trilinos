#ifndef PANZER_EQUATIONSET_DEFAULT_IMPL_T_H
#define PANZER_EQUATIONSET_DEFAULT_IMPL_T_H

// ***********************************************************************
template <typename EvalT>
panzer::EquationSet_DefaultImpl<EvalT>::
EquationSet_DefaultImpl(const panzer::InputEquationSet& ies,
			const panzer::CellData& cell_data) :
  m_input_eq_set(ies),
  m_cell_data(cell_data)
{ }

// ***********************************************************************
/*
template <typename EvalT>
void panzer::EquationSet<EvalT>::
buildAndRegisterMaterialModelEvaluators(int physics_id, 
					PHX::FieldManager<panzer::Traits>& fm,
					const std::vector<std::pair<std::string,Teuchos::RCP<panzer::Basis> > > & dofs) const
{
  
  Teuchos::RCP< std::vector< Teuchos::RCP<PHX::Evaluator<panzer::Traits> > > > evaluators;
  
  panzer::MaterialModelEvaluatorFactoryHandle<EvalT> 
    factory(m_input_eq_set.model_factory);
  
    evaluators = factory.buildEvaluators(m_input_eq_set, entries, 
					 *default_param_list);
    
    // Loop over evaluators and register them with field manager
    for (std::size_t i=0; i < evaluators->size(); ++i)
      fm.template registerEvaluator<EvalT>((*evaluators)[i]);
    
}
*/

// ***********************************************************************

#endif
