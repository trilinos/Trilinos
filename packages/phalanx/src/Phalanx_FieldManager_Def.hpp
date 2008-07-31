#ifndef PHX_FIELD_MANAGER_DEF_HPP
#define PHX_FIELD_MANAGER_DEF_HPP

#include "Sacado_mpl_size.hpp"
#include "boost/mpl/at.hpp"
#include "Phalanx_EvaluationContainer_TemplateBuilder.hpp"

// **************************************************************
template<typename Traits>
inline
PHX::FieldManager<Traits>::FieldManager() :
  m_max_num_cells(-1)
{
  m_num_evaluation_types = 
    Sacado::mpl::size<typename Traits::EvalTypes>::value;
  PHX::EvaluationContainer_TemplateBuilder<Traits> builder;
  m_eval_containers.buildObjects(builder);
}

// **************************************************************
template<typename Traits>
inline
PHX::FieldManager<Traits>::~FieldManager()
{ }

// **************************************************************
template<typename Traits>
template<typename DataT, typename EvalT> 
inline
void PHX::FieldManager<Traits>::
getFieldData(PHX::Field<DataT>& h)
{
  h.setFieldData(m_eval_containers.template 
    getAsObject<EvalT>()->template getFieldData<DataT>(h.fieldTag()) );
}
    
// **************************************************************
template<typename Traits>
template<typename DataT, typename EvalT> 
inline
void PHX::FieldManager<Traits>::
getFieldData(const PHX::FieldTag& v, Teuchos::ArrayRCP<DataT>& d)
{
  d = m_eval_containers.template 
    getAsObject<EvalT>()->template getFieldData<DataT>(v);
}

// **************************************************************
template<typename Traits>
inline
void PHX::FieldManager<Traits>::
requireFieldForAllTypes(const PHX::FieldTag& v)
{
  typedef PHX::EvaluationContainer_TemplateManager<Traits> SCTM;
  
  typename SCTM::iterator it = m_eval_containers.begin();
  for (; it != m_eval_containers.end(); ++it) {
    it->requireField(v);
  }
}

// **************************************************************
template<typename Traits>
template<typename EvalT>
inline
void PHX::FieldManager<Traits>::
requireFieldForScalarType(const PHX::FieldTag& v)
{
  m_eval_containers.template getAsBase<EvalT>()->template requireField(v);
}
    
// **************************************************************
template<typename Traits>
inline
void PHX::FieldManager<Traits>::
registerEvaluatorForAllTypes(const Teuchos::RCP<PHX::Evaluator<Traits> >& p)
{
  typedef PHX::EvaluationContainer_TemplateManager<Traits> SCTM;
  
  typename SCTM::iterator it = m_eval_containers.begin();
  for (; it != m_eval_containers.end(); ++it) {
    it->registerEvaluator(p);
  }
}

// **************************************************************
template<typename Traits>
template<typename EvalT>
inline
void PHX::FieldManager<Traits>::
registerEvaluatorForScalarType(const Teuchos::RCP<PHX::Evaluator<Traits> >& p)
{
  m_eval_containers.template getAsBase<EvalT>()->template registerEvaluator(p);
}

// **************************************************************
template<typename Traits>
inline
void PHX::FieldManager<Traits>::
registerEvaluatorForScalarType(FieldManager::iterator it,
			const Teuchos::RCP<PHX::Evaluator<Traits> >& e)
{
  it->registerEvaluator(e);
}

// **************************************************************
template<typename Traits>
inline
void PHX::FieldManager<Traits>::
postRegistrationSetup(std::size_t max_num_cells)
{
  m_max_num_cells = max_num_cells;

  typedef PHX::EvaluationContainer_TemplateManager<Traits> SCTM;
  typename SCTM::iterator it = m_eval_containers.begin();
  for (; it != m_eval_containers.end(); ++it)
    it->postRegistrationSetup(m_max_num_cells, *this);
}

// **************************************************************
template<typename Traits>
template<typename EvalT>
inline
void PHX::FieldManager<Traits>::
evaluateFields(typename Traits::EvalData d)
{
  m_eval_containers.template getAsBase<EvalT>()->template evaluateFields(d);
}

// **************************************************************
template<typename Traits>
template<typename EvalT>
inline
void PHX::FieldManager<Traits>::
preEvaluate(typename Traits::PreEvalData d)
{
  m_eval_containers.template getAsBase<EvalT>()->template preEvaluate(d);
}

// **************************************************************
template<typename Traits>
template<typename EvalT>
inline
void PHX::FieldManager<Traits>::
postEvaluate(typename Traits::PostEvalData d)
{
  m_eval_containers.template getAsBase<EvalT>()->template postEvaluate(d);
}

// **************************************************************
template<typename Traits>
inline
std::size_t PHX::FieldManager<Traits>::getMaxNumCells() const
{
  return m_max_num_cells;
}

// **************************************************************
template<typename Traits>
inline
typename PHX::FieldManager<Traits>::iterator 
PHX::FieldManager<Traits>::begin()
{
  return m_eval_containers.begin();
}

// **************************************************************
template<typename Traits>
inline
typename PHX::FieldManager<Traits>::iterator 
PHX::FieldManager<Traits>::end()
{
  return m_eval_containers.end();
}

// **************************************************************
template<typename Traits>
inline
void PHX::FieldManager<Traits>::print(std::ostream& os) const
{
  typedef PHX::EvaluationContainer_TemplateManager<Traits> SCTM;
  typename SCTM::const_iterator it = m_eval_containers.begin();
  for (; it != m_eval_containers.end(); ++it)
    os << (*it);
}

// **************************************************************
template<typename Traits>
inline
std::ostream& PHX::operator<<(std::ostream& os, 
			      const PHX::FieldManager<Traits>& vm)
{
  vm.print(os);
  return os;
}

// **************************************************************

#endif
