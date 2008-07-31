#ifndef PHX_FIELD_MANAGER_DEF_HPP
#define PHX_FIELD_MANAGER_DEF_HPP

#include "Sacado_mpl_size.hpp"
#include "boost/mpl/at.hpp"
#include "Phalanx_EvaluationContainer_TemplateBuilder.hpp"

// **************************************************************
template<typename Traits>
inline
PHX::FieldManager<Traits>::FieldManager() :
  max_num_cells_(-1)
{
  num_scalar_types_ = Sacado::mpl::size<typename Traits::ScalarTypes>::value;
  PHX::ScalarContainer_TemplateBuilder<Traits> builder;
  scalar_containers_.buildObjects(builder);
}

// **************************************************************
template<typename Traits>
inline
PHX::FieldManager<Traits>::~FieldManager()
{ }

// **************************************************************
template<typename Traits>
template<typename DataT> 
inline
void PHX::FieldManager<Traits>::
getFieldData(PHX::Field<DataT>& h)
{
  h.setFieldData(scalar_containers_.template 
    getAsObject< typename boost::mpl::at<typename Traits::DataToScalarMap,
	    DataT>::type >()->template getFieldData<DataT>(h.fieldTag()) );
}
    
// **************************************************************
template<typename Traits>
template<typename DataT> 
inline
void PHX::FieldManager<Traits>::
getFieldData(const PHX::FieldTag& v, Teuchos::ArrayRCP<DataT>& d)
{
  d = scalar_containers_.template 
    getAsObject< typename boost::mpl::at<typename Traits::DataToScalarMap,
    DataT>::type >()->template getFieldData<DataT>(v);
}

// **************************************************************
template<typename Traits>
inline
void PHX::FieldManager<Traits>::
requireFieldForAllTypes(const PHX::FieldTag& v)
{
  typedef PHX::ScalarContainer_TemplateManager<Traits> SCTM;
  
  typename SCTM::iterator it = scalar_containers_.begin();
  for (; it != scalar_containers_.end(); ++it) {
    it->requireField(v);
  }
}

// **************************************************************
template<typename Traits>
template<typename ScalarT>
inline
void PHX::FieldManager<Traits>::
requireFieldForScalarType(const PHX::FieldTag& v)
{
  scalar_containers_.template getAsBase<ScalarT>()->template requireField(v);
}
    
// **************************************************************
template<typename Traits>
inline
void PHX::FieldManager<Traits>::
registerEvaluatorForAllTypes(const Teuchos::RCP<PHX::Evaluator<Traits> >& p)
{
  typedef PHX::ScalarContainer_TemplateManager<Traits> SCTM;
  
  typename SCTM::iterator it = scalar_containers_.begin();
  for (; it != scalar_containers_.end(); ++it) {
    it->registerEvaluator(p);
  }
}

// **************************************************************
template<typename Traits>
template<typename ScalarT>
inline
void PHX::FieldManager<Traits>::
registerEvaluatorForScalarType(const Teuchos::RCP<PHX::Evaluator<Traits> >& p)
{
  scalar_containers_.template getAsBase<ScalarT>()->template registerEvaluator(p);
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
  max_num_cells_ = max_num_cells;

  typedef PHX::ScalarContainer_TemplateManager<Traits> SCTM;
  typename SCTM::iterator it = scalar_containers_.begin();
  for (; it != scalar_containers_.end(); ++it)
    it->postRegistrationSetup(max_num_cells_, *this);
}

// **************************************************************
template<typename Traits>
template<typename ScalarT>
inline
void PHX::FieldManager<Traits>::
evaluateFields(typename Traits::EvalData d)
{
  scalar_containers_.template getAsBase<ScalarT>()->template evaluateFields(d);
}

// **************************************************************
template<typename Traits>
template<typename ScalarT>
inline
void PHX::FieldManager<Traits>::
preEvaluate(typename Traits::PreEvalData d)
{
  scalar_containers_.template getAsBase<ScalarT>()->template preEvaluate(d);
}

// **************************************************************
template<typename Traits>
template<typename ScalarT>
inline
void PHX::FieldManager<Traits>::
postEvaluate(typename Traits::PostEvalData d)
{
  scalar_containers_.template getAsBase<ScalarT>()->template postEvaluate(d);
}

// **************************************************************
template<typename Traits>
inline
std::size_t PHX::FieldManager<Traits>::getMaxNumCells() const
{
  return max_num_cells_;
}

// **************************************************************
template<typename Traits>
inline
typename PHX::FieldManager<Traits>::iterator 
PHX::FieldManager<Traits>::begin()
{
  return scalar_containers_.begin();
}

// **************************************************************
template<typename Traits>
inline
typename PHX::FieldManager<Traits>::iterator 
PHX::FieldManager<Traits>::end()
{
  return scalar_containers_.end();
}

// **************************************************************
template<typename Traits>
inline
void PHX::FieldManager<Traits>::print(std::ostream& os) const
{
  typedef PHX::ScalarContainer_TemplateManager<Traits> SCTM;
  typename SCTM::const_iterator it = scalar_containers_.begin();
  for (; it != scalar_containers_.end(); ++it)
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
