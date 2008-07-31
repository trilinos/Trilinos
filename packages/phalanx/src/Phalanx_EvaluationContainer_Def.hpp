#ifndef PHX_SCALAR_CONTAINER_DEF_HPP
#define PHX_SCALAR_CONTAINER_DEF_HPP

#include "Phalanx_Evaluator.hpp"
#include "Phalanx_DebugStrings.hpp"

// *************************************************************************
template <typename ScalarT, typename Traits>
PHX::ScalarContainer<ScalarT, Traits>::ScalarContainer()
{
  this->vp_manager_.setScalarTypeName( PHX::getTypeString<ScalarT, Traits>() );
  this->data_container_template_manager_.buildObjects();
}

// *************************************************************************
template <typename ScalarT, typename Traits> 
PHX::ScalarContainer<ScalarT, Traits>::~ScalarContainer()
{

}

// *************************************************************************
template <typename ScalarT, typename Traits>
void PHX::ScalarContainer<ScalarT, Traits>::
requireField(const PHX::FieldTag& f)
{
  this->vp_manager_.requireField(f);
}

// *************************************************************************
template <typename ScalarT, typename Traits>
void PHX::ScalarContainer<ScalarT, Traits>::
registerEvaluator(const Teuchos::RCP<PHX::Evaluator<Traits> >& p)
{
  this->vp_manager_.registerEvaluator(p);
}

// *************************************************************************
template <typename ScalarT, typename Traits> 
void PHX::ScalarContainer<ScalarT, Traits>::
postRegistrationSetup(std::size_t max_num_cells,
		      PHX::FieldManager<Traits>& fm)
{
  // Figure out all evaluator dependencies
  if ( !(this->vp_manager_.sortingCalled()) )
    this->vp_manager_.sortAndOrderEvaluators();

  // Determine total amount of memory for all variables
  allocator_.reset();
  
  const std::vector<PHX::FieldTag>& var_list = 
    this->vp_manager_.getFieldTags();

  std::vector<PHX::FieldTag>::const_iterator  var = var_list.begin();
  for (; var != var_list.end(); ++var) {

    typename DCTM::iterator it = data_container_template_manager_.begin();
    for (; it != data_container_template_manager_.end(); ++it) {
      
      if (var->dataLayout()->getAlgebraicTypeInfo() == it->getAlgebraicTypeInfo()) {
	std::size_t size_of_data_type = it->getSizeOfDataType();
	int total_size = max_num_cells * size_of_data_type * var->dataLayout()->size();
	allocator_.addRequiredBytes(total_size);
      }
    }
  }

  allocator_.setup();

  // Allocate field data arrays
  //std::vector<PHX::FieldTag>::const_iterator  var = var_list.begin();
  for (var = var_list.begin(); var != var_list.end(); ++var) {

    typename DCTM::iterator it = data_container_template_manager_.begin();
    for (; it != data_container_template_manager_.end(); ++it) {
      
      if (var->dataLayout()->getAlgebraicTypeInfo() == it->getAlgebraicTypeInfo()) {
	it->allocateField(*var, max_num_cells, allocator_);
      }
    }
  }

  // Allow field evaluators to grab pointers to relevant field data
  this->vp_manager_.postRegistrationSetup(fm);
}

// *************************************************************************
template <typename ScalarT, typename Traits>
void PHX::ScalarContainer<ScalarT, Traits>::
evaluateFields(typename Traits::EvalData d)
{
  this->vp_manager_.evaluateFields(d);
}

// *************************************************************************
template <typename ScalarT, typename Traits>
void PHX::ScalarContainer<ScalarT, Traits>::
preEvaluate(typename Traits::PreEvalData d)
{
  this->vp_manager_.preEvaluate(d);
}

// *************************************************************************
template <typename ScalarT, typename Traits>
void PHX::ScalarContainer<ScalarT, Traits>::
postEvaluate(typename Traits::PostEvalData d)
{
  this->vp_manager_.postEvaluate(d);
}

// *************************************************************************
template <typename ScalarT, typename Traits> template <typename DataT>
Teuchos::ArrayRCP<DataT> 
PHX::ScalarContainer<ScalarT, Traits>::getFieldData(const PHX::FieldTag& f)
{
  Teuchos::ArrayRCP<DataT> r = 
    data_container_template_manager_.template getAsObject<DataT>()->
    getFieldData(f);
  return r;
}


// *************************************************************************
template <typename ScalarT, typename Traits>
void PHX::ScalarContainer<ScalarT, Traits>::print(std::ostream& os) const
{
  std::string type = PHX::getTypeString<ScalarT, Traits>();

  os << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << std::endl;
  os << "Starting PHX::ScalarContainer Output" << std::endl;
  os << "Scalar Type = " << type << std::endl;
  os << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << std::endl;
  os << this->vp_manager_ << std::endl;
  typename DCTM::const_iterator it = data_container_template_manager_.begin();
  for (; it != data_container_template_manager_.end(); ++it)
    os << *it << std::endl;
  os << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << std::endl;
  os << "Finished PHX::ScalarContainer Output" << std::endl;
  os << "Scalar Type = " << type << std::endl;
  os << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << std::endl;
  os << std::endl;
}

// *************************************************************************

#endif 
