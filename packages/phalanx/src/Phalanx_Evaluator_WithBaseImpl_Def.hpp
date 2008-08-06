#ifndef PHX_EVALUATOR_WITHBASEIMPL_DEF_H
#define PHX_EVALUATOR_WITHBASEIMPL_DEF_H

#include "Phalanx_ConfigDefs.hpp"
#include "Phalanx_FieldTag_STL_Functors.hpp"

//**********************************************************************
template<typename Traits>
PHX::EvaluatorWithBaseImpl<Traits>::
EvaluatorWithBaseImpl(const std::string& evaluator_name) :
  name_(evaluator_name)
{ }

//**********************************************************************
template<typename Traits>
PHX::EvaluatorWithBaseImpl<Traits>::EvaluatorWithBaseImpl() :
  name_("???")
{ }

//**********************************************************************
template<typename Traits>
PHX::EvaluatorWithBaseImpl<Traits>::~EvaluatorWithBaseImpl()
{ }

//**********************************************************************
template<typename Traits>
void PHX::EvaluatorWithBaseImpl<Traits>::
addEvaluatedField(const PHX::FieldTag& v)
{ 
  PHX::FTPredRef pred(v);
  std::vector< Teuchos::RCP<FieldTag> >::iterator test = 
    std::find_if(evaluated_.begin(), evaluated_.end(), pred);
  
  if ( test == evaluated_.end() )
    evaluated_.push_back(v.clone());
}

//**********************************************************************
template<typename Traits>
template<typename DataT>
void PHX::EvaluatorWithBaseImpl<Traits>::
addEvaluatedField(const PHX::Field<DataT>& f)
{ 
  this->template addEvaluatedField(f.fieldTag());
}

//**********************************************************************
template<typename Traits>
void PHX::EvaluatorWithBaseImpl<Traits>::
addDependentField(const PHX::FieldTag& v)
{
  PHX::FTPredRef pred(v);
  std::vector< Teuchos::RCP<FieldTag> >::iterator test = 
    std::find_if(required_.begin(), required_.end(), pred);
  
  if ( test == required_.end() )
    required_.push_back(v.clone());
}

//**********************************************************************
template<typename Traits>
template<typename DataT>
void PHX::EvaluatorWithBaseImpl<Traits>::
addDependentField(const PHX::Field<DataT>& v)
{
  this->template addDependentField(v.fieldTag());
}

//**********************************************************************
template<typename Traits>
void PHX::EvaluatorWithBaseImpl<Traits>::
setName(const std::string& name)
{ name_ = name; }

//**********************************************************************
template<typename Traits>
const std::vector< Teuchos::RCP<PHX::FieldTag> >&
PHX::EvaluatorWithBaseImpl<Traits>::evaluatedFields() const
{ return evaluated_; }

//**********************************************************************
template<typename Traits>
const std::vector< Teuchos::RCP<PHX::FieldTag> >&
PHX::EvaluatorWithBaseImpl<Traits>::dependentFields() const
{ return required_; }

//**********************************************************************
template<typename Traits>
void PHX::EvaluatorWithBaseImpl<Traits>::
preEvaluate(typename Traits::PreEvalData d)
{ }

//**********************************************************************
template<typename Traits>
void PHX::EvaluatorWithBaseImpl<Traits>::
postEvaluate(typename Traits::PostEvalData d)
{ }

//**********************************************************************
template<typename Traits>
const std::string& PHX::EvaluatorWithBaseImpl<Traits>::
getName() const
{return name_;}

//**********************************************************************

#endif
