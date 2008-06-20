#ifndef PHX_FIELD_EVALUATOR_UTILITIES_DEF_H
#define PHX_FIELD_EVALUATOR_UTILITIES_DEF_H

#include <algorithm>
#include <vector>

//**********************************************************************
template<typename Traits>
PHX::FieldEvaluatorUtilities<Traits>::FieldEvaluatorUtilities() :
  name_("???")
{ }

//**********************************************************************
template<typename Traits>
PHX::FieldEvaluatorUtilities<Traits>::~FieldEvaluatorUtilities()
{ }

//**********************************************************************
template<typename Traits>
void PHX::FieldEvaluatorUtilities<Traits>::
addEvaluatedField(const PHX::FieldTag& v)
{ 
  std::vector<FieldTag>::iterator test = 
    std::find(evaluated_.begin(), evaluated_.end(), v);
  
  if ( test == evaluated_.end() )
    evaluated_.push_back(v);
}

//**********************************************************************
template<typename Traits>
template<typename DataT>
void PHX::FieldEvaluatorUtilities<Traits>::
addEvaluatedField(const PHX::Field<DataT>& f)
{ 
  this->template addEvaluatedField(f.fieldTag());
}

//**********************************************************************
template<typename Traits>
void PHX::FieldEvaluatorUtilities<Traits>::
addDependentField(const PHX::FieldTag& v)
{
  std::vector<FieldTag>::iterator test = 
    std::find(required_.begin(), required_.end(), v);
  
  if ( test == required_.end() )
    required_.push_back(v);
}

//**********************************************************************
template<typename Traits>
template<typename DataT>
void PHX::FieldEvaluatorUtilities<Traits>::
addDependentField(const PHX::Field<DataT>& v)
{
  this->template addDependentField(v.fieldTag());
}

//**********************************************************************
template<typename Traits>
void PHX::FieldEvaluatorUtilities<Traits>::
setName(const std::string& name)
{ name_ = name; }

//**********************************************************************
template<typename Traits>
const std::vector<PHX::FieldTag>&
PHX::FieldEvaluatorUtilities<Traits>::evaluatedFields() const
{ return evaluated_; }

//**********************************************************************
template<typename Traits>
const std::vector<PHX::FieldTag>&
PHX::FieldEvaluatorUtilities<Traits>::dependentFields() const
{ return required_; }

//**********************************************************************
template<typename Traits>
void PHX::FieldEvaluatorUtilities<Traits>::
preEvaluate(typename Traits::PreEvalData d)
{ }

//**********************************************************************
template<typename Traits>
void PHX::FieldEvaluatorUtilities<Traits>::
postEvaluate(typename Traits::PostEvalData d)
{ }

//**********************************************************************
template<typename Traits>
const std::string& PHX::FieldEvaluatorUtilities<Traits>::
getName() const
{return name_;}

//**********************************************************************

#endif
