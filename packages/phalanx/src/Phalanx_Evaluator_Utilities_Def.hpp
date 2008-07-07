#ifndef PHX_FIELD_EVALUATOR_UTILITIES_DEF_H
#define PHX_FIELD_EVALUATOR_UTILITIES_DEF_H

#include <algorithm>
#include <vector>

//**********************************************************************
template<typename Traits>
PHX::EvaluatorUtilities<Traits>::EvaluatorUtilities() :
  name_("???")
{ }

//**********************************************************************
template<typename Traits>
PHX::EvaluatorUtilities<Traits>::~EvaluatorUtilities()
{ }

//**********************************************************************
template<typename Traits>
void PHX::EvaluatorUtilities<Traits>::
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
void PHX::EvaluatorUtilities<Traits>::
addEvaluatedField(const PHX::Field<DataT>& f)
{ 
  this->template addEvaluatedField(f.fieldTag());
}

//**********************************************************************
template<typename Traits>
void PHX::EvaluatorUtilities<Traits>::
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
void PHX::EvaluatorUtilities<Traits>::
addDependentField(const PHX::Field<DataT>& v)
{
  this->template addDependentField(v.fieldTag());
}

//**********************************************************************
template<typename Traits>
void PHX::EvaluatorUtilities<Traits>::
setName(const std::string& name)
{ name_ = name; }

//**********************************************************************
template<typename Traits>
const std::vector<PHX::FieldTag>&
PHX::EvaluatorUtilities<Traits>::evaluatedFields() const
{ return evaluated_; }

//**********************************************************************
template<typename Traits>
const std::vector<PHX::FieldTag>&
PHX::EvaluatorUtilities<Traits>::dependentFields() const
{ return required_; }

//**********************************************************************
template<typename Traits>
void PHX::EvaluatorUtilities<Traits>::
preEvaluate(typename Traits::PreEvalData d)
{ }

//**********************************************************************
template<typename Traits>
void PHX::EvaluatorUtilities<Traits>::
postEvaluate(typename Traits::PostEvalData d)
{ }

//**********************************************************************
template<typename Traits>
const std::string& PHX::EvaluatorUtilities<Traits>::
getName() const
{return name_;}

//**********************************************************************

#endif
