#ifndef PHX_FIELD_DEF_H
#define PHX_FIELD_DEF_H

#include "Teuchos_TestForException.hpp"

//**********************************************************************
#ifdef PHX_DEBUG
template<typename DataT>
const std::string PHX::Field<DataT>::field_tag_error_msg = 
    "Error - PHX::Field::fieldTag() - No tag has been set!";
template<typename DataT>
const std::string PHX::Field<DataT>::field_data_error_msg = "Error - PHX::Field::operator[] - No data has been set!  Please call setFieldData(this) on all PHX::Field objects in providers!";
#endif

//**********************************************************************
template<typename DataT>
PHX::Field<DataT>::Field(const std::string& name, 
			 const Teuchos::RCP<PHX::DataLayout>& t) :
  tag(name,t)
#ifdef PHX_DEBUG
  , tag_set(true),
  data_set(false)
#endif
{ }

//**********************************************************************
template<typename DataT>
PHX::Field<DataT>::Field(const PHX::FieldTag& v) :
  tag(v)
#ifdef PHX_DEBUG
  ,tag_set(true),
  data_set(false)
#endif
{ }

//**********************************************************************
template<typename DataT>
PHX::Field<DataT>::Field() :
  tag("???", Teuchos::null)
#ifdef PHX_DEBUG
  ,tag_set(false),
  data_set(false)
#endif
{ }

//**********************************************************************
template<typename DataT>
PHX::Field<DataT>::~Field()
{ }

//**********************************************************************
template<typename DataT>
inline
const PHX::FieldTag& PHX::Field<DataT>::fieldTag() const
{ 
#ifdef PHX_DEBUG
  TEST_FOR_EXCEPTION(!tag_set, std::logic_error, field_tag_error_msg);
#endif
  return tag;
}

//**********************************************************************
template<typename DataT>
inline
DataT& PHX::Field<DataT>::operator[](int index)
{ 
#ifdef PHX_DEBUG
  TEST_FOR_EXCEPTION(!data_set, std::logic_error, field_data_error_msg);
#endif
  return field_data[index];
}

//**********************************************************************
template<typename DataT>
inline
typename Teuchos::ArrayRCP<DataT>::Ordinal PHX::Field<DataT>::size()
{ 
#ifdef PHX_DEBUG
  TEST_FOR_EXCEPTION(!data_set, std::logic_error, field_data_error_msg);
#endif
  return field_data.size();
}

//**********************************************************************
template<typename DataT>
void PHX::Field<DataT>::setFieldTag(const PHX::FieldTag& v)
{  
#ifdef PHX_DEBUG
  tag_set = true;
#endif
  tag = v;
}

//**********************************************************************
template<typename DataT>
void PHX::Field<DataT>::setFieldData(const Teuchos::ArrayRCP<DataT>& d)
{ 
#ifdef PHX_DEBUG
  data_set = true;
#endif
  field_data = d;
}

//**********************************************************************
template<typename DataT>
void PHX::Field<DataT>::print(std::ostream& os) const
{
  os << "Printing Field: \n" << tag << std::endl;
  typedef typename Teuchos::ArrayRCP<DataT>::Ordinal size_type;
  for (size_type i = 0; i < field_data.size(); ++i)
    os << "value[" << i << "] = " << field_data[i] << std::endl;
}

//**********************************************************************
template<typename DataT>
std::ostream& PHX::operator<<(std::ostream& os, const PHX::Field<DataT>& f)
{
  f.print(os);
  return os;
}

//**********************************************************************

#endif
