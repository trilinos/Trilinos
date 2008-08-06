// @HEADER
// @HEADER

#ifndef PHX_FIELDTAG_TAG_DEF_HPP
#define PHX_FIELDTAG_TAG_DEF_HPP

#include <sstream>
#include "Teuchos_TestForException.hpp"

//**********************************************************************
template<typename DataT>
PHX::Tag<DataT>::Tag(const std::string& name,
		     const Teuchos::RCP<PHX::DataLayout>& dl) :
  m_name(name),
  m_data_layout(dl)
{ }

//**********************************************************************
template<typename DataT>
PHX::Tag<DataT>::~Tag()
{ }

//**********************************************************************
template<typename DataT>
Teuchos::RCP<PHX::FieldTag> PHX::Tag<DataT>::clone() const
{
  Teuchos::RCP<PHX::FieldTag> tag = 
    Teuchos::rcp(new PHX::Tag<DataT>(m_name, m_data_layout));
  return (tag);
}

//**********************************************************************
template<typename DataT>
void PHX::Tag<DataT>::operator=(const PHX::Tag<DataT>& t)
{
  m_name = t.m_name;
  m_data_layout = t.m_data_layout;
}

//**********************************************************************
template<typename DataT>
bool PHX::Tag<DataT>::operator==(const PHX::FieldTag& t) const
{
  return (  (this->name() == t.name()) &&
	    (this->dataLayout() == t.dataLayout()) &&
	    (this->dataTypeInfo() == t.dataTypeInfo()) );
}

//**********************************************************************
template<typename DataT>
const std::string& PHX::Tag<DataT>::name() const
{ return m_name; }

//**********************************************************************
template<typename DataT>
const PHX::DataLayout& PHX::Tag<DataT>::dataLayout() const
{ return *m_data_layout; }

//**********************************************************************
template<typename DataT>
const std::type_info& PHX::Tag<DataT>::dataTypeInfo() const
{ 
  DataT tmp;
  return typeid(tmp);
}

//**********************************************************************
template<typename DataT>
const std::string PHX::Tag<DataT>::identifier() const
{
  std::ostringstream ost;

  ost << this->name() 
      << this->dataTypeInfo().name() 
      << this->dataLayout().identifier();
  return ost.str(); 
}

//**********************************************************************
template<typename DataT>
void PHX::Tag<DataT>::print(std::ostream& os) const
{
  DataT tmp;
  os << "Tag:  " << m_name << ", " << typeid(tmp).name()
     << ", DataLayout: " << *m_data_layout;

}

//**********************************************************************
template<typename DataT>
std::ostream& PHX::operator<<(std::ostream& os, const PHX::Tag<DataT>& v)
{
  v.print(os);
  return os;
}

//**********************************************************************

#endif
