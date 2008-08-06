// @HEADER
// @HEADER

#include <sstream>
#include <typeinfo>
#include "Phalanx_DataLayout_Generic.hpp"

//**********************************************************************
PHX::Generic::Generic(const std::string& unique_identifier, std::size_t size) :
  m_name(unique_identifier),
  m_size(size)
{ }

//**********************************************************************
PHX::Generic::~Generic()
{ }

//**********************************************************************
bool PHX::Generic::operator==(const PHX::DataLayout& right) const
{
  const PHX::Generic* tmp = 0;
  tmp = dynamic_cast< const PHX::Generic* >(&right);

  if (tmp == 0)
    return false;

  return (  (this->name() == tmp->name()) &&
	    (this->size() == tmp->size()) );
}

//**********************************************************************
const std::string& PHX::Generic::name() const
{ return m_name; }

//**********************************************************************
std::size_t PHX::Generic::size() const
{ return m_size; }

//**********************************************************************
const std::string PHX::Generic::identifier() const
{ 
  std::ostringstream ost;
  ost << this->name() << this->size();
  return ost.str(); 
}

//**********************************************************************
void PHX::Generic::print(std::ostream& os, int indent) const
{
  std::ostringstream s;
  for (int i = 0; i < indent; i++)
    s << " ";

  os << s.str() << m_name << ", size = " << m_size;
}

//**********************************************************************
std::ostream& PHX::operator<<(std::ostream& os, const PHX::Generic& v)
{
  v.print(os);
  return os;
}

//**********************************************************************
