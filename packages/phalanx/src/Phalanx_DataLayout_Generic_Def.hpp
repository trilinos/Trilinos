#ifndef PHX_DATA_LAYOUT_GENERIC_DEF
#define PHX_DATA_LAYOUT_GENERIC_DEF

#include <sstream>
#include <typeinfo>

//**********************************************************************
template<typename Entity>
PHX::Generic<Entity>::
Generic(const std::string& unique_identifier, std::size_t size) :
  name_(unique_identifier),
  size_(size)
{ }

//**********************************************************************
template<typename Entity>
PHX::Generic<Entity>::~Generic()
{ }

//**********************************************************************
template<typename Entity>
bool PHX::Generic<Entity>::operator==(const PHX::DataLayout& right) const
{
  const PHX::Generic<Entity>* tmp = 0;
  tmp = dynamic_cast< const PHX::Generic<Entity>* >(&right);

  if (tmp == 0)
    return false;

  return (  (this->name() == tmp->name()) &&
	    (this->size() == tmp->size()) &&
	    (this->getAlgebraicTypeInfo() == tmp->getAlgebraicTypeInfo()) );
}

//**********************************************************************
template<typename Entity>
const std::string& PHX::Generic<Entity>::name() const
{ return name_; }

//**********************************************************************
template<typename Entity>
std::size_t PHX::Generic<Entity>::size() const
{ return size_; }

//**********************************************************************
template<typename Entity>
const std::type_info& PHX::Generic<Entity>::getAlgebraicTypeInfo() const
{ 
  Entity tmp;
  return typeid(tmp);
}

//**********************************************************************
template<typename Entity>
void PHX::Generic<Entity>::print(std::ostream& os, int indent) const
{
  std::ostringstream s;
  for (int i = 0; i < indent; i++)
    s << " ";

  os << s.str() << name_ << ", size = " << size_ << ", " 
     << Entity::name;
}

//**********************************************************************
template<typename Entity>
std::ostream& PHX::operator<<(std::ostream& os, const PHX::Generic<Entity>& v)
{
  v.print(os);
  return os;
}

//**********************************************************************

#endif
