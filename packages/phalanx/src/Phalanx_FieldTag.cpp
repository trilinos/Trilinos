#include "Phalanx_FieldTag.hpp"
#include "Teuchos_TestForException.hpp"

//**********************************************************************
PHX::FieldTag::FieldTag(const std::string& name,
			const Teuchos::RCP<PHX::DataLayout>& dl) :
  name_(name),
  data_layout_(dl)
{ }

//**********************************************************************
PHX::FieldTag::~FieldTag()
{ }

//**********************************************************************
PHX::FieldTag& PHX::FieldTag::operator=(const PHX::FieldTag& a)
{
  name_ = a.name();
  data_layout_ = a.dataLayout();
  return *this;
}

//**********************************************************************
bool PHX::FieldTag::operator==(const PHX::FieldTag& a) const
{
  return (  (this->name() == a.name()) &&
	    (*(this->dataLayout()) == *(a.dataLayout()))  );
}

//**********************************************************************
bool PHX::FieldTag::operator<(const PHX::FieldTag& a) const
{
  std::ostringstream left;
  std::ostringstream right;
  left << name() << *dataLayout(); 
  right << a.name() << *(a.dataLayout());

  return (left.str() < right.str());
}

//**********************************************************************
const std::string& PHX::FieldTag::name() const
{ return name_; }

//**********************************************************************
const Teuchos::RCP<PHX::DataLayout> PHX::FieldTag::dataLayout() const
{ return data_layout_; }

//**********************************************************************
void PHX::FieldTag::print(std::ostream& os, int indent) const
{
  std::ostringstream s;
  for (int i = 0; i < indent; i++)
    s << " ";

  os << s.str() << "FieldTag:  " << name_ << ", DataLayout: " << *data_layout_;

}

//**********************************************************************
std::ostream& PHX::operator<<(std::ostream& os, const PHX::FieldTag& v)
{
  v.print(os);
  return os;
}

//**********************************************************************
