
#include "Teuchos_ParameterEntry.hpp" // class definition
#include "Teuchos_ParameterList.hpp"	 // for sublists

using namespace Teuchos;

ParameterEntry::ParameterEntry() : 
  isUsed_(false),
  isDefault_(false),
  isList_(false)
{
}

ParameterEntry::ParameterEntry(const ParameterEntry& source)
{
  operator=(source);
}

ParameterEntry& ParameterEntry::operator=(const ParameterEntry& source)
{
  if (&source == this)
    return *this;

  val_ = source.val_;
  isUsed_ = source.isUsed_;
  isDefault_ = source.isDefault_;
  isList_ = source.isList_;

  return *this;
}

ParameterList& ParameterEntry::setList(bool isDefault)
{
  val_ = ParameterList();
  isDefault_ = isDefault;
  isUsed_ = true;
  isList_ = true;
  return any_cast<ParameterList>( val_ );
}

void ParameterEntry::reset()
{
  //delete val_;
  isUsed_ = false;
  isDefault_ = false;
}

ostream& ParameterEntry::leftshift(ostream& os) const
{
  if( !isList_ ) os << val_;

  if (isDefault_)
    os << "   [default]";
  else if (!isUsed_)
    os << "   [unused]";

  return os;
}



