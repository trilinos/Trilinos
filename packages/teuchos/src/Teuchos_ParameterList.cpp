#include "Teuchos_ParameterList.hpp"	// class definition

/* NOTE: ASCI Red (TFLOP) does not support the i-> function for iterators 
 * in the STL.  Therefore when compiling for the TFLOP we must redefine the 
 * iterator from i-> to (*i). This slows things down on other platforms 
 * so we switch between the two when necessary.
 */
using namespace Teuchos;

ParameterList::ParameterList() {}

ParameterList::ParameterList(const ParameterList& source) 
{
  params_ = source.params_;
}

ParameterList& ParameterList::operator=(const ParameterList& source) 
{
  if (&source == this)
    return *this;

  params_ = source.params_;
  return *this;
}

ParameterList::~ParameterList() 
{
}

void ParameterList::unused() const
{
  for (ConstIterator i = params_.begin(); i != params_.end(); ++i) {
    if (!(entry(i).isUsed())) {
      cout << "WARNING: Parameter \"" << name(i) << "\" " << entry(i)
	   << " is unused" << endl;
    }
  }
}

bool ParameterList::isParameterSublist(const string& name) const
{
  ConstIterator i = params_.find(name);

  if (i != params_.end())
    return (entry(i).isList());

  return false;
}

bool ParameterList::isParameter(const string& name) const
{
  return (params_.find(name) != params_.end());
}

ParameterList& ParameterList::sublist(const string& name)
{
  // Find name in list, if it exists.
  Iterator i = params_.find(name);

  // If it does exist and is a list, return the list value.
  // Otherwise, throw an error.
  if (i != params_.end()) {
    if (entry(i).isList())
      return getValue<ParameterList>(entry(i));
    else
    {
      cerr << "ERROR: Parameter " << name << " is not a valid list." << endl;
      throw "Teuchos Error";
    } 
  }

  // If it does not exist, create a new empty list and return a reference
  return params_[name].setList(true);
}

const ParameterList& ParameterList::sublist(const string& name) const
{
  // Find name in list, if it exists.
  ConstIterator i = params_.find(name);

  // If it does not exist, throw an error
  if (i == params_.end()) {
    cerr << "ERROR: Parameter " << name << " is not a valid list." << endl;
    throw "Teuchos Error";
  }

  // If it does exist and is a list, return the list value.
  if (entry(i).isList())
    return getValue<ParameterList>(entry(i));
  else {
  // Otherwise, the parameter exists but is not a list. Throw an error.
  cerr << "ERROR: Parameter " << name << " is not a list." << endl;
  throw "Teuchos Error";
  }
}
  
ostream& ParameterList::print(ostream& os, int indent) const
{
  if (params_.begin() == params_.end()) 
  {
    for (int j = 0; j < indent; j ++)
      os << ' ';
    os << "[empty list]" << endl;
  }
  else 
    for (ConstIterator i = params_.begin(); i != params_.end(); ++i) 
    {
      for (int j = 0; j < indent; j ++)
	os << ' ';
      if (entry(i).isList()) 
      {
	os << name(i) << " -> " << endl;
	getValue<ParameterList>(entry(i)).print(os, indent + 2);
      }
      else
	os << name(i) << " = " << entry(i) << endl;
    }
  return os;
}


#if defined(TFLOP)

const string& ParameterList::name(ConstIterator i) const
{
  return ((*i).first);
}

ParameterEntry& ParameterList::entry(Iterator i)
{
  return ((*i).second);
}

const ParameterEntry& ParameterList::entry(ConstIterator i) const
{
  return ((*i).second);
}

#else

const string& ParameterList::name(ConstIterator i) const
{
  return (i->first);
}

ParameterEntry& ParameterList::entry(Iterator i)
{
  return (i->second);
}

const ParameterEntry& ParameterList::entry(ConstIterator i) const
{
  return (i->second);
}

#endif

