
#ifndef TEUCHOS_PARAMETER_LIST_H
#define TEUCHOS_PARAMETER_LIST_H

#include "Teuchos_ParameterEntry.hpp" // class data element 
#include "Teuchos_TestForException.hpp"
#include "Teuchos_map.hpp"

namespace Teuchos {

//! Manipulating lists of parameters.
class ParameterList {

  //! Parameter container typedef
  typedef Teuchos::map<string, ParameterEntry> Map;

  //! Parameter container const iterator typedef
  typedef Map::const_iterator ConstIterator;

  //! Parameter container iterator typedef
  typedef Map::iterator Iterator;

public:

  //! Constructor
  ParameterList();

  //! Copy Constructor
  ParameterList(const ParameterList& source);

  //! Copy
  ParameterList& operator=(const ParameterList& source);

  //! Deconstructor
  ~ParameterList();

  //! ParameterList %unused parameters
  void unused() const;

  //! Creates and empty sublist and returns a reference to the
  //! sublist. If the list already exists, returns reference to that
  //! sublist. If the name exists but is not a sublist, throws an error.
  ParameterList& sublist(const string& name);

  //! Returns a const reference to the sublist
  /*! 
    If the list does not already exist, throws an error. If the name
    exists but is not a sublist, throws an error.
  */
  const ParameterList& sublist(const string& name) const;

  /** @name Setting Parameters 
   
    Sets different types of parameters. The type depends on the second
    entry. Be sure to use static_cast<type>() when the type is
    ambiguous. Both char* and string map to are stored as strings
    internally. Sets the parameter as "unused".
  */
  //@{
  template<typename T>
  void setParameter(const string& name, T value);

  // Handles the case when a user sets the parameter with a character
  // string in parenthesis.
  void setParameter(const string& name, const char value[])
	{ setParameter(name, std::string(value) ); }
  //@}

  /** @name Getting Parameters 
   
    Get different types of parameters. The type depends on the second
    entry. Returns the nominal value if that parameter has not been
    specified. The non-const version adds the (name, nominal) pair to
    the list if it's not already specified. Be sure to use
    static_cast<type>() when the type is ambiguous. Both char* and
    string map return string values. Sets the parameters as "used".
  */
  //@{
  template<typename T>
  T& getParameter(const string& name, T def_value);

  std::string& getParameter(const string& name, const char def_value[])
  	{ return getParameter(name, std::string(def_value)); }

  template<typename T>
  const T& getParameter(const string& name, T def_value) const;

  const std::string& getParameter(const string& name, const char def_value[]) const
	{ return getParameter(name, std::string(def_value)); }

  // The parameter should already be set, so these methods return it or throw
  // an exception if the parameter doesn't exist.
  template<typename T>
  T& getParameter(const string& name);

  template<typename T>
  const T& getParameter(const string& name) const;
  //@}

  //! Return true if a parameter with this name exists.
  bool isParameter(const string& name) const;

  bool isParameterSublist(const string& name) const;

  // Needs to be called like PL.template isParameterType<int>( "Iters" )
  template<typename T>
  bool isParameterType(const string& name) const;

  template<typename T>
  bool isParameterType(const string& name, T* ptr) const;

  //! Printing 
  ostream& print(ostream& os, int indent = 0) const;

private:

  //! Access to name (i.e., returns i->first)
  const string& name(ConstIterator i) const;

  //! Access to ParameterEntry (i.e., returns i->second)
  ParameterEntry& entry(Iterator i);

  //! Access to ParameterEntry (i.e., returns i->second)
  const ParameterEntry& entry(ConstIterator i) const;

private:

  //! Parameter list
  Map params_;
 
  //! Used to create a string when the getParameter is called with a
  //! char* nominal value. A new string is created for each such
  //! argument. The whole group of strings is destroyed when this object
  //! is destroyed. This is really annoying, but I don't know a better 
  //! way.
  mutable vector<string> tmpstrings;
};


template<typename T>
void ParameterList::setParameter(const string& name, T value)
{
  params_[name].setValue(value);
}

template<typename T>
T& ParameterList::getParameter(const string& name, T def_value)
{
  ConstIterator i = params_.find(name);

  // The parameter was not found, add it to the list
  if (i == params_.end()) {
    params_[name].setValue(def_value, true);
    i = params_.find(name);
  }
  // Return the value of the parameter
  return getValue<T>(entry(i));
}

template<typename T>
const T& ParameterList::getParameter(const string& name, T def_value) const
{
  ConstIterator i = params_.find(name);

  // This parameter was not found, add it to the list
  if (i != params_.end() && isParameterType(name, &def_value)) {
	return getValue<T>(entry(i));
  }

  // Return the value of the default parameter if not found.
  return def_value;
}

template<typename T>
T& ParameterList::getParameter(const string& name) 
{
  ConstIterator i = params_.find(name);

  // This parameter was not found or is wrong type, throw an exception
  TEST_FOR_EXCEPTION( i == params_.end(), std::runtime_error,
	"getParameter ( " << name << " ) failed -- parameter does not exist! " );
  TEST_FOR_EXCEPTION( !isParameterType( name, (T*)NULL ), std::runtime_error,
	"getParameter ( " << name << " ) failed -- parameter is wrong type! " );
	 
  // Return the value.
  return getValue<T>(entry(i));
}

template<typename T>
const T& ParameterList::getParameter(const string& name) const
{
  ConstIterator i = params_.find(name);

  // This parameter was not found, throw and exception
  TEST_FOR_EXCEPTION( i == params_.end(), std::runtime_error,
	"getParameter ( " << name << " ) failed -- parameter does not exist! " );
  TEST_FOR_EXCEPTION( !isParameterType( name, (T*)NULL ), std::runtime_error,
	"getParameter ( " << name << " ) failed -- parameter is wrong type! " );

  // Return the default value for this type
  return getValue<T>(entry(i));
}

template<typename T>
bool ParameterList::isParameterType(const string& name, T* ptr) const
{
  ConstIterator i = params_.find(name);
  
  // If parameter doesn't exist, return false.
  if (i == params_.end()) 
    return false;
  // Try to cast the parameter to the type we think it should be.
  try {
    getValue<T>(entry(i));
  }
  catch( exception& e ) {
    return false;
  }
  // If no exception was thrown, we should be OK.
  return true;
}

template<typename T>
bool ParameterList::isParameterType(const string& name) const
{
  ConstIterator i = params_.find(name);
  
  // If parameter doesn't exist, return false.
  if (i == params_.end()) 
    return false;
  // Try to cast the parameter to the type we think it should be.
  try {
    getValue<T>(entry(i));
  }
  catch( exception& e ) {
    return false;
  }
  // If no exception was thrown, we should be OK.
  return true;
}
 
template<typename T>
bool isParameterType( ParameterList& l, const string& name )
{
  return l.isParameterType( name, (T*)NULL );
}

inline ostream& operator<<(ostream& os, const ParameterList& l)
{
  return l.print(os);
}

} // end of Teuchos namespace

#endif


