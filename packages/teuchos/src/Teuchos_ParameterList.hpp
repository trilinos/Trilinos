
#ifndef TEUCHOS_PARAMETER_LIST_H
#define TEUCHOS_PARAMETER_LIST_H

#include "Teuchos_ParameterEntry.hpp" // class data element 

namespace Teuchos {

//! Manipulating lists of parameters.
class ParameterList {

  //! Parameter container typedef
  typedef map<string, ParameterEntry> Map;

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

  template<typename T>
  const T& getParameter(const string& name, T def_value) const;
  //@}

  //! Return true if a parameter with this name exists.
  bool isParameter(const string& name) const;

  bool isParameterSublist(const string& name) const;

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

  if (i == params_.end()) {
    params_[name].setValue(def_value, true);
    i = params_.find(name);
  }

  if ( i != params_.end() )
    return getValue<T>(entry(i));
}

template<typename T>
const T& ParameterList::getParameter(const string& name, T def_value) const
{
  ConstIterator i = params_.find(name);

  if (i == params_.end()) {
    params_[name].setValue(def_value, true);
    i = params_.find(name);
  }

  if ( i != params.end() )
    //return entry(i).getValue(def_value);
    return getValue<T>(entry(i));
}

inline ostream& operator<<(ostream& os, const ParameterList& l)
{
  return l.print(os);
}

} // end of Teuchos namespace

#endif


