// @HEADER
// ***********************************************************************
// 
//                    Teuchos: Common Tools Package
//                 Copyright (2004) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// 
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//  
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//  
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
// Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
// 
// ***********************************************************************
// @HEADER


#ifndef TEUCHOS_PARAMETER_LIST_H
#define TEUCHOS_PARAMETER_LIST_H

/*! \file Teuchos_ParameterList.hpp
    \brief Templated Parameter List class
*/  

#include "Teuchos_ParameterEntry.hpp" // class data element 
#include "Teuchos_TestForException.hpp"
#include "Teuchos_ConfigDefs.hpp"
#include "Teuchos_map.hpp"
#include <typeinfo>

/*! \class Teuchos::ParameterList
    \brief The Teuchos::ParameterList class provides a templated parameter list.
  
    Parameters can be added and retreived with the templated "get" and "set" methods.
    These parameters can be the standard types (double, float, int, ...) and parameter lists, 
    allowing for a hierarchy of parameter lists.  These parameters can also be 
    pointers to vectors or functions.  

    \note <ul>
	  <li> Use static_cast<T>() when the type is ambiguous.  
          <li> Both char* and string map to are stored as strings internally. 
	  </ul>
*/

namespace Teuchos {

class ParameterList {

  //! Parameter container typedef
  typedef Teuchos::map<string, ParameterEntry> Map;

  //! Parameter container const iterator typedef
  typedef Map::const_iterator ConstIterator;

  //! Parameter container iterator typedef
  typedef Map::iterator Iterator;
  
public:
  
  //@{ \name Constructors/Destructor.
  //! Constructor
  ParameterList();
  
  //! Copy Constructor
  ParameterList(const ParameterList& source);
  
  //! Deconstructor
  ~ParameterList();
  //@}
  
  //@{ \name Set Methods 
  
  //! Replace the current parameter list with \c source.
  ParameterList& operator=(const ParameterList& source);
  
  /*! \brief Sets different types of parameters. The type depends on the second entry.  
    
    \note <ul>
    <li> Use static_cast<T>() when the type is ambiguous. 
    <li> Both char* and string map to are stored as strings internally. 
    <li> Sets the parameter as "unused".
    </ul>
  */
  template<typename T>
  void set(const string& name, T value);

  /*! \brief Template specialization for the case when a user sets the parameter with a character
    string in parenthesis.
  */
  void set(const string& name, char* value) 
  { set( name, std::string(value) ); }

  /*! \brief Template specialization for the case when a user sets the parameter with a character
    string in parenthesis.
  */
  void set(const string& name, const char* value) 
  { set( name, std::string(value) ); }

  /*! \brief Template specialization for the case when a user sets the parameter with a ParameterList.
   */
  void set(const string& name, ParameterList value)
  { sublist(name) = value; }
  
  //@}
  
  //@{ \name Get Methods
  
  /*! \brief Retrieves parameter \c name of type \c T from list, if it exists, else the \c def_value is
    used to enter a new parameter into the list.
    
    \note <ul> 
    <li> Use the static_cast<T>() when the type is ambiguous.
    <li> Both char* and string map to are stored as strings internally. 
    <li> Sets the parameter as "used".
    </ul>
  */
  template<typename T>
  T& get(const string& name, T def_value);
  
  /*! \brief Template specialization of get, where the nominal value is a character string in parenthesis.
    Both char* and string are stored as strings and return string values.
  */
  std::string& get(const string& name, char* def_value)
  { return get(name, std::string(def_value)); }

  /*! \brief Template specialization of get, where the nominal value is a character string in parenthesis.
    Both char* and string are stored as strings and return string values.
  */
  std::string& get(const string& name, const char* def_value)
  { return get(name, std::string(def_value)); }
  
  /*! \brief Retrieves parameter \c name of type \c T from a list, an exception is thrown if this parameter doesn't exist.
    \note The syntax for calling this method is:  <tt> list.template get<int>( "Iters" ) </tt>
  */
  template<typename T>
  T& get(const string& name);
  
  /*! \brief Retrieves parameter \c name of type \c T from a constant list, an exception is thrown if this parameter doesn't exist.
    \note The syntax for calling this method is:  <tt> list.template get<int>( "Iters" ) </tt>
  */
  template<typename T>
  const T& get(const string& name) const;  
  
  //@}
  
  //@{ \name Sublist Methods
  /*! \brief Creates an empty sublist and returns a reference to the sublist \c name. If the list already exists, returns reference to that sublist. If the name exists but is not a sublist, an exception is thrown.
   */
  ParameterList& sublist(const string& name);
  
  /*! \brief Return a const reference to an existing sublist \c name.  If the list does not already exist or the name exists but is not a sublist, an exception is thrown.
   */
  const ParameterList& sublist(const string& name) const;
  //@}
  
  //@{ \name Attribute Methods
  /*! \brief Query the existence of a parameter.
    \return "true" if a parameter with this \c name exists, else "false".
  */
  bool isParameter(const string& name) const;
  
  /*! \brief Query the existence of a parameter and whether it is a parameter list.
    \return "true" if a parameter with this \c name exists and is itself a parameter list, else "false".
  */
  bool isSublist(const string& name) const;
  
  /*! \brief Query the existence and type of a parameter.
    \return "true" is a parameter with this \c name exists and is of type \c T, else "false".
    \note The syntax for calling this method is:  <tt> list.template isType<int>( "Iters" ) </tt>
  */
  template<typename T>
  bool isType(const string& name) const;

#ifndef DOXYGEN_SHOULD_SKIP_THIS  
  /*! \brief Query the existence and type of a parameter.
    \return "true" is a parameter with this \c name exists and is of type \c T, else "false".
    \note <b>It is not recommended that this method be used directly!</b>  
    Please use either the helper function <b>isParameterType</b> or non-nominal <b>isType</b> method. 
  */
  template<typename T>
  bool isType(const string& name, T* ptr) const;
#endif
  //@}
  
  //@{ \name I/O Methods
  //! Printing method for parameter lists.  Indenting is used to indicate parameter list hierarchies. 
  ostream& print(ostream& os, int indent = 0) const;
  
  //! Print out unused parameters in the ParameterList.
  void unused(ostream& os) const;
  //@}
  
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
};
  
  
  template<typename T>
  void ParameterList::set(const string& name, T value)
  {
    params_[name].setValue(value);
  }
  
  template<typename T>
  T& ParameterList::get(const string& name, T def_value)
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
  T& ParameterList::get(const string& name) 
  {
    ConstIterator i = params_.find(name);
    
    // This parameter was not found or is wrong type, throw an exception
    TEST_FOR_EXCEPTION( i == params_.end(), std::runtime_error,
			"get ( " << name << " ) failed -- parameter does not exist! " );
    TEST_FOR_EXCEPTION( !isType( name, (T*)NULL ), std::runtime_error,
			"get ( " << name << " ) failed -- parameter is wrong type! " );
    
    // Return the value.
    return getValue<T>(entry(i));
  }
  
  template<typename T>
  const T& ParameterList::get(const string& name) const
  {
    ConstIterator i = params_.find(name);
    
    // This parameter was not found, throw and exception
    TEST_FOR_EXCEPTION( i == params_.end(), std::runtime_error,
			"get ( " << name << " ) failed -- parameter does not exist! " );
    TEST_FOR_EXCEPTION( !isType( name, (T*)NULL ), std::runtime_error,
			"get ( " << name << " ) failed -- parameter is wrong type! " );
    
    // Return the default value for this type
    return getValue<T>(entry(i));
  }
  
#ifndef DOXYGEN_SHOULD_SKIP_THIS
  template<typename T>
  bool ParameterList::isType(const string& name, T* ptr) const
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
#endif
  
  template<typename T>
  bool ParameterList::isType(const string& name) const
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
  
  
  /*! \relates ParameterList
    \brief A templated helper function for getting a parameter from a non-const list.
    This helper function prevents the need for giving a nominal value of the specific template type.
    
    \note The syntax for calling this function is:  <tt> getParameter<int>( list, "Iters" ) </tt>
  */
  template<typename T>
  T& getParameter( ParameterList& l, const string& name )
  {
    T temp_var;
    // CAREFUL:  We need to be sure the parameter exists before using the temp variable to get the parameter.
    // This parameter was not found or is wrong type, throw an exception
    TEST_FOR_EXCEPTION( !l.isParameter( name ), std::runtime_error,
			"getParameter ( " << name << " ) failed -- parameter does not exist! " );
    TEST_FOR_EXCEPTION( !l.isType( name, &temp_var ), std::runtime_error,
			"getParameter ( " << name << " ) failed -- parameter is wrong type! " );
    //
    // This parameter exists and is of the right type, so we can retrieve it safely.    
    //
    return l.get( name, temp_var );
  }
  
  /*! \relates ParameterList
    \brief A templated helper function for getting a parameter from a const list.
    This helper function prevents the need for giving a nominal value of the specific template type.
    
    \note The syntax for calling this function is:  <tt> getParameter<int>( list, "Iters" ) </tt>    
  */
  template<typename T>
  const T& getParameter( const ParameterList& l, const string& name )
  {
    T temp_var;
    // CAREFUL:  We need to be sure the parameter exists before using the temp variable to get the parameter.
    // This parameter was not found or is wrong type, throw an exception
    TEST_FOR_EXCEPTION( !l.isParameter( name ), std::runtime_error,
			"getParameter ( " << name << " ) failed -- parameter does not exist! " );
    TEST_FOR_EXCEPTION( !l.isType( name, &temp_var ), std::runtime_error,
			"getParameter ( " << name << " ) failed -- parameter is wrong type! " );
    //
    // This parameter exists and is of the right type, so we can retrieve it safely.    
    //
    return l.get( name, temp_var );
  }
  
  /*! \relates ParameterList
    \brief A templated helper function for determining the type of a parameter entry for a non-const list.  
    This helper function avoids the need for giving a nominal value of the specific template type.
    
    \note The syntax for calling this function is:  <tt> isParameterType<int>( list, "Iters" ) </tt>
  */
  template<typename T>
  bool isParameterType( ParameterList& l, const string& name )
  {
    return l.isType( name, (T*)NULL );
  }
  
  /*! \relates ParameterList
    \brief A templated helper function for determining the type of a parameter entry for a const list.  
    This helper function avoids the need for giving a nominal value of the specific template type.
    
    \note The syntax for calling this function is:  <tt> isParameterType<int>( list, "Iters" ) </tt>
  */
  template<typename T>
  bool isParameterType( const ParameterList& l, const string& name )
  {
    return l.isType( name, (T*)NULL );
  }
  
  /*! \relates ParameterList
    \brief Output stream operator for handling the printing of the parameter list.
  */
  inline ostream& operator<<(ostream& os, const ParameterList& l)
  {
    return l.print(os);
  }
  
} // end of Teuchos namespace

#endif


