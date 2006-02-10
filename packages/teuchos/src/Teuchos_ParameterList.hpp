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

/** \example ParameterList/cxx_main.cpp
    This is an example of how to use the Teuchos::ParameterList class.
*/

namespace Teuchos {

/** \brief .
 * \relates ParameterList
 */
enum EValidateUsed {
  VALIDATE_USED_ENABLED   /*< Validate that parameters in <tt>*this</tt> list
                              set using the default value are present in
                              the validation list */
  ,VALIDATE_USED_DISABLED /*< Do not validate that parameters in <tt>*this</tt> list
                              set using the default value are present in
                              the validation list */
};

/** \brief .
 * \relates ParameterList
 */
enum EValidateDefaults {
  VALIDATE_DEFAULTS_ENABLED   /*< Validate that parameters in <tt>*this</tt> list
                                  set using the default value are present in
                                   the validation list */
  ,VALIDATE_DEFAULTS_DISABLED /*< Do not validate that parameters in <tt>*this</tt> list
                                  set using the default value are present in
                                  the validation list */
};

namespace Exceptions {

/** \brief .
 * \relates ParameterList
 */
class InvalidParameter : public std::logic_error
{public: InvalidParameter(const std::string& what_arg) : std::logic_error(what_arg) {}};

} // namespace Exceptions

/*! \brief Templated parameter list.
  
    Parameters can be added and retreived with the templated "get" and "set"
    functions.  These parameters can any data type which uses value sementics
    (e.g. double, float, int, *double, *float, *int, ...) which includes other
    parameter lists, allowing for a hierarchy of parameter lists.  These
    parameters can also be pointers to vectors or functions.

    \note <ul>
	  <li> Use static_cast<T>() when the type is ambiguous.  
          <li> Both char* and string map to are stored as strings internally. 
	  </ul>
*/
class ParameterList {

  //! Parameter container typedef
  typedef Teuchos::map<string, ParameterEntry> Map;

  //! Parameter container iterator typedef
  typedef Map::iterator Iterator;
  
public:

  /** \name Public types */
  //@{

  //! Parameter container const iterator typedef
  typedef Map::const_iterator ConstIterator;

  //@}
  
  /** \name Constructors/Destructor. */
  //@{

  //! Constructor
  ParameterList();
  
  //! Copy Constructor
  ParameterList(const ParameterList& source);
  
  //! Deconstructor
  ~ParameterList();

  //@}
  
  /** \name Set Functions */
  //@{
  
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
  void set(const string& name, char* value);

  /*! \brief Template specialization for the case when a user sets the parameter with a character
    string in parenthesis.
  */
  void set(const string& name, const char* value);

  /*! \brief Template specialization for the case when a user sets the parameter with a ParameterList.
   */
  void set(const string& name, ParameterList value);

  /*! \brief Set a parameter directly as a ParameterEntry. 
   * \note This is required to preserve the isDefault value when reading back
   * from XML. KL 7 August 2004 
   */
  void setEntry(const string& name, const ParameterEntry& entry);

  //@}
  
  /** \name Get Functions */
  //@{
  
  /*! \brief Retrieves parameter \c name of type \c T from list, if it exists, else the \c def_value is
    used to enter a new parameter into the list.
    
    \note <ul> 
    <li> Use the static_cast<T>() when the type is ambiguous.
    <li> Both char* and string map to are stored as strings internally. 
    <li> Sets the parameter as "used".
    <li> Exception is thrown if \c name exists, but is not of type \c T.
    </ul>
  */
  template<typename T>
  T& get(const string& name, T def_value);
  
  /*! \brief Template specialization of get, where the nominal value is a character string in parenthesis.
    Both char* and string are stored as strings and return string values.
  */
  std::string& get(const string& name, char* def_value);

  /*! \brief Template specialization of get, where the nominal value is a character string in parenthesis.
    Both char* and string are stored as strings and return string values.
  */
  std::string& get(const string& name, const char* def_value);
  
  /*! \brief Retrieves parameter \c name of type \c T from a list, an
    exception is thrown if this parameter doesn't exist or is the wrong type.
    \note The syntax for calling this method is: <tt> list.template get<int>(
    "Iters" ) </tt>
  */
  template<typename T>
  T& get(const string& name);
  
  /*! \brief Retrieves parameter \c name of type \c T from a constant list, an
    exception is thrown if this parameter doesn't exist or is the wrong type.
    \note The syntax for calling this method is: <tt> list.template get<int>(
    "Iters" ) </tt>
  */
  template<typename T>
  const T& get(const string& name) const;  
  
  /*! \brief Retrieves the pointer for parameter \c name of type \c T from a
    list.  A null pointer is returned if this parameter doesn't exist or is
    the wrong type.  \note The syntax for calling this method is: <tt>
    list.template getPtr<int>( "Iters" ) </tt>
  */
  template<typename T>
  T* getPtr(const string& name);
  
  /*! \brief Retrieves the pointer for parameter \c name of type \c T from a
    constant list.  A null pointer is returned if this parameter doesn't exist
    or is the wrong type.  \note The syntax for calling this method is: <tt>
    list.template getPtr<int>( "Iters" ) </tt>
  */
  template<typename T>
  const T* getPtr(const string& name) const;  
  
  /*! \brief Retrieves the pointer for an entry with the name <tt>name</tt> if
   *  it exists. */
  const ParameterEntry* getEntryPtr(const string& name) const;  

  //@}
  
  /** \name Sublist Functions */
  //@{

  /*! \brief Creates an empty sublist and returns a reference to the sublist
   *  \c name. If the list already exists, returns reference to that
   *  sublist. If the name exists but is not a sublist, an exception is
   *  thrown.
   */
  ParameterList& sublist(const string& name);
  
  /*! \brief Return a const reference to an existing sublist \c name.  If the
   *  list does not already exist or the name exists but is not a sublist, an
   *  exception is thrown.
   */
  const ParameterList& sublist(const string& name) const;

  //@}
  
  /** \name Attribute Functions */
  //@{

  /*! \brief Query the existence of a parameter.  \return "true" if a
    parameter with this \c name exists, else "false".
  */
  bool isParameter(const string& name) const;
  
  /*! \brief Query the existence of a parameter and whether it is a parameter
    list.  \return "true" if a parameter with this \c name exists and is
    itself a parameter list, else "false".
  */
  bool isSublist(const string& name) const;
  
  /*! \brief Query the existence and type of a parameter.  \return "true" is a
    parameter with this \c name exists and is of type \c T, else "false".
    \note The syntax for calling this method is: <tt> list.template
    isType<int>( "Iters" ) </tt>
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
  
  /** \name I/O Functions */
  //@{

  /*! \brief Printing method for parameter lists.  Indenting is used to indicate
    parameter list hierarchies. */
  ostream& print(ostream& os, int indent = 0, bool showTypes = false) const;
  
  //! Print out unused parameters in the ParameterList.
  void unused(ostream& os) const;

  //@}

  /** \name Read-only access to the iterator */
  //@{

  //! An iterator pointing to the first entry
  ConstIterator begin() const ;

  //! An iterator pointing beyond the last entry
  ConstIterator end() const ;

  //! Access to ParameterEntry (i.e., returns i->second)
  const ParameterEntry& entry(ConstIterator i) const;
  
  //! Access to name (i.e., returns i->first)
  const string& name(ConstIterator i) const;

  //@}

  /** \name Validation Functions */
  //@{

  /** \brief Validate the parameters is this list given valid selections in
   * the input list.
   *
   * \param  paramListName
   *              [in] The name of <tt>*this</tt> parameter list.  This is used to generate
   *              the exception message in case of a failure.
   * \param  validParamList
   *              [in] This is the list that the parameters and sublist in <tt>*this</tt>
   *              are compared against.
   * \param  depth
   *              [in] Determines the number of levels of depth that the validation will
   *              recurse into.  A value of <tt>dpeth=0</tt> means that only the top level
   *              parameters and sublists will be checked.  Default: <tt>depth = large number</tt>.
   * \param  validateUsed
   *              [in] Determines if parameters that have been used are checked
   *              against those in <tt>validParamList</tt>.  Default:
   *              <tt>validateDefaults = VALIDATE_DEFAULTS_ENABLED</tt>.
   * \param  validateDefaults
   *              [in] Determines if parameters set at their default values using <tt>get(name,defaultVal)</tt>
   *              are checked against those in <tt>validParamList</tt>.  Default:
   *              <tt>validateDefaults = VALIDATE_DEFAULTS_ENABLED</tt>.
   *
   * If a parameter in <tt>*this</tt> is not found in
   * <tt>validParamList</tt> then an exception of type
   * <tt>Exceptions::InvalidParameter</tt> will be thrown which contains an
   * excellent error message returned by <tt>excpt.what()</tt>.
   *
   * A depth-first search is performed to validate parameters since this is
   * the most effiicent way to perform the validation.  This could be changed
   * to a breath-first search but it would take more work.
   */
  void validateParameters(
    const std::string          &paramListName   
    ,const ParameterList       &validParamList
    ,const int                 depth            = 1000
    ,const EValidateUsed       validateUsed     = VALIDATE_USED_ENABLED
    ,const EValidateDefaults   validateDefaults = VALIDATE_DEFAULTS_ENABLED
    ) const;

  //@}
  
private: // Functions
  
  //! Access to ParameterEntry (i.e., returns i->second)
  ParameterEntry& entry(Iterator i);
  
private: // Data members
  
  //! Parameter list
  Map params_;
};

// /////////////////////////////////////////////////////
// Inline and Template Function Definitions

// Set functions
 
template<typename T>
inline
void ParameterList::set(const string& name, T value)
{
  params_[name].setValue(value);
}

inline
void ParameterList::set(const string& name, char* value) 
{ set( name, std::string(value) ); }

inline
void ParameterList::set(const string& name, const char* value) 
{ set( name, std::string(value) ); }

inline
void ParameterList::set(const string& name, ParameterList value)
{ sublist(name) = value; }

inline
void ParameterList::setEntry(const string& name, const ParameterEntry& entry)
{params_[name] = entry;}

// Get functions

template<typename T>
T& ParameterList::get(const string& name, T def_value)
{
  ConstIterator i = params_.find(name);
    
  // The parameter was not found, add it to the list
  if (i == params_.end()) {
    params_[name].setValue(def_value, true);
    i = params_.find(name);
  } else {
    // The parameter was found, make sure it is the same type as T.
    TEST_FOR_EXCEPTION( entry(i).getAny().type() != typeid(T), std::runtime_error,
                        "get ( " << name << ", T def_value ) failed -- parameter is wrong type! " );
  }

  // Return the value of the parameter
  return getValue<T>(entry(i));
}

inline
std::string& ParameterList::get(const string& name, char* def_value)
{ return get(name, std::string(def_value)); }

inline
std::string& ParameterList::get(const string& name, const char* def_value)
{ return get(name, std::string(def_value)); }

template<typename T>
T& ParameterList::get(const string& name) 
{
  ConstIterator i = params_.find(name);
  TEST_FOR_EXCEPTION( i == params_.end(), std::runtime_error,
                      "get ( " << name << " ) failed -- parameter does not exist! " );
  TEST_FOR_EXCEPTION( entry(i).getAny().type() != typeid(T), std::runtime_error,
                      "get ( " << name << " ) failed -- parameter is wrong type! " );
  return getValue<T>(entry(i));
}
  
template<typename T>
const T& ParameterList::get(const string& name) const
{
  ConstIterator i = params_.find(name);
  TEST_FOR_EXCEPTION( i == params_.end(), std::runtime_error,
                      "get ( " << name << " ) failed -- parameter does not exist! " );
  TEST_FOR_EXCEPTION( entry(i).getAny().type() != typeid(T), std::runtime_error,
                      "get ( " << name << " ) failed -- parameter is wrong type! " );
  return getValue<T>(entry(i));
}

template<typename T>
inline
T* ParameterList::getPtr(const string& name) 
{
  ConstIterator i = params_.find(name);
  if ( i == params_.end() || entry(i).getAny().type() != typeid(T) )
    return NULL;
  return &getValue<T>(entry(i));
}
  
template<typename T>
inline
const T* ParameterList::getPtr(const string& name) const
{
  ConstIterator i = params_.find(name);
  if ( i == params_.end() || entry(i).getAny().type() != typeid(T) )
    return NULL;
  return &getValue<T>(entry(i));
}

inline
const ParameterEntry*
ParameterList::getEntryPtr(const string& name) const
{
  ConstIterator i = params_.find(name);
  if ( i == params_.end() )
    return NULL;
  return &entry(i);
}

// Attribute Functions
  
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
  catch( std::exception& e ) {
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
  catch( std::exception& e ) {
    return false;
  }
  // If no exception was thrown, we should be OK.
  return true;
}

// //////////////////////////////////////
// Helper functions
  
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


