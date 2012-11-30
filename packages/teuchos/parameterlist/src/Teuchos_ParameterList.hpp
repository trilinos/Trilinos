// @HEADER
// ***********************************************************************
//
//                    Teuchos: Common Tools Package
//                 Copyright (2004) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
//
// ***********************************************************************
// @HEADER


#ifndef TEUCHOS_PARAMETER_LIST_H
#define TEUCHOS_PARAMETER_LIST_H

/*! \file Teuchos_ParameterList.hpp
    \brief Templated Parameter List class
*/  

#include "Teuchos_ParameterListExceptions.hpp"
#include "Teuchos_ParameterEntry.hpp"
#include "Teuchos_StringIndexedOrderedValueObjectContainer.hpp"
#include "Teuchos_Assert.hpp"
#include "Teuchos_RCP.hpp"
#include "Teuchos_Array.hpp"
#include "Teuchos_map.hpp"


/** \example ParameterList/cxx_main.cpp
    This is an example of how to use the Teuchos::ParameterList class.
*/

namespace Teuchos {

/** \brief Validation used enum.
 * \relates ParameterList
 */
enum EValidateUsed {
  VALIDATE_USED_ENABLED   /*< \brief Validate that parameters in <tt>*this</tt> list
                              set using the default value are present in
                              the validation list */
  ,VALIDATE_USED_DISABLED /*< \brief Do not validate that parameters in <tt>*this</tt> list
                              set using the default value are present in
                              the validation list */
};

/** \brief Validation defaults enum.
 * \relates ParameterList
 */
enum EValidateDefaults {
  VALIDATE_DEFAULTS_ENABLED   /*< \brief Validate that parameters in <tt>*this</tt> list
                                  set using the default value are present in
                                   the validation list */
  ,VALIDATE_DEFAULTS_DISABLED /*< \brief Do not validate that parameters in <tt>*this</tt> list
                                  set using the default value are present in
                                  the validation list */
};

/*! \brief Templated parameter list.
  
    Parameters can be added and retreived with the templated "get" and "set"
    functions.  These parameters can any data type which uses value sementics
    (e.g. double, float, int, *double, *float, *int, ...) which includes other
    parameter lists, allowing for a hierarchy of parameter lists.  These
    parameters can also be pointers to vectors or functions.

    \note <ul>
	  <li> Use static_cast<T>() when the type is ambiguous.  
          <li> Both char* and std::string std::map to are stored as strings internally. 
	  </ul>
*/
class TEUCHOSPARAMETERLIST_LIB_DLL_EXPORT ParameterList {

  //! Internal data-structure
  typedef StringIndexedOrderedValueObjectContainer<ParameterEntry> params_t;

  //! Parameter container iterator typedef
  typedef params_t::Iterator Iterator;
  
public:

  //! @name Public types 
  //@{

  //! Parameter container const iterator typedef
  typedef params_t::ConstIterator ConstIterator;

  /** \brief Utility class for setting and passing in print options. */
  class PrintOptions {
  public:
    PrintOptions() : indent_(0), showTypes_(false), showFlags_(false), showDoc_(false) {}
    PrintOptions& indent(int _indent)        { indent_ = _indent; return *this; }
    PrintOptions& showTypes(bool _showTypes) { showTypes_ = _showTypes; return *this; }
    PrintOptions& showFlags(bool _showFlags) { showFlags_ = _showFlags; return *this; }
    PrintOptions& showDoc(bool _showDoc)     { showDoc_ = _showDoc; return *this; }
    PrintOptions& incrIndent(int indents)    { indent_ += indents; return *this; }
    int indent() const { return indent_; }
    bool showTypes() const { return showTypes_; }
    bool showFlags() const { return showFlags_; }
    bool showDoc() const { return showDoc_; }
    PrintOptions copy() const { return PrintOptions(*this); }
  private:
    int    indent_;
    bool   showTypes_;
    bool   showFlags_;
    bool   showDoc_;
  };

  //@}
  
  //! @name Constructors/Destructor/Info. 
  //@{

  //! Constructor
  ParameterList();

  //! Constructor
  ParameterList(const std::string &name);
  
  //! Copy Constructor
  ParameterList(const ParameterList& source);
  
  //! Deconstructor
  virtual ~ParameterList();

  //! Get the number of stored parameters.
  Ordinal numParams() const;

  //@}
  
  //! @name Set Functions 
  //@{

  /** \brief Set the name of <tt>*this</tt> list.
   */
  ParameterList& setName( const std::string &name );
  
  /** Replace the current parameter list with \c source.
   * \note This also replaces the name returned by <tt>this->name()</tt>
   */
  ParameterList& operator=(const ParameterList& source);
  
  /** Set the parameters in <tt>source</tt>.
   *
   * Note, this function will set the parameters and sublists from
   * <tt>source</tt> into <tt>*this</tt> but will not result in parameters
   * being removed from <tt>*this</tt>.  Parameters in <tt>*this</tt> with the
   * same names as those in <tt>source</tt> will be overwritten.
   */
  ParameterList& setParameters(const ParameterList& source);
  
  /** Set the parameters in <tt>source</tt> that are not already set in
   * <tt>*this</tt>.
   *
   * Note, this function will set the parameters and sublists from
   * <tt>source</tt> into <tt>*this</tt> but will not result in parameters
   * being removed from <tt>*this</tt> or in parameters already set in
   * <tt>*this</tt> being overrided.  Parameters in <tt>*this</tt> with the
   * same names as those in <tt>source</tt> will not be overwritten.
   */
  ParameterList& setParametersNotAlreadySet(const ParameterList& source);

  /** Disallow recusive validation when this sublist is used in a valid
   * parameter list.
   *
   * This function should be called when setting a sublist in a valid
   * parameter list which is broken off to be passed to another object.
   * The other object should validate its own list.
   */
  ParameterList& disableRecursiveValidation();
  
  /*! \brief Sets different types of parameters. The type depends on the second entry.  
    
    \note <ul>
    <li> Use static_cast<T>() when the type is ambiguous. 
    <li> Both char* and std::string std::map to are stored as strings internally. 
    <li> Sets the parameter as "unused".
    </ul>
  */
  template<typename T>
  ParameterList& set(
    std::string const& name, T const& value, std::string const& docString = "",
    RCP<const ParameterEntryValidator> const& validator = null
    );

  /*! \brief Specialization for the case when a user sets the parameter with a character
    std::string in parenthesis.
  */
  ParameterList& set(
    std::string const& name, char value[], std::string const& docString = "",
    RCP<const ParameterEntryValidator> const& validator = null
    );

  /*! \brief Specialization for the case when a user sets the parameter with a character
    std::string in parenthesis.
  */
  ParameterList& set(
    std::string const& name, const char value[], std::string const& docString = "",
    RCP<const ParameterEntryValidator> const& validator = null
    );

  /*! \brief Template specialization for the case when a user sets the parameter with a ParameterList.
   */
  ParameterList& set(
    std::string const& name, ParameterList const& value, std::string const& docString = ""
    );

  /*! \brief Set a parameter directly as a ParameterEntry. 
   * \note This is required to preserve the isDefault value when reading back
   * from XML. KL 7 August 2004 
   */
  ParameterList& setEntry(const std::string& name, const ParameterEntry& entry);

  //@}
  
  //! @name Get Functions 
  //@{
  
  /*! \brief Retrieves parameter \c name of type \c T from list, if it exists, else the \c def_value is
    used to enter a new parameter into the list.
    
    \note <ul> 
    <li> Use the static_cast<T>() when the type is ambiguous.
    <li> Both char* and std::string std::map to are stored as strings internally. 
    <li> Sets the parameter as "used".
    <li> Exception is thrown if \c name exists, but is not of type \c T.
    </ul>
  */
  template<typename T>
  T& get(const std::string& name, T def_value);

  /*! \brief Specialization of get, where the nominal value is a character std::string in parenthesis.
    Both char* and std::string are stored as strings and return std::string values.
  */
  std::string& get(const std::string& name, char def_value[]);
  
  /*! \brief Specialization of get, where the nominal value is a character std::string in parenthesis.
    Both char* and std::string are stored as strings and return std::string values.
  */
  std::string& get(const std::string& name, const char def_value[]);
  
  /*! \brief Retrieves parameter \c name of type \c T from a list, an
    <tt>Exceptions::InvalidParameter</tt> std::exception is thrown if this
    parameter doesn't exist (<tt>Exceptions::InvalidParameterName</tt>) or is
    the wrong type (<tt>Exceptions::InvalidParameterType</tt>).  \note The
    syntax for calling this method is: <tt> list.template get<int>( "Iters" )
    </tt>
  */
  template<typename T>
  T& get(const std::string& name);
  
  /*! \brief Retrieves parameter \c name of type \c T from a constant list, an
    <tt>Exceptions::InvalidParameter</tt> std::exception is thrown if this
    parameter doesn't exist (<tt>Exceptions::InvalidParameterName</tt>) or is
    the wrong type (<tt>Exceptions::InvalidParameterType</tt>).  \note The
    syntax for calling this method is: <tt> list.template get<int>( "Iters" )
    </tt>
  */
  template<typename T>
  const T& get(const std::string& name) const;  
  
  /*! \brief Retrieves the pointer for parameter \c name of type \c T from a
    list.  A null pointer is returned if this parameter doesn't exist or is
    the wrong type.  \note The syntax for calling this method is: <tt>
    list.template getPtr<int>( "Iters" ) </tt>
  */
  template<typename T>
  inline
  T* getPtr(const std::string& name);
  
  /*! \brief Retrieves the pointer for parameter \c name of type \c T from a
    constant list.  A null pointer is returned if this parameter doesn't exist
    or is the wrong type.  \note The syntax for calling this method is: <tt>
    list.template getPtr<int>( "Iters" ) </tt>
  */
  template<typename T>
  inline
  const T* getPtr(const std::string& name) const;  

  // ToDo: Add getSafePtr() functions to return Ptr<T> instead of raw T*
  
  /*! \brief Retrieves an entry with the name <tt>name</tt>.
   *
   * Throws <tt>Exceptions::InvalidParameterName</tt> if this parameter does
   * not exist.
   */
  ParameterEntry& getEntry(const std::string& name);  
  
  /*! \brief Retrieves a const entry with the name <tt>name</tt>.
   *
   * Throws <tt>Exceptions::InvalidParameterName</tt> if this parameter does
   * not exist.
   */
  inline
  const ParameterEntry& getEntry(const std::string& name) const;  
  
  /*! \brief Retrieves the pointer for an entry with the name <tt>name</tt> if
   *  it exists. */
  inline
  ParameterEntry* getEntryPtr(const std::string& name);  

  // ToDo: Add function called getEntrySafePtr() to return Ptr<> as the main
  // implementation and deprecate getEntryPtr()
  
  /*! \brief Retrieves the pointer for a constant entry with the name <tt>name</tt> if
   *  it exists. */
  inline
  const ParameterEntry* getEntryPtr(const std::string& name) const;  

  /*! \brief Retrieves the RCP for an entry with the name <tt>name</tt> if
   *  it exists. */
  inline RCP<ParameterEntry> getEntryRCP(const std::string& name);  
  
  /*! \brief Retrieves the RCP for a constant entry with the name <tt>name</tt> if
   *  it exists. */
  inline RCP<const ParameterEntry> getEntryRCP(const std::string& name) const;

  //@}

  //! @name Parameter removal functions
  //@{
 
  /** \brief Remove a parameter (does not depend on the type of the
   * parameter).
   *
   * \param name [in] The name of the parameter to remove
   *
   * \param throwIfNotExists [in] If <tt>true</tt> then if the parameter with
   * the name <tt>name</tt> does not exist then a std::exception will be
   * thrown!
   *
   * \returns Returns <tt>true</tt> if the parameter was removed, and
   * <tt>false</tt> if the parameter was not removed (<tt>false</tt> return
   * value possible only if <tt>throwIfExists==false</tt>).
   */
  bool remove(
    std::string const& name, bool throwIfNotExists = true
    );

  //@}
  
  //! @name Sublist Functions 
  //@{

  /*! \brief Creates an empty sublist and returns a reference to the sublist
   *  \c name. If the list already exists, returns reference to that
   *  sublist. If the name exists but is not a sublist, an std::exception is
   *  thrown.
   */
  ParameterList& sublist(
    const std::string& name, bool mustAlreadyExist = false
    ,const std::string& docString = ""
    );
  
  /*! \brief Return a const reference to an existing sublist \c name.  If the
   *  list does not already exist or the name exists but is not a sublist, an
   *  std::exception is thrown.
   */
  const ParameterList& sublist(const std::string& name) const;

  //@}
  
  //! @name Attribute Functions 
  //@{

  /*! \brief Query the name of this parameter list. */
  const std::string& name() const;

  /*! \brief Query the existence of a parameter.  \return "true" if a
      parameter with this \c name exists, else "false".  Warning, this
      function should almost never be used!  Instead, consider using
      getEntryPtr() instead. \todo Depreicate this function and instead use
      index lookup.
  */
  bool isParameter(const std::string& name) const;
  
  /*! \brief Query the existence of a parameter and whether it is a parameter
      list.  \return "true" if a parameter with this \c name exists and is
      itself a parameter list, else "false".  Warning, this function should
      almost never be used!  Instead, consider using getEntryPtr() instead.
  */
  bool isSublist(const std::string& name) const;
  
  /*! \brief Query the existence and type of a parameter.  \return "true" is a
    parameter with this \c name exists and is of type \c T, else "false".
    \note The syntax for calling this method is: <tt> list.template
    isType<int>( "Iters" ) </tt>.  Warning, this function should almost never
    be used!  Instead, consider using getEntryPtr() instead.
  */
  template<typename T>
  bool isType(const std::string& name) const;

#ifndef DOXYGEN_SHOULD_SKIP_THIS  
  /*! \brief Query the existence and type of a parameter.
   *
   * \return "true" is a parameter with this \c name exists and is of type \c
   * T, else "false".
   *
   * \note <b>It is not recommended that this method be used directly!</b>
   *
   * Please use either the helper function <b>isParameterType</b> or
   * non-nominal <b>isType</b> method.
  */
  template<typename T>
  bool isType(const std::string& name, T* ptr) const;
#endif

  //@}
  
  //! @name I/O Functions 
  //@{

  /*! \brief Print function to use in debugging in a debugger.
   *
   * Prints to *VerboseObjectBase::getDefaultOStream() so it will print well
   * in parallel.
   */
  void print() const;

  /*! \brief Printing method for parameter lists which takes an print options
   *  object.*/
  std::ostream& print(std::ostream& os, const PrintOptions &printOptions) const;

  /*! \brief Printing method for parameter lists.  Indenting is used to indicate
    parameter list hierarchies. */
  std::ostream& print(std::ostream& os, int indent = 0, bool showTypes = false, bool showFlags = true ) const;
  
  //! Print out unused parameters in the ParameterList.
  void unused(std::ostream& os) const;

  //! Create a single formated std::string of all of the zero-level parameters in this list
  std::string currentParametersString() const;

  //@}

  //! @name Read-only access to the iterator 
  //@{

  //! An iterator pointing to the first entry
  inline ConstIterator begin() const;

  //! An iterator pointing beyond the last entry
  inline ConstIterator end() const;
  
  //! Access to name (i.e., returns i->first)
  inline const std::string& name(ConstIterator i) const;

  //! Access to ParameterEntry (i.e., returns i->second)
  inline const ParameterEntry& entry(ConstIterator i) const;

  //@}

  //! @name Validation Functions 
  //@{

  /** \brief Validate the parameters in this list given valid selections in
   * the input list.
   *
   * \param validParamList [in] This is the list that the parameters and
   * sublist in <tt>*this</tt> are compared against.
   *
   * \param depth [in] Determines the number of levels of depth that the
   * validation will recurse into.  A value of <tt>depth=0</tt> means that
   * only the top level parameters and sublists will be checked.  Default:
   * <tt>depth = large number</tt>.
   *
   * \param validateUsed [in] Determines if parameters that have been used are
   * checked against those in <tt>validParamList</tt>.  Default:
   * <tt>validateDefaults = VALIDATE_DEFAULTS_ENABLED</tt>.
   *
   * \param validateDefaults [in] Determines if parameters set at their
   * default values using <tt>get(name,defaultVal)</tt> are checked against
   * those in <tt>validParamList</tt>.  Default: <tt>validateDefaults =
   * VALIDATE_DEFAULTS_ENABLED</tt>.
   *
   * If a parameter in <tt>*this</tt> is not found in <tt>validParamList</tt>
   * then an <tt>std::exception</tt> of type
   * <tt>Exceptions::InvalidParameterName</tt> will be thrown which will
   * contain an excellent error message returned by <tt>excpt.what()</tt>.  If
   * the parameter exists but has the wrong type, then an std::exception type
   * <tt>Exceptions::InvalidParameterType</tt> will be thrown.  If the
   * parameter exists and has the right type, but the value is not valid then
   * an std::exception type <tt>Exceptions::InvalidParameterValue</tt> will be
   * thrown.
   *
   * Recursive validation stops when:<ul>
   *
   * <li>The maxinum <tt>depth</tt> is reached
   *
   * <li>A sublist note in <tt>validParamList</tt> has been marked with the
   * <tt>disableRecursiveValidation()</tt> function, or
   *
   * <li>There are not more parameters or sublists left in <tt>*this</tt>
   *
   * </ul>
   *
   * A breath-first search is performed to validate all of the parameters in
   * one sublist before moving into nested subslist.
   */
  void validateParameters(
    ParameterList const& validParamList,
    int const depth = 1000,
    EValidateUsed const validateUsed = VALIDATE_USED_ENABLED,
    EValidateDefaults const validateDefaults = VALIDATE_DEFAULTS_ENABLED
    ) const;

  /** \brief Validate the parameters in this list given valid selections in
   * the input list and set defaults for those not set.
   *
   * \param validParamList [in] This is the list that the parameters and
   * sublist in <tt>*this</tt> are compared against.
   *
   * \param depth [in] Determines the number of levels of depth that the
   * validation will recurse into.  A value of <tt>depth=0</tt> means that
   * only the top level parameters and sublists will be checked.  Default:
   * <tt>depth = large number</tt>.
   *
   * If a parameter in <tt>*this</tt> is not found in <tt>validParamList</tt>
   * then an std::exception of type <tt>Exceptions::InvalidParameterName</tt> will
   * be thrown which will contain an excellent error message returned by
   * <tt>excpt.what()</tt>.  If the parameter exists but has the wrong type,
   * then an std::exception type <tt>Exceptions::InvalidParameterType</tt> will be
   * thrown.  If the parameter exists and has the right type, but the value is
   * not valid then an std::exception type
   * <tt>Exceptions::InvalidParameterValue</tt> will be thrown.  If a
   * parameter in <tt>validParamList</tt> does not exist in <tt>*this</tt>,
   * then it will be set at its default value as determined by
   * <tt>validParamList</tt>.
   *
   * Recursive validation stops when:<ul>
   *
   * <li>The maxinum <tt>depth</tt> is reached
   *
   * <li>A sublist note in <tt>validParamList</tt> has been marked with the
   * <tt>disableRecursiveValidation()</tt> function, or
   *
   * <li>There are not more parameters or sublists left in <tt>*this</tt>
   *
   * </ul>
   *
   * A breath-first search is performed to validate all of the parameters in
   * one sublist before moving into nested subslist.
   */
  void validateParametersAndSetDefaults(
    ParameterList const& validParamList,
    int const depth = 1000
    );

  //@}
  
private: // Functions

  //! An iterator pointing to the first entry
  inline Iterator nonconstBegin();
  //! An iterator pointing beyond the last entry
  inline Iterator nonconstEnd();
  //! Access to ParameterEntry (i.e., returns i->second)
  inline ParameterEntry& nonconstEntry(Iterator i);
  //! Validate that a parameter exists
  void validateEntryExists(const std::string &funcName, const std::string &name,
    const ParameterEntry *entry) const;
  // ToDo: Change above function to take Ptr<> instead of raw pointer.
  //! Validate that a type is the same
  template<typename T>
  void validateEntryType(const std::string &funcName, const std::string &name,
    const ParameterEntry &entry ) const;
  //! Validate a sublist param is indeed a sublist
  void validateEntryIsList(const std::string &name, const ParameterEntry &entry) const;
  //! Throw a sublist does not exist exception
  void validateMissingSublistMustExist(const std::string &baselist_name,
    const std::string &sublist_name, const bool mustAlreadyExist) const;
  //! Update sublist names recursively
  void updateSubListNames(int depth = 0);
  
private: // Data members

  //! Name of the (sub)list
  std::string name_;

  //! Parameter list
//use pragmas to disable some false-positive warnings for windows sharedlibs export
//#ifdef _MSC_VER
//#pragma warning(push)
//#pragma warning(disable:4251)
//#endif
  params_t params_;
//#ifdef _MSC_VER
//#pragma warning(pop)
//#endif

  //! Validate into list or not
  bool disableRecursiveValidation_;

};


/** \brief Nonmember constructor.
 *
 * \relates ParameterList
 */
inline
RCP<ParameterList> parameterList()
{
  return rcp(new ParameterList);
}


/** \brief Nonmember constructor.
 *
 * \relates ParameterList
 */
inline
RCP<ParameterList> parameterList(const std::string &name)
{
  return rcp(new ParameterList(name));
}
  

/** \brief Nonmember constructor.
 *
 * \relates ParameterList
 */
inline
RCP<ParameterList> parameterList(const ParameterList& source)
{
  return rcp(new ParameterList(source));
}


/** \brief Nonmember constructor.
 *
 * \relates ParameterList
 */
inline
RCP<ParameterList> createParameterList()
{
  return rcp(new ParameterList);
}


/** \brief Nonmember constructor.
 *
 * \relates ParameterList
 */
inline
RCP<ParameterList> createParameterList(const std::string &name)
{
  return rcp(new ParameterList(name));
}


/** \brief Traits specialization.
 *
 * \relates ParameterList
 */
template<>
class TEUCHOSPARAMETERLIST_LIB_DLL_EXPORT TypeNameTraits<ParameterList> {
public:
  static std::string name() { return "ParameterList"; }
  static std::string concreteName( const ParameterList& /*t2*/ )
    { return name(); }
};


/** \brief Returns true if two parameter lists are the same.
 *
 * \relates ParameterList
 */
TEUCHOSPARAMETERLIST_LIB_DLL_EXPORT bool operator==( const ParameterList& list1, const ParameterList& list2 );


/** \brief Returns true if two parameter lists are <b>not</tt> the same.
 *
 * \relates ParameterList
 */
inline
bool operator!=( const ParameterList& list1, const ParameterList& list2 )
{
  return !( list1 == list2 );
}


/** \brief Returns true if two parameter lists have the same values.
 *
 * Two parameter lists may have the same values but may not be identical.  For
 * example, two parameters can have the same values but not have the same
 * documentation strings or the same validators.
 *
 * \note This test respects ordering of the ParameterList entries; the same values in a different 
 *       order will result in \false.
 *
 * \relates ParameterList
 */
TEUCHOSPARAMETERLIST_LIB_DLL_EXPORT bool haveSameValues( const ParameterList& list1, const ParameterList& list2 );


// /////////////////////////////////////////////////////
// Inline and Template Function Definitions


inline
ParameterList& ParameterList::setName( const std::string &name_in )
{
  name_ = name_in;
  return *this;
}


// Set functions


template<typename T>
inline
ParameterList& ParameterList::set(
  std::string const& name_in, T const& value_in, std::string const& docString_in,
  RCP<const ParameterEntryValidator> const& validator_in
  )
{
  typedef StringIndexedOrderedValueObjectContainerBase SIOVOCB;
  const Ordinal param_idx = params_.getObjOrdinalIndex(name_in);
  if (param_idx != SIOVOCB::getInvalidOrdinal()) {
    Ptr<ParameterEntry> param = params_.getNonconstObjPtr(param_idx);
    const std::string docString =
      (docString_in.length() ? docString_in : param->docString());
    const RCP<const ParameterEntryValidator> validator =
      (nonnull(validator_in) ? validator_in : param->validator());
     // Create temp param to validate before setting
    ParameterEntry param_new(value_in, false, false, docString, validator );
    if (nonnull(validator)) {
      validator->validate(param_new, name_in, this->name());
    }
    // Strong guarantee: (if exception is thrown, the value is not changed)
    *param = param_new;
  }
  else {
    ParameterEntry param_new(value_in, false, false, docString_in, validator_in);
    if (nonnull(param_new.validator())) {
      param_new.validator()->validate(param_new, name_in, this->name());
    }
    params_.setObj(name_in, param_new);
  }
  return *this;
}


inline
ParameterList& ParameterList::set(
  std::string const& name_in, char value[], std::string const& docString
  ,RCP<const ParameterEntryValidator> const& validator
  ) 
{ return set(name_in, std::string(value), docString, validator); }


inline
ParameterList& ParameterList::set(
  const std::string& name_in, const char value[], const std::string &docString
  ,RCP<const ParameterEntryValidator> const& validator
  ) 
{ return set( name_in, std::string(value), docString, validator ); }


inline
ParameterList& ParameterList::set(
  std::string const& name_in, ParameterList const& value, std::string const& /*docString*/
  )
{
  sublist(name_in) = value;
  return *this;
}


inline
ParameterList& ParameterList::setEntry(std::string const& name_in, ParameterEntry const& entry_in)
{
  params_.setObj(name_in, entry_in);
  return *this;
}


// Get functions


template<typename T>
T& ParameterList::get(const std::string& name_in, T def_value)
{
  typedef StringIndexedOrderedValueObjectContainerBase SIOVOCB;
  Ordinal param_idx = params_.getObjOrdinalIndex(name_in); 
  if (param_idx == SIOVOCB::getInvalidOrdinal()) {
    // Param does not exist
    param_idx = params_.setObj(name_in, ParameterEntry(def_value, true));
  }
  Ptr<ParameterEntry> param = params_.getNonconstObjPtr(param_idx);
  this->template validateEntryType<T>("get", name_in, *param);
  return getValue<T>(*param);
}


inline
std::string& ParameterList::get(const std::string& name_in, char def_value[])
{ return get(name_in, std::string(def_value)); }


inline
std::string& ParameterList::get(const std::string& name_in, const char def_value[])
{ return get(name_in, std::string(def_value)); }


template<typename T>
T& ParameterList::get(const std::string& name_in) 
{
  ParameterEntry *foundEntry = this->getEntryPtr(name_in);
  validateEntryExists("get",name_in,foundEntry);
  this->template validateEntryType<T>("get",name_in,*foundEntry);
  return getValue<T>(*foundEntry);
}

  
template<typename T>
const T& ParameterList::get(const std::string& name_in) const
{
  const ParameterEntry *foundEntry = this->getEntryPtr(name_in);
  validateEntryExists("get",name_in,foundEntry);
  this->template validateEntryType<T>("get",name_in,*foundEntry);
  return getValue<T>(*foundEntry);
}


template<typename T>
inline
T* ParameterList::getPtr(const std::string& name_in) 
{
  typedef StringIndexedOrderedValueObjectContainerBase SIOVOCB;
  const Ordinal param_idx = params_.getObjOrdinalIndex(name_in);
  if (param_idx != SIOVOCB::getInvalidOrdinal()) {
    const Ptr<ParameterEntry> param_ptr = params_.getNonconstObjPtr(param_idx);
    if (param_ptr->isType<T>()) {
      return &param_ptr->getValue<T>(0);
    }
    // Note: The above is inefficinet.  You have to do the dynamic_cast twice
    // (once to see if it is the type and once to do the cast).  This could be
    // made more efficinet by upgrading Teuchos::any to add a any_cast_ptr()
    // function but I don't think anyone actually uses this function.
    return 0;
  }
  return 0;
}

  
template<typename T>
inline
const T* ParameterList::getPtr(const std::string& name_in) const
{
  typedef StringIndexedOrderedValueObjectContainerBase SIOVOCB;
  const Ordinal param_idx = params_.getObjOrdinalIndex(name_in);
  if (param_idx != SIOVOCB::getInvalidOrdinal()) {
    const Ptr<const ParameterEntry> param_ptr = params_.getObjPtr(param_idx);
    if (param_ptr->isType<T>()) {
      return &param_ptr->getValue<T>(0);
    }
    // Note: The above is inefficinet, see above non-const getPtr() function.
    return 0;
  }
  return 0;
}


inline
ParameterEntry& ParameterList::getEntry(const std::string& name_in)
{
  ParameterEntry *foundEntry = this->getEntryPtr(name_in);
  validateEntryExists("get", name_in, foundEntry);
  return *foundEntry;
}

  
inline
const ParameterEntry& ParameterList::getEntry(const std::string& name_in) const
{
  const ParameterEntry *foundEntry = this->getEntryPtr(name_in);
  validateEntryExists("get", name_in, foundEntry);
  return *foundEntry;
}


inline
ParameterEntry*
ParameterList::getEntryPtr(const std::string& name_in)
{
  typedef StringIndexedOrderedValueObjectContainerBase SIOVOCB;
  const Ordinal param_idx = params_.getObjOrdinalIndex(name_in);
  if (param_idx != SIOVOCB::getInvalidOrdinal()) {
    return &*params_.getNonconstObjPtr(param_idx);
  }
  return 0;
}


inline
const ParameterEntry*
ParameterList::getEntryPtr(const std::string& name_in) const
{
  typedef StringIndexedOrderedValueObjectContainerBase SIOVOCB;
  const Ordinal param_idx = params_.getObjOrdinalIndex(name_in);
  if (param_idx != SIOVOCB::getInvalidOrdinal()) {
    return &*params_.getObjPtr(param_idx);
  }
  return 0;
}


inline RCP<ParameterEntry>
ParameterList::getEntryRCP(const std::string& name_in)
{
  typedef StringIndexedOrderedValueObjectContainerBase SIOVOCB;
  const Ordinal param_idx = params_.getObjOrdinalIndex(name_in);
  if (param_idx != SIOVOCB::getInvalidOrdinal()) {
    return rcpFromPtr(params_.getNonconstObjPtr(param_idx));
  }
  return null;
}


inline RCP<const ParameterEntry>
ParameterList::getEntryRCP(const std::string& name_in) const
{
  typedef StringIndexedOrderedValueObjectContainerBase SIOVOCB;
  const Ordinal param_idx = params_.getObjOrdinalIndex(name_in);
  if (param_idx != SIOVOCB::getInvalidOrdinal()) {
    return rcpFromPtr(params_.getObjPtr(param_idx));
  }
  return null;
}


// Attribute Functions


inline
const std::string& ParameterList::name() const
{
  return name_;
}

  
#ifndef DOXYGEN_SHOULD_SKIP_THIS
template<typename T>
bool ParameterList::isType(const std::string& name_in, T* /*ptr*/) const
{
  typedef StringIndexedOrderedValueObjectContainerBase SIOVOCB;
  const Ordinal param_idx = params_.getObjOrdinalIndex(name_in);
  if (param_idx != SIOVOCB::getInvalidOrdinal()) {
    return params_.getObjPtr(param_idx)->isType<T>();
  }
  return false;
}
#endif

  
template<typename T>
bool ParameterList::isType(const std::string& name_in) const
{
  return this->isType(name_in, static_cast<T*>(0));
}


// Read-only access to the iterator


inline ParameterList::ConstIterator ParameterList::begin() const
{
  return params_.begin();
}


inline ParameterList::ConstIterator ParameterList::end() const
{
  return params_.end();
}


inline const std::string& ParameterList::name(ConstIterator i) const
{
  return (i->first);
}


inline const ParameterEntry& ParameterList::entry(ConstIterator i) const
{
  return (i->second);
}


// private


inline ParameterList::Iterator ParameterList::nonconstBegin()
{
  return params_.nonconstBegin();
}


inline ParameterList::Iterator ParameterList::nonconstEnd()
{
  return params_.nonconstEnd();
}


inline ParameterEntry& ParameterList::nonconstEntry(Iterator i)
{
  return (i->second);
}


template<typename T>
inline
void ParameterList::validateEntryType(
  const std::string &/*funcName*/, const std::string &name_in,
  const ParameterEntry &entry_in
  ) const
{
  TEUCHOS_TEST_FOR_EXCEPTION_PURE_MSG(
    entry_in.getAny().type() != typeid(T), Exceptions::InvalidParameterType
    ,"Error!  An attempt was made to access parameter \""<<name_in<<"\""
    " of type \""<<entry_in.getAny().typeName()<<"\""
    "\nin the parameter (sub)list \""<<this->name()<<"\""
    "\nusing the incorrect type \""<<TypeNameTraits<T>::name()<<"\"!"
    );
}


// //////////////////////////////////////
// Helper functions

  
/*! \relates ParameterList
  \brief A templated helper function for getting a parameter from a non-const list.
  This helper function prevents the need for giving a nominal value of the specific template type.
    
  \note The syntax for calling this function is:  <tt> getParameter<int>( list, "Iters" ) </tt>
*/
template<typename T>
T& getParameter( ParameterList& l, const std::string& name )
{
  return l.template get<T>(name);
}

  
/*! \relates ParameterList
  \brief A shorter name for <tt>getParameter()</tt>.
    
  \note The syntax for calling this function is:  <tt> get<int>( list, "Iters" ) </tt>
*/
template<typename T>
inline
T& get( ParameterList& l, const std::string& name )
{
  return getParameter<T>(l,name);
}

  
/*! \relates ParameterList
  \brief A templated helper function for getting a parameter from a const list.
  This helper function prevents the need for giving a nominal value of the specific template type.
    
  \note The syntax for calling this function is:  <tt> getParameter<int>( list, "Iters" ) </tt>    
*/
template<typename T>
const T& getParameter( const ParameterList& l, const std::string& name )
{
  return l.template get<T>(name);
}

  
/*! \relates ParameterList
  \brief A templated helper function for getting a pointer to a parameter from
  a non-const list, if it exists.  This helper function prevents the need for
  giving a nominal value of the specific template type.
  \note The syntax for calling this function is:
  <tt>getParameterPtr<int>(list,"Iters")</tt>
*/
template<typename T>
inline
T* getParameterPtr( ParameterList& l, const std::string& name )
{
  return l.template getPtr<T>(name);
}

  
/*! \relates ParameterList
  \brief A templated helper function for getting a pointer to a parameter from
  a non-const list, if it exists.  This helper function prevents the need for
  giving a nominal value of the specific template type.
  \note The syntax for calling this function is:
  <tt>getParameterPtr<int>(list,"Iters")</tt>
*/
template<typename T>
inline
const T* getParameterPtr( const ParameterList& l, const std::string& name )
{
  return l.template getPtr<T>(name);
}

  
/*! \relates ParameterList
  \brief A templated helper function for determining the type of a parameter entry for a non-const list.  
  This helper function avoids the need for giving a nominal value of the specific template type.
    
  \note The syntax for calling this function is:  <tt> isParameterType<int>( list, "Iters" ) </tt>
*/
template<typename T>
inline
bool isParameterType( ParameterList& l, const std::string& name )
{
  return l.isType( name, (T*)NULL );
}

  
/*! \relates ParameterList
  \brief A templated helper function for determining the type of a parameter entry for a const list.  
  This helper function avoids the need for giving a nominal value of the specific template type.
    
  \note The syntax for calling this function is:  <tt> isParameterType<int>( list, "Iters" ) </tt>
*/
template<typename T>
inline
bool isParameterType( const ParameterList& l, const std::string& name )
{
  return l.isType( name, (T*)NULL );
}

  
/** \brief Set a std::string parameter representation of an array.
 *
 * \param paramName [in] The name of the parameter containing the std::string
 * representation of the array.
 *
 * \param array [in] The array that will be set as a std::string parameter.
 *
 * \param paramList [in/out] The parameter list that the array will be set on.
 *
 * \relates ParameterList
 */
template<typename T>
void setStringParameterFromArray(
  const std::string          &paramName
  ,const Array<T>       &array
  ,ParameterList        *paramList
  )
{
  TEUCHOS_TEST_FOR_EXCEPT(!paramList);
  paramList->set(paramName,toString(array));
}

  
/** \brief Get an Array object (with entries of type <tt>T</tt>) from a
 * parameter holding a std::string representation of the array.
 *
 * \param paramList [in] The parameter list to extract the parameter array
 * from.
 *
 * \param paramName [in] The name of the parameter containing the std::string
 * representation of the array.
 *
 * \param arrayDim [in] If <tt>arrayDim >= 0</tt>, then the read in array must
 * be equal to this dimension, or an std::exception will be thrown.  If
 * <tt>arrayDim < 0</tt>, then an array of any dimension will be returned.
 * The default is <tt>-1</tt> and therefore no array length validation will be
 * performed.
 *
 * \param mustExist [in] If <tt>mustExist==true</tt>, then the parameter
 * <tt>paramName</tt> must exist and must contain a valid array, or an
 * std::exception is thrown.  If <tt>mustExist==false</tt>, and if the
 * parameter <tt>paramName</tt> does not exist or contains an empty array
 * std::string value, then an empty array object will be returned.
 *
 * \returns an array object if an std::exception is not thrown.  If
 * <tt>mustExist==false</tt> and the parameter does not exist, then an empty
 * array object will be returned.  If <tt>mustExist==true</tt> and
 * <tt>arrayDim < 0</tt>, then if the parameter <tt>paramName</tt> exists and
 * its array value is valid, then the converted array, of any size, will be
 * returned.  If <tt>mustExist==true</tt> and <tt>arrayDim >= 0</tt> then an
 * array of dimension <tt>arrayDim</tt> will be returned if an std::exception is
 * not thrown.
 *
 * <b>Exceptions:</b>
 *
 * <ul>
 *
 * <li><tt>Exceptions::InvalidParameterName</tt> will be thrown if
 * <tt>mustExist==true</tt> and the parameter <tt>paramName</tt> does not
 * exist in <tt>paramList</tt>
 *
 * <li><tt>Exceptions::InvalidParameterType</tt> will be thrown if the
 * parameter exists but does not have a value type of <tt>std::string</tt>.
 *
 * <li><tt>Exceptions::InvalidParameterValue</tt> will be thrown in the following cases:
 *
 *   <ul>
 *   <li>If the parameter <tt>paramName</tt> exists but the array in std::string form
 *       is not formated correctly.
 *   <li>If <tt>arrayDim >= 0</tt> and the read in array dimension dies not equal
 *       <tt>arrayDim</tt>
 *   </ul>
 *
 * </ul>
 *
 * <b>Detailed Description:</b>
 *
 * This function allows <tt>Array<T></tt> objects to be read in from a
 * parameter with a std::string representation of the array.  The templated function
 * <tt>Teuchos::fromStringToArray()</tt> (see documentation for
 * <tt>Teuchos::Array</tt>) is used to parse the std::string representation and
 * return the array object (see this function's documentation for details on
 * what the formatting of the array std::string must be and what can be handled and
 * what can not be handled.
 *
 * \relates ParameterList
 */
template<typename T>
Array<T> getArrayFromStringParameter(
  const ParameterList   &paramList
  ,const std::string         &paramName
  ,const int            arrayDim        = -1
  ,const bool           mustExist       = true
  )
{
  std::string arrayStr;
  if(mustExist) {
    arrayStr = getParameter<std::string>(paramList,paramName);
  }
  else {
    const std::string
      *arrayStrPtr = getParameterPtr<std::string>(paramList,paramName);
    if(arrayStrPtr) {
      arrayStr = *arrayStrPtr;
    }
    else {
      return Array<T>(); // Return an empty array
    }
  }
  Array<T> a;
  try {
    a = fromStringToArray<T>(arrayStr);
  }
  catch( const InvalidArrayStringRepresentation&) {
    TEUCHOS_TEST_FOR_EXCEPTION_PURE_MSG(
      true, Exceptions::InvalidParameterValue
      ,"Error!  The parameter \""<<paramName<<"\"\n"
      "in the sublist \""<<paramList.name()<<"\"\n"
      "exists, but the std::string value:\n"
      "----------\n"
      <<arrayStr<<
      "\n----------\n"
      "is not a valid array represntation!"
      );
  }
  TEUCHOS_TEST_FOR_EXCEPTION_PURE_MSG(
    ( ( a.size()>0 && arrayDim>=0 ) && static_cast<int>(a.size())!=arrayDim )
    ,Exceptions::InvalidParameterValue
    ,"Error!  The parameter \""<<paramName<<"\"\n"
    "in the sublist \""<<paramList.name()<<"\"\n"
    "exists and is a valid array, but the dimension of\n"
    "the read in array a.size() = " << a.size() << "\n"
    "was not equal to the expected size arrayDim = " << arrayDim << "!"
    );
  return a;
}


/*! \relates ParameterList
  \brief Return a RCP to a sublist in another RCP-ed parameter list.
*/
inline
RCP<ParameterList> sublist(
  const RCP<ParameterList> &paramList, const std::string& name,
  bool mustAlreadyExist = false, const std::string& docString = ""
  )
{
  return rcpWithEmbeddedObjPostDestroy(
    &paramList->sublist(name, mustAlreadyExist, docString), paramList, false );
}


/*! \relates ParameterList
  \brief Return a RCP to a sublist in another RCP-ed parameter list.
*/
inline
RCP<const ParameterList> sublist(
  const RCP<const ParameterList> &paramList, const std::string& name
  )
{
  return rcpWithEmbeddedObjPostDestroy(
    &paramList->sublist(name), paramList, false );
}

  
/*! \relates ParameterList
  \brief Output stream operator for handling the printing of the parameter list.
*/
inline std::ostream& operator<<(std::ostream& os, const ParameterList& l)
{
  return l.print(os);
}

  
} // end of Teuchos namespace


#endif
