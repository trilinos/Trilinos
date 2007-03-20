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

#ifndef TEUCHOS_STANDARD_PARAMETER_ENTRY_VALIDATORS_H
#define TEUCHOS_STANDARD_PARAMETER_ENTRY_VALIDATORS_H

#include "Teuchos_ParameterEntryValidator.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_ParameterListExceptions.hpp"
#include "Teuchos_VerbosityLevel.hpp"
#include "Teuchos_Array.hpp"
#include "Teuchos_StrUtils.hpp"


namespace Teuchos {


/** \brief Standard implementation of a ParameterEntryValidator that maps from
 * a list of strings to some integral type value.
 *
 * Objects of this type are meant to be used as both abstract objects passed
 * to <tt>Teuchos::ParameterList</tt> objects to be used to validate parameter
 * types and values, and to be used by the code that reads parameter values.
 * Having a single definition for the types of valids input and outputs for a
 * parameter value makes it easier to write error free validated code.
 */
template<class IntegralType>
class StringToIntegralParameterEntryValidator : public ParameterEntryValidator {
public:

  /** \name Constructors */
  //@{

  /** \brief Construct with a mapping from strings to ordinals <tt>0</tt> to
   * </tt>n-1</tt>.
   *
   * \param  strings
   *             [in] Array of unique string names.
   * \param  defaultParameterName
   *             [in] The default name of the parameter (used in error messages)
   */
  StringToIntegralParameterEntryValidator(
    Array<std::string> const& strings,
    std::string const& defaultParameterName
    );

  /** \brief Construct with a mapping from strings to aribitrary typed
   * integral values.
   *
   * \param  strings
   *             [in] Array of unique string names.
   * \param  integralValues
   *            [in] Array that gives the integral values associated with
   *            <tt>strings[]</tt>
   * \param  defaultParameterName
   *             [in] The default name of the parameter (used in error messages)
   *
   * <b>Preconditions:</b><ul>
   * <li> <tt>strings.size() == integralValues.size()</tt>
   * </ul>
   */
  StringToIntegralParameterEntryValidator(
    Array<std::string> const& strings,
    Array<IntegralType> const& integralValues, 
    std::string const& defaultParameterName
    );

  /** \brief Construct with a mapping from strings (with documentation) to
   * aribitrary typed integral values.
   *
   * \param  strings
   *             [in] Array of unique string names.
   * \param  stringsDocs
   *             [in] Array of documentation strings for each string value.
   * \param  integralValues
   *            [in] Array that gives the integral values associated with
   *            <tt>strings[]</tt>
   * \param  defaultParameterName
   *             [in] The default name of the parameter (used in error messages)
   *
   * <b>Preconditions:</b><ul>
   * <li> <tt>strings.size() == stringDocs.size()</tt>
   * <li> <tt>strings.size() == integralValues.size()</tt>
   * </ul>
   */
  StringToIntegralParameterEntryValidator(
    Array<std::string> const& strings,
    Array<std::string> const& stringsDocs,
    Array<IntegralType> const& integralValues, 
    std::string const& defaultParameterName
    );

  //@}

  /** \name Local non-virtual validated lookup functions */
  //@{

  /** \brief Perform a mapping from a string value to its integral value.
   *
   * \param  str  [in] String that is being used to lookup the corresponding
   *              integral value.
   * \param  paramName
   *              [in] Optional name that will be used to generate error messages.
   *
   * If the string name <tt>str</tt> does not exist, the an exception will be
   * thrown with a very descriptive error message.
   */
  IntegralType getIntegralValue(
    const std::string &str, const std::string &paramName = "",
    const std::string &sublistName = ""
    ) const;

  /** \brief Perform a mapping from a string value embedded in a
   * <tt>ParameterEntry</tt> object and return its associated integral value.
   *
   * \param  str  [in] String that is being used to lookup the corresponding
   *              integral value.
   * \param  paramName
   *              [in] Optional name that will be used to generate error messages.
   * \param  activeQuery
   *              [in] If true, then this lookup will be recored as an active query
   *              which will turn the <tt>isUsed</tt> bool to <tt>true</tt>.
   *
   * If the string name <tt>str</tt> does not exist, the an exception will be
   * thrown with a very descriptive error message.
   */
  IntegralType getIntegralValue(
    const ParameterEntry &entry, const std::string &paramName = "",
    const std::string &sublistName = "", const bool activeQuery = true
    ) const;

  /** \brief Lookup a parameter from a parameter list, perform a mapping from
   * a string value embedded in the <tt>ParameterEntry</tt> object and return
   * its associated integral value.
   */
  IntegralType getIntegralValue(
    ParameterList &paramList, const std::string &paramName,
    const std::string &defaultValue
    ) const;

  /** \brief Lookup a parameter from a parameter list, validate the string
   * value, and return the string value.
   */
  std::string getStringValue(
    ParameterList &paramList, const std::string &paramName,
    const std::string &defaultValue
    ) const;

  /** \brief Validate the string and pass it on..
   *
   * \param  str  [in] String that is being used to lookup the corresponding
   *              integral value.
   * \param  name [in] Optional name that will be used to generate error messages.
   *
   * If the string name <tt>str</tt> does not exist, the an exception will be
   * thrown with a very descriptive error message.
   */
  std::string validateString(
    const std::string &str, const std::string &paramName = "",
    const std::string &sublistName = ""
    ) const;

  //@}

  /** \name Overridden from ParameterEntryValidator */
  //@{

  /** \brief . */
  void printDoc(
    std::string const& docString,
    std::ostream & out
    ) const;

  /** \brief . */
  Teuchos::RefCountPtr<const Array<std::string> >
  validStringValues() const;

  /** \brief . */
  void validate(
    ParameterEntry const& entry,
    std::string const& paramName,
    std::string const& sublistName
    ) const;

  //@}

private:

  typedef std::map<std::string,IntegralType> map_t;
  std::string defaultParameterName_;
  std::string validValues_;
  RefCountPtr<const Array<std::string> > validStringValues_;
  RefCountPtr<const Array<std::string> > validStringValuesDocs_;
  map_t map_;

  void setValidValues(
    Array<std::string> const& strings,
    Array<std::string> const* stringsDocs = NULL
    );

  // Not defined and not to be called.
  StringToIntegralParameterEntryValidator();

};

/** \brief Return a validator for <tt>EVerbosityLevel</tt>. */
RefCountPtr<StringToIntegralParameterEntryValidator<EVerbosityLevel> >
verbosityLevelParameterEntryValidator(std::string const& defaultParameterName);

/** \brief Standard implementation of a ParameterEntryValidator that accepts
 * numbers from a number of different formats and converts them to numbers
 * from another format.
 *
 * Objects of this type are meant to be used as both abstract objects passed
 * to <tt>Teuchos::ParameterList</tt> objects to be used to validate parameter
 * types and values, and to be used by the code that reads parameter values.
 * Having a single definition for the types of valids input and outputs for a
 * parameter value makes it easier to write error free validated code.
 */
class AnyNumberParameterEntryValidator : public ParameterEntryValidator {
public:

  /** \name Public types */
  //@{

  /** \brief Determines the types that are accepted.
   */
  class AcceptedTypes {
  public:
    /** \brief Allow all types or not on construction. */
    AcceptedTypes( bool allowAllTypesByDefault = true )
      :allowInt_(allowAllTypesByDefault),allowDouble_(allowAllTypesByDefault),
       allowString_(allowAllTypesByDefault)
      {}
    /** \brief Set allow an <tt>int</tt> value or not */
    AcceptedTypes& allowInt( bool _allowInt )
      { allowInt_ = _allowInt; return *this; }
    /** \brief Set allow a <tt>double</tt> value or not */
    AcceptedTypes& allowDouble( bool _allowDouble )
      { allowDouble_ = _allowDouble; return *this; }
    /** \brief Set allow an <tt>string</tt> value or not */
    AcceptedTypes& allowString( bool _allowString )
      { allowString_ = _allowString; return *this; }
    /** \brief Allow an <tt>int</tt> value? */
    bool allowInt() const { return allowInt_; }
    /** \brief Allow an <tt>double</tt> value? */
    bool allowDouble() const { return allowDouble_; }
    /** \brief Allow an <tt>string</tt> value? */
    bool allowString() const { return allowString_; }
  private:
    bool  allowInt_;
    bool  allowDouble_;
    bool  allowString_;
  };

  //@}

  /** \name Constructors */
  //@{

  /** \brief Construct with allowed input and output types. */
  AnyNumberParameterEntryValidator(
    AcceptedTypes const& acceptedTypes
    );

  //@}

  /** \name Local non-virtual validated lookup functions */
  //@{

  /** \brief Get an integer value from a parameter entry. */
  int getInt(
    const ParameterEntry &entry, const std::string &paramName = "",
    const std::string &sublistName = "", const bool activeQuery = true
    ) const;

  /** \brief Get a double value from a parameter entry. */
  double getDouble(
    const ParameterEntry &entry, const std::string &paramName = "",
    const std::string &sublistName = "", const bool activeQuery = true
    ) const;

  /** \brief Get a string value from a parameter entry. */
  std::string getString(
    const ParameterEntry &entry, const std::string &paramName = "",
    const std::string &sublistName = "", const bool activeQuery = true
    ) const;

  /** \brief Lookup parameter from a parameter list and return as an int
   * value.
   */
  int getInt(
    ParameterList &paramList, const std::string &paramName,
    const int defaultValue
    ) const;

  /** \brief Lookup parameter from a parameter list and return as an double
   * value.
   */
  double getDouble(
    ParameterList &paramList, const std::string &paramName,
    const double defaultValue
    ) const;

  /** \brief Lookup parameter from a parameter list and return as an string
   * value.
   */
  std::string getString(
    ParameterList &paramList, const std::string &paramName,
    const std::string &defaultValue
    ) const;
  
  //@}

  /** \name Overridden from ParameterEntryValidator */
  //@{

  /** \brief . */
  void printDoc(
    std::string const& docString,
    std::ostream & out
    ) const;

  /** \brief . */
  Teuchos::RefCountPtr<const Array<std::string> >
  validStringValues() const;

  /** \brief . */
  void validate(
    ParameterEntry const& entry,
    std::string const& paramName,
    std::string const& sublistName
    ) const;

  //@}

private:

  const AcceptedTypes acceptedTypes_;
  std::string acceptedTypesString_;

  void throwTypeError(
    ParameterEntry const& entry,
    std::string const& paramName,
    std::string const& sublistName
    ) const;

  // Not defined and not to be called.
  AnyNumberParameterEntryValidator();

};


// ///////////////////////////
// Implementations


//
// StringToIntegralParameterEntryValidator
//


// Constructors


template<class IntegralType>
StringToIntegralParameterEntryValidator<IntegralType>::StringToIntegralParameterEntryValidator(
  Array<std::string> const& strings, std::string const& defaultParameterName
  )
  :defaultParameterName_(defaultParameterName)
{
  typedef typename map_t::value_type val_t;
  for( int i = 0; i < static_cast<int>(strings.size()); ++i ) {
    const bool unique = map_.insert( val_t( strings[i], i ) ).second;
    TEST_FOR_EXCEPTION(
      !unique, std::logic_error
      ,"Error, the string \"" << strings[i] << "\" is a duplicate for parameter \""
      << defaultParameterName_ << "\"."
      );
  }
  setValidValues(strings);
}


template<class IntegralType>
StringToIntegralParameterEntryValidator<IntegralType>::StringToIntegralParameterEntryValidator(
  Array<std::string> const& strings, Array<IntegralType> const& integralValues 
  ,std::string const& defaultParameterName
  )
  :defaultParameterName_(defaultParameterName)
{
#ifdef TEUCHOS_DEBUG
  TEST_FOR_EXCEPT( strings.size() != integralValues.size() );
#endif
  typedef typename map_t::value_type val_t;
  for( int i = 0; i < static_cast<int>(strings.size()); ++i ) {
    const bool unique = map_.insert( val_t( strings[i], integralValues[i] ) ).second;
    TEST_FOR_EXCEPTION(
      !unique, std::logic_error
      ,"Error, the string \"" << strings[i] << "\" is a duplicate for parameter \""
      << defaultParameterName_ << "\""
      );
  }
  setValidValues(strings);
}


template<class IntegralType>
StringToIntegralParameterEntryValidator<IntegralType>::StringToIntegralParameterEntryValidator(
  Array<std::string>    const& strings
  ,Array<std::string>   const& stringsDocs
  ,Array<IntegralType>  const& integralValues 
  ,std::string          const& defaultParameterName
  )
  :defaultParameterName_(defaultParameterName)
{
#ifdef TEUCHOS_DEBUG
  TEST_FOR_EXCEPT( strings.size() != stringsDocs.size() );
  TEST_FOR_EXCEPT( strings.size() != integralValues.size() );
#endif
  typedef typename map_t::value_type val_t;
  for( int i = 0; i < static_cast<int>(strings.size()); ++i ) {
    const bool unique = map_.insert( val_t( strings[i], integralValues[i] ) ).second;
    TEST_FOR_EXCEPTION(
      !unique, std::logic_error
      ,"Error, the string \"" << strings[i] << "\" is a duplicate for parameter \""
      << defaultParameterName_ << "\""
      );
  }
  setValidValues(strings,&stringsDocs);
}


// Lookup functions


template<class IntegralType>
IntegralType
StringToIntegralParameterEntryValidator<IntegralType>::getIntegralValue(
  const std::string &str, const std::string &paramName
  ,const std::string &sublistName
  ) const
{
  typename map_t::const_iterator itr = map_.find(str);
  TEST_FOR_EXCEPTION_PURE_MSG(
    itr == map_.end(), Exceptions::InvalidParameterValue
    ,"Error, the value \"" << str << "\" is not recognized for the parameter \""
    << ( paramName.length() ? paramName : defaultParameterName_ ) << "\""
    << "\nin the sublist \"" << sublistName << "\"."
    << "\n\nValid values include:"
    << "\n  {\n"
    << validValues_
    << "  }"
    );
  return (*itr).second;	
}


template<class IntegralType>
IntegralType
StringToIntegralParameterEntryValidator<IntegralType>::getIntegralValue(
  const ParameterEntry &entry, const std::string &paramName
  ,const std::string &sublistName, const bool activeQuery
  ) const
{
  const bool validType = ( entry.getAny(activeQuery).type() == typeid(std::string) );
  TEST_FOR_EXCEPTION_PURE_MSG(
    !validType, Exceptions::InvalidParameterType
    ,"Error, the parameter {paramName=\""<<(paramName.length()?paramName:defaultParameterName_)
    << "\",type=\""<<entry.getAny(activeQuery).typeName()<<"\"}"
    << "\nin the sublist \"" << sublistName << "\""
    << "\nhas the wrong type."
    << "\n\nThe correct type is \"string\"!"
    );
  const std::string
    &strValue = any_cast<std::string>(entry.getAny(activeQuery)); // This cast should not fail!
  return getIntegralValue(strValue,paramName,sublistName); // This will validate the value and throw!
}


template<class IntegralType>
IntegralType
StringToIntegralParameterEntryValidator<IntegralType>::getIntegralValue(
  ParameterList &paramList, const std::string &paramName
  ,const std::string &defaultValue
  ) const
{
  const std::string
    &strValue = paramList.get(paramName,defaultValue);
  return getIntegralValue(strValue,paramName,paramList.name());
}


template<class IntegralType>
std::string
StringToIntegralParameterEntryValidator<IntegralType>::getStringValue(
  ParameterList &paramList, const std::string &paramName
  ,const std::string &defaultValue
  ) const
{
  const std::string
    &strValue = paramList.get(paramName,defaultValue);
  getIntegralValue(strValue,paramName,paramList.name()); // Validate!
  return strValue;
}


template<class IntegralType>
std::string
StringToIntegralParameterEntryValidator<IntegralType>::validateString(
  const std::string &str, const std::string &paramName
  ,const std::string &sublistName
  ) const
{
  getIntegralValue(str,paramName,sublistName); // Validate!
  return str;
}


// Overridden from ParameterEntryValidator


template<class IntegralType>
void StringToIntegralParameterEntryValidator<IntegralType>::printDoc(
  std::string         const& docString
  ,std::ostream            & out
  ) const
{
  StrUtils::printLines(out,"# ",docString);
  out << "#   Valid string values:\n";
  out << "#     {\n";
  if(validStringValuesDocs_.get()) {
    for( int i = 0; i < static_cast<int>(validStringValues_->size()); ++i ) {
      out << "#       \"" << (*validStringValues_)[i] << "\"\n";
      StrUtils::printLines(out,"#          ",(*validStringValuesDocs_)[i] );
    }
  }
  else {
    StrUtils::printLines(out,"#   ",validValues_);
    // Note: Above validValues_ has for initial spaces already so indent should
    // be correct!
  }
  out << "#     }\n";
}


template<class IntegralType>
Teuchos::RefCountPtr<const Array<std::string> >
StringToIntegralParameterEntryValidator<IntegralType>::validStringValues() const
{
  return validStringValues_;
}


template<class IntegralType>
void StringToIntegralParameterEntryValidator<IntegralType>::validate(
  ParameterEntry  const& entry
  ,std::string    const& paramName
  ,std::string    const& sublistName
  ) const
{
  this->getIntegralValue(entry,paramName,sublistName,false);
}


// private


template<class IntegralType>
void StringToIntegralParameterEntryValidator<IntegralType>::setValidValues(
  Array<std::string>   const& strings
  ,Array<std::string>  const* stringsDocs
  )
{
  validStringValues_ = rcp(new Array<std::string>(strings));
  if(stringsDocs)
    validStringValuesDocs_ = rcp(new Array<std::string>(*stringsDocs));
  // Here I build the list of valid values in the same order as passed in by
  // the client!
  std::ostringstream oss;
  typename map_t::const_iterator itr = map_.begin();
  for( int i = 0; i < static_cast<int>(strings.size()); ++i ) {
    oss << "    \""<<strings[i]<<"\"\n";
  }
  // Note: Above four spaces is designed for the error output above.
  validValues_ = oss.str();
}


} // namespace Teuchos


#endif // TEUCHOS_STANDARD_PARAMETER_ENTRY_VALIDATORS_H
