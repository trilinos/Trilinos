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
#include "Teuchos_Array.hpp"
#include "Teuchos_StrUtils.hpp"

namespace Teuchos {

/** \brief Standard implementation of a ParameterEntryValidator that maps from
 * a list of strings to some integral type.
 */
template<class IntegralType>
class StringToIntegralParameterEntryValidator : public ParameterEntryValidator {
public:

  //@}

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
    Array<std::string>   const& strings
    ,std::string         const& defaultParameterName
    );

  /** \brief Construct with a mapping from strings to ordinals <tt>0</tt> to
   * </tt>n-1</tt>.
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
    Array<std::string>    const& strings
    ,Array<IntegralType>  const& integralValues 
    ,std::string          const& defaultParameterName
    );

  //@}

  /** \name Local non-virtual lookup functions */
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
    const std::string &str, const std::string &paramName = ""
    ,const std::string &sublistName = ""
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
    const ParameterEntry &entry, const std::string &paramName = ""
    ,const std::string &sublistName = "", const bool activeQuery = true
    ) const;

  /** \brief Lookup a parameter from a parameter list, perform a mapping from
   * a string value embedded in the <tt>ParameterEntry</tt> object and return
   * its associated integral value.
   */
  IntegralType getIntegralValue(
    ParameterList &paramList, const std::string &paramName = ""
    ,const std::string &defaultValue = ""
    ) const;

  /** \brief Lookup a parameter from a parameter list, validate the string
   * value, and return the string value.
   */
  std::string getStringValue(
    ParameterList &paramList, const std::string &paramName = ""
    ,const std::string &defaultValue = ""
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
    const std::string &str, const std::string &paramName = ""
    ,const std::string &sublistName = ""
    ) const;

  //@}

  /** \name Overridden from ParameterEntryValidator */
  //@{

  /** \brief . */
  void printDoc(
    std::string         const& docString
    ,std::ostream            & out
    ) const;

  /** \brief . */
  Teuchos::RefCountPtr<const Array<std::string> >
  validStringValues() const;

  /** \brief . */
  void validate(
    ParameterEntry  const& entry
    ,std::string    const& paramName
    ,std::string    const& sublistName
    ) const;

  //@}

private:

  typedef std::map<std::string,IntegralType> map_t;
  std::string                             defaultParameterName_;
  std::string                             validValues_;
  RefCountPtr<const Array<std::string> >  validStringValues_;
  map_t                                   map_;

  void setValidValues( Array<std::string> const& strings );

  // Not defined and not to be called.
  StringToIntegralParameterEntryValidator();

};

// ///////////////////////////
// Implementations

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
      ,"Teuchos::StringToIntegralParameterEntryValidator::StringToIntegralParameterEntryValidator(...):"
      << "\n\nError, the string \"" << strings[i] << "\" is a duplicate for parameter \""
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
      ,"Teuchos::StringToIntegralParameterEntryValidator::StringToIntegralParameterEntryValidator(...):"
      << "\n\nError, the string \"" << strings[i] << "\" is a duplicate for parameter \""
      << defaultParameterName_ << "\""
      );
  }
  setValidValues(strings);
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
  TEST_FOR_EXCEPTION(
    itr == map_.end(), Exceptions::InvalidParameterValue
    ,"Teuchos::StringToIntegralParameterEntryValidator::getIntegralValue(\""<<str<<"\",...):"
    << "\n\nError, the value \"" << str << "\" is not recognized for the parameter \""
    << ( paramName.length() ? paramName : defaultParameterName_ )
    << "\" in the sublist \"" << sublistName << "\"."
    << "\n\nValid selections include: " << validValues_  << "."
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
  TEST_FOR_EXCEPTION(
    !validType, Exceptions::InvalidParameterType
    ,"Teuchos::StringToIntegralParameterEntryValidator::getIntegralValue(...):"
    << "\n\nError, the parameter {paramName=\""<<(paramName.length()?paramName:defaultParameterName_)
    << "\",type=\""<<entry.getAny(activeQuery).typeName()<<"\"}"
    << "\nin the sublist \"" << sublistName << "\""
    << "\nhas the wrong type.  The correct type is \"string\"!"
    );
  const std::string
    &strValue = any_cast<std::string>(entry.getAny(activeQuery)); // This cast should not fail!
  return getIntegralValue(strValue); // This will validate the value!
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
  out << "#   Valid string values: " << validValues_ << ".\n";
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
  Array<std::string> const& strings
  )
{
  validStringValues_ = rcp(new Array<std::string>(strings));
  // Here I build the list of valid values in the same order as passed in by
  // the client!
  std::ostringstream oss;
  typename map_t::const_iterator itr = map_.begin();
  for( int i = 0; i < static_cast<int>(strings.size()); ++i ) {
    if(i > 0) oss << ", ";
    oss << "\""<<strings[i]<<"\"";
  }
  validValues_ = oss.str();
}

} // namespace Teuchos

#endif // TEUCHOS_STANDARD_PARAMETER_ENTRY_VALIDATORS_H
