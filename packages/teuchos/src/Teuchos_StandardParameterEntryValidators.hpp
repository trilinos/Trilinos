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

#ifndef TEUCHOS_STANDARD_PARAMETER_ENTRY_VALIDATORS_H
#define TEUCHOS_STANDARD_PARAMETER_ENTRY_VALIDATORS_H

#include "Teuchos_ParameterEntryValidator.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_ParameterListExceptions.hpp"
#include "Teuchos_VerbosityLevel.hpp"
#include "Teuchos_TwoDArray.hpp"
#include "Teuchos_Assert.hpp"
#include "Teuchos_StrUtils.hpp"
#include "Teuchos_TypeNameTraits.hpp"
#include "Teuchos_DummyObjectGetter.hpp"

namespace Teuchos {


/** \brief Standard implementation of a ParameterEntryValidator that maps from
 * a list of strings to some integral type value.
 *
 * Objects of this type are meant to be used as both abstract objects passed
 * to <tt>Teuchos::ParameterList</tt> objects to be used to validate parameter
 * types and values, and to be used by the code that reads parameter values.
 * Having a single definition for the types of valids input and outputs for a
 * parameter value makes it easier to write error free validated code.
 *
 * Please see <tt>StringToIntegralValidatorXMLConverter</tt> for documenation 
 * regarding the XML representation of this validator.
 */
template<class IntegralType>
class StringToIntegralParameterEntryValidator : public ParameterEntryValidator {
public:
  /** \name Constructors */
  //@{

  /** \brief Construct with a mapping from strings to ordinals <tt>0</tt> to
   * </tt>n-1</tt>.
   *
   * \param strings [in] Array of unique std::string names.
   *
   * \param defaultParameterName [in] The default name of the parameter (used
   * in error messages)
   */
  StringToIntegralParameterEntryValidator(
    ArrayView<const std::string> const& strings,
    std::string const& defaultParameterName
    );

  /** \brief Construct with a mapping from strings to aribitrary typed
   * integral values.
   *
   * \param strings [in] Array of unique std::string names.
   *
   * \param integralValues [in] Array that gives the integral values
   * associated with <tt>strings[]</tt>
   *
   * \param defaultParameterName [in] The default name of the parameter (used
   * in error messages)
   *
   * <b>Preconditions:</b><ul>
   * <li> <tt>strings.size() == integralValues.size()</tt>
   * </ul>
   */
  StringToIntegralParameterEntryValidator(
    ArrayView<const std::string> const& strings,
    ArrayView<const IntegralType> const& integralValues, 
    std::string const& defaultParameterName
    );

  /** \brief Construct with a mapping from strings (with documentation) to
   * aribitrary typed integral values.
   *
   * \param strings [in] Array of unique std::string names.
   *
   * \param stringsDocs [in] Array of documentation strings for each
   * std::string value.
   *
   * \param integralValues [in] Array that gives the integral values
   * associated with <tt>strings[]</tt>
   *
   * \param defaultParameterName [in] The default name of the parameter (used
   * in error messages)
   *
   * <b>Preconditions:</b><ul>
   * <li> <tt>strings.size() == stringDocs.size()</tt>
   * <li> <tt>strings.size() == integralValues.size()</tt>
   * </ul>
   */
  StringToIntegralParameterEntryValidator(
    ArrayView<const std::string> const& strings,
    ArrayView<const std::string> const& stringsDocs,
    ArrayView<const IntegralType> const& integralValues, 
    std::string const& defaultParameterName
    );

  //@}

  /** \name Local non-virtual validated lookup functions */
  //@{

  /** \brief Perform a mapping from a std::string value to its integral value.
   *
   * \param str [in] String that is being used to lookup the corresponding
   * integral value.
   *
   * \param paramName [in] Optional name that will be used to generate error
   * messages.
   *
   * If the std::string name <tt>str</tt> does not exist, the an std::exception will be
   * thrown with a very descriptive error message.
   */
  IntegralType getIntegralValue(
    const std::string &str, const std::string &paramName = "",
    const std::string &sublistName = ""
    ) const;

  /** \brief Perform a mapping from a std::string value embedded in a
   * <tt>ParameterEntry</tt> object and return its associated integral value.
   *
   * \param entry [in] The std::string entry.
   *
   * \param paramName [in] Optional name that will be used to generate error
   * messages.
   *
   * \param sublistName [in] The name of the sublist.
   *
   * \param activeQuery [in] If true, then this lookup will be recored as an
   * active query which will turn the <tt>isUsed</tt> bool to <tt>true</tt>.
   */
  IntegralType getIntegralValue(
    const ParameterEntry &entry, const std::string &paramName = "",
    const std::string &sublistName = "", const bool activeQuery = true
    ) const;

  /** \brief Get and validate a std::string value embedded in a
   * <tt>ParameterEntry</tt> object.
   *
   *
   * \param entry [in] The std::string entry.
   *
   * \param paramName [in] Optional name that will be used to generate error
   * messages.
   *
   * \param sublistName [in] The name of the sublist.
   *
   * \param activeQuery [in] If true, then this lookup will be recored as an
   * active query which will turn the <tt>isUsed</tt> bool to <tt>true</tt>.
   */
  std::string getStringValue(
    const ParameterEntry &entry, const std::string &paramName = "",
    const std::string &sublistName = "", const bool activeQuery = true
    ) const;

  /** \brief Lookup a parameter from a parameter list, perform a mapping from
   * a std::string value embedded in the <tt>ParameterEntry</tt> object and return
   * its associated integral value.
   */
  IntegralType getIntegralValue(
    ParameterList &paramList, const std::string &paramName,
    const std::string &defaultValue
    ) const;

  /** \brief Lookup a parameter from a parameter list, validate the std::string
   * value, and return the std::string value.
   */
  std::string getStringValue(
    ParameterList &paramList, const std::string &paramName,
    const std::string &defaultValue
    ) const;

  /** \brief Get a pointer to the array containing all the documentation
   * strings.
   *
   * \return A point to the array containing all the documentation strings.
   */
  ValidStringsList getStringDocs() const;

  /** \brief Get the name of the default parameter for the validator.
   *
   * \return The name of the default parameter for the validator.
   */
  const std::string& getDefaultParameterName() const;

  /** \brief Validate the std::string and pass it on.
   *
   * \param str [in] String that is being used to lookup the corresponding
   * integral value.
   *
   * \param name [in] Optional name that will be used to generate error
   * messages.
   *
   * If the std::string name <tt>str</tt> does not exist, the an
   * std::exception will be thrown with a very descriptive error message.
   */
  std::string validateString(
    const std::string &str, const std::string &paramName = "",
    const std::string &sublistName = ""
    ) const;

  //@}

  /** \name Overridden from ParameterEntryValidator */
  //@{

  /** \brief . */
  const std::string getXMLTypeName() const;

  /** \brief . */
  void printDoc(
    std::string const& docString,
    std::ostream & out
    ) const;

  /** \brief . */
  ValidStringsList
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
  ValidStringsList validStringValues_;
  ValidStringsList validStringValuesDocs_;
  map_t map_;

  void setValidValues(
    ArrayView<const std::string> const& strings,
    ArrayView<const std::string> const* stringsDocs = NULL
    );

  // Not defined and not to be called.
  StringToIntegralParameterEntryValidator();

};


/** \brief Nonmember constructor (see implementation).
 *
 * \relates StringToIntegralParameterEntryValidator
 */
template<class IntegralType>
RCP<StringToIntegralParameterEntryValidator<IntegralType> >
stringToIntegralParameterEntryValidator(
  ArrayView<const std::string> const& strings,
  std::string const& defaultParameterName
  );


/** \brief Nonmember constructor (see implementation).
 *
 * \relates StringToIntegralParameterEntryValidator
 */
template<class IntegralType>
RCP<StringToIntegralParameterEntryValidator<IntegralType> >
stringToIntegralParameterEntryValidator(
  ArrayView<const std::string> const& strings,
  ArrayView<const IntegralType> const& integralValues, 
  std::string const& defaultParameterName
  );


/** \brief Nonmember constructor (see implementation).
 *
 * \relates StringToIntegralParameterEntryValidator
 */
template<class IntegralType>
RCP<StringToIntegralParameterEntryValidator<IntegralType> >
stringToIntegralParameterEntryValidator(
  ArrayView<const std::string> const& strings,
  ArrayView<const std::string> const& stringsDocs,
  ArrayView<const IntegralType> const& integralValues, 
  std::string const& defaultParameterName
  );

/** \brief Set up a std::string parameter that will use an embedded validator
 * to allow the extraction of an integral value.
 *
 * The function <tt>getIntegralValue()</tt> can then be used to extract the
 * integral value of the std::string parameter.  In this case, the integral
 * value return will just be the zero-based index of the std::string value in
 * the list <tt>strings</tt>.
 *
 * \relates ParameterList
 */
template<class IntegralType>
void setStringToIntegralParameter(
  std::string const& paramName,
  std::string const& defaultValue,
  std::string const& docString,
  ArrayView<const std::string> const& strings,
  ParameterList * paramList
  );


/** \brief Set up a std::string parameter that will use an embedded validator
 * to allow the extraction of an integral value from a list of integral
 * values.
 *
 * The function <tt>getIntegralValue()</tt> can then be used to extract the
 * integral value of the std::string parameter.  In this case, the integral
 * value return will just be the zero-based index of the std::string value in
 * the list <tt>strings</tt>.
 *
 * \relates ParameterList
 */
template<class IntegralType>
void setStringToIntegralParameter(
  std::string const& paramName,
  std::string const& defaultValue,
  std::string const& docString,
  ArrayView<const std::string> const& strings,
  ArrayView<const IntegralType> const& integralValues, 
  ParameterList * paramList
  );


/** \brief Set up a std::string parameter with documentation strings for each
 * valid value that will use an embedded validator to allow the extraction of
 * an integral value from a list of integral values.
 *
 * The function <tt>getIntegralValue()</tt> can then be used to extract the
 * integral value of the std::string parameter.  In this case, the integral
 * value return will just be the zero-based index of the std::string value in
 * the list <tt>strings</tt>.
 *
 * \relates ParameterList
 */
template<class IntegralType>
void setStringToIntegralParameter(
  std::string const& paramName,
  std::string const& defaultValue,
  std::string const& docString,
  ArrayView<const std::string> const& strings,
  ArrayView<const std::string> const& stringsDocs,
  ArrayView<const IntegralType> const& integralValues, 
  ParameterList * paramList
  );


/** \brief Get an integral value for a parameter that is assumed to already be
 * set.
 *
 * This function does a dynamic cast to get the underlying valiator of type
 * StringToIntegralParameterEntryValidator<IntegralType>.  If this dynamic
 * cast failes then an <tt>Exceptions::InvalidParameterType</tt>
 * std::exception is thrown with an excellent error message.
 *
 * \relates ParameterList
 */
template<class IntegralType>
IntegralType getIntegralValue(
  ParameterList const& paramList, std::string const& paramName
  );


/** \brief Get a std::string value for a parameter that is assumed to already
 * be set.
 *
 * This function does a dynamic cast to get the underlying valiator of type
 * StringToIntegralParameterEntryValidator<IntegralValue>.  The default type
 * for IntegralValue is int.  If this dynamic cast failes then an
 * <tt>Exceptions::InvalidParameterType</tt> std::exception is thrown with an
 * excellent error message.
 *
 * \relates ParameterList
 */
template<class IntegralType>
std::string getStringValue(
  ParameterList const& paramList, std::string const& paramName
  );


/** \brief Get a StringToIntegralParameterEntryValidator<IntegralType> object
 * out of a ParameterEntry object.
 *
 * This function with thrown of the validator does not exist.
 */
template<class IntegralType>
RCP<const StringToIntegralParameterEntryValidator<IntegralType> >
getStringToIntegralParameterEntryValidator(
  ParameterEntry const& entry, ParameterList const& paramList,
  std::string const& paramName
  );


/** \brief Return the std::string name of the verbosity level as it is
 * accepted by the verbosity level parameter.
 *
 * \relates EVerbosityLevel
 */
std::string getVerbosityLevelParameterValueName(
  const EVerbosityLevel verbLevel
  );


/** \brief Return a validator for <tt>EVerbosityLevel</tt>.
 *
 * \relates EVerbosityLevel
 */
RCP<StringToIntegralParameterEntryValidator<EVerbosityLevel> >
verbosityLevelParameterEntryValidator(std::string const& defaultParameterName);

/** \brief Speicialized class for retrieving a dummy object of type
 * StringToIntegralParameterEntryValidator<IntegralType>.
 *
 * \relates StringToIntegralParameterEntryValidator
 */
template<class IntegralType>
class DummyObjectGetter<StringToIntegralParameterEntryValidator<IntegralType> >{

public:

  /** \name Getter Functions */
  //@{

  /** \brief Retrieves a dummy object of type
  * StringToIntegralParameterEntryValidator<IntegralType>.
  */
  static RCP<StringToIntegralParameterEntryValidator<IntegralType> >
    getDummyObject();
  
  //@}
};

template<class IntegralType>
RCP<StringToIntegralParameterEntryValidator<IntegralType> >
  DummyObjectGetter<StringToIntegralParameterEntryValidator<IntegralType> >::getDummyObject()
{
  return stringToIntegralParameterEntryValidator<IntegralType>(
    tuple<std::string>(""), tuple<std::string>(""), 
    tuple<IntegralType>((IntegralType)1), "");
}



/** \brief Standard implementation of a ParameterEntryValidator that accepts
 * numbers from a number of different formats and converts them to numbers in
 * another format.
 *
 * Objects of this type are meant to be used as both abstract objects passed
 * to <tt>Teuchos::ParameterList</tt> objects to be used to validate parameter
 * types and values, and to be used by the code that reads parameter values.
 * Having a single definition for the types of valids input and outputs for a
 * parameter value makes it easier to write error-free validated code.
 *
 * Please see <tt>AnyNumberValidatorXMLConverter</tt> for documenation 
 * regarding the XML representation of this validator.
 */
class TEUCHOS_LIB_DLL_EXPORT AnyNumberParameterEntryValidator
 : public ParameterEntryValidator
{
public:

  /** \name Public types */
  //@{

  /** \brief Determines what type is the preferred type. */
  enum EPreferredType { PREFER_INT, PREFER_DOUBLE, PREFER_STRING };


  /** \brief Determines the types that are accepted.  */
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
    /** \brief Set allow an <tt>std::string</tt> value or not */
    AcceptedTypes& allowString( bool _allowString )
      { allowString_ = _allowString; return *this; }
    /** \brief Allow an <tt>int</tt> value? */
    bool allowInt() const { return allowInt_; }
    /** \brief Allow an <tt>double</tt> value? */
    bool allowDouble() const { return allowDouble_; }
    /** \brief Allow an <tt>std::string</tt> value? */
    bool allowString() const { return allowString_; }
  private:
    bool  allowInt_;
    bool  allowDouble_;
    bool  allowString_;
  };

  //@}

  /** \name Constructors*/
  //@{

  /** \brief Construct with a preferrded type of double and accept all
   * types.
   */
  AnyNumberParameterEntryValidator();

  /** \brief Construct with allowed input and output types and the preferred
   * type.
   *
   * \param preferredType [in] Determines the preferred type.  This enum value
   * is used to set the default value in the override
   * <tt>validateAndModify()</tt>.
   *
   * \param acceptedType [in] Determines the types that are allowed in the
   * parameter list.
   */
  AnyNumberParameterEntryValidator(
    EPreferredType const preferredType,
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

  /** \brief Get a std::string value from a parameter entry. */
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

  /** \brief Lookup parameter from a parameter list and return as an std::string
   * value.
   */
  std::string getString(
    ParameterList &paramList, const std::string &paramName,
    const std::string &defaultValue
    ) const;

  /** \brief Lookup whether or not Doubles are allowed.
   */
  bool isDoubleAllowed() const;

  /** \brief Lookup whether or not ints are allowed.
   */
  bool isIntAllowed() const;
  
  /** \brief Lookup whether or not strings are allowed.
   */
  bool isStringAllowed() const;

  /** \brief Lookup the preferred type
   * */
  EPreferredType getPreferredType() const;

  /** \brief Gets the string representation of a given preferred type enum. */
  static const std::string& 
    getPrefferedTypeString(EPreferredType enumValue)
  {
    switch(enumValue){
      case PREFER_INT:
        return getIntEnumString();
        break;
      case PREFER_DOUBLE:
        return getDoubleEnumString();
        break;
      case PREFER_STRING:
        return getStringEnumString();
        break;
      default:
        static const std::string typeString(toString(enumValue));
        throw std::runtime_error("Cannot convert enumValue: " + typeString + " to a string");
        //Should never get here. This code is here so that a warning is not generated.
        return typeString;
    }
  }

  /** \brief Gets the preferred type enum associated with a give string. */
  static EPreferredType getPrefferedTypeStringEnum(const std::string& enumString)
  {
    if(enumString == getIntEnumString()){
      return PREFER_INT;
    }
    else if(enumString == getDoubleEnumString()){
      return PREFER_DOUBLE;
    }
    else if(enumString == getStringEnumString()){
      return PREFER_STRING;
    }
    else{
      throw std::runtime_error("Cannot convert enumString: " + enumString + " to an enum");
    }
    //Should never get here. This code is here so that a warning is not generated.
    return (EPreferredType)-1;  
  }
  
  //@}

  /** \name Overridden from ParameterEntryValidator */
  //@{

  /** \brief . */
  const std::string getXMLTypeName() const;

  /** \brief . */
  void printDoc(
    std::string const& docString,
    std::ostream & out
    ) const;

  /** \brief . */
  ValidStringsList 
  validStringValues() const;

  /** \brief . */
  void validate(
    ParameterEntry const& entry,
    std::string const& paramName,
    std::string const& sublistName
    ) const;

  /** \brief . */
  void validateAndModify(
    std::string const& paramName,
    std::string const& sublistName,
    ParameterEntry * entry
    ) const;


  //@}

private:

  // ////////////////////////////
  // Private data members

  EPreferredType preferredType_;
  std::string acceptedTypesString_;

//use pragmas to disable some false-positive warnings for windows sharedlibs export
#ifdef _MSC_VER
#pragma warning(push)
#pragma warning(disable:4251)
#endif
  const AcceptedTypes acceptedTypes_;
#ifdef _MSC_VER
#pragma warning(pop)
#endif

  // ////////////////////////////
  // Private member functions

  /* \brief Gets the string representing the "int" preferred type enum */
  static const std::string& getIntEnumString(){
    static const std::string intEnumString_ = TypeNameTraits<int>::name();
    return intEnumString_;
  }

  /* \brief Gets the string representing the "double" preferred type enum */
  static const std::string& getDoubleEnumString(){
    static const std::string doubleEnumString_ = TypeNameTraits<double>::name();
    return doubleEnumString_;
  }

  /* \brief Gets the string representing the "string" preferred type enum */
  static const std::string& getStringEnumString(){
    static const std::string stringEnumString_ = TypeNameTraits<std::string>::name();
    return stringEnumString_;
  }


  void finishInitialization();

  void throwTypeError(
    ParameterEntry const& entry,
    std::string const& paramName,
    std::string const& sublistName
    ) const;

};


// Nonmember helper functions


/** \brief Nonmember constructor AnyNumberParameterEntryValidator.
 *
 * \relates AnyNumberParameterEntryValidator
 */
TEUCHOS_LIB_DLL_EXPORT RCP<AnyNumberParameterEntryValidator>
anyNumberParameterEntryValidator();


/** \brief Nonmember constructor AnyNumberParameterEntryValidator.
 *
 * \relates AnyNumberParameterEntryValidator
 */
TEUCHOS_LIB_DLL_EXPORT RCP<AnyNumberParameterEntryValidator>
anyNumberParameterEntryValidator(
  AnyNumberParameterEntryValidator::EPreferredType const preferredType,
  AnyNumberParameterEntryValidator::AcceptedTypes const& acceptedTypes
  );

/** \brief Set an integer parameter that allows for (nearly) any input
 * parameter type that is convertible to an int.
 *
 * \relates ParameterList
 */
TEUCHOS_LIB_DLL_EXPORT void setIntParameter(
  std::string const& paramName,
  int const value, std::string const& docString,
  ParameterList *paramList,
  AnyNumberParameterEntryValidator::AcceptedTypes const& acceptedTypes
  = AnyNumberParameterEntryValidator::AcceptedTypes()
  );


/** \brief Set an double parameter that allows for (nearly) any input
 * parameter type that is convertible to a double.
 *
 * \relates ParameterList
 */
TEUCHOS_LIB_DLL_EXPORT void setDoubleParameter(
  std::string const& paramName,
  double const& value, std::string const& docString,
  ParameterList *paramList,
  AnyNumberParameterEntryValidator::AcceptedTypes const& acceptedTypes
  = AnyNumberParameterEntryValidator::AcceptedTypes()
  );


/** \brief Set an numeric parameter preferred as a std::string that allows for
 * (nearly) any input parameter type that is convertible to a std::string.
 *
 * \relates ParameterList
 */
TEUCHOS_LIB_DLL_EXPORT void setNumericStringParameter(
  std::string const& paramName,
  std::string const& value, std::string const& docString,
  ParameterList *paramList,
  AnyNumberParameterEntryValidator::AcceptedTypes const& acceptedTypes
  = AnyNumberParameterEntryValidator::AcceptedTypes()
  );


/** \brief Get an integer parameter.
 *
 * If the underlying parameter type is already an integer, then all is good.
 * However, if it is not, then a AnyNumberParameterEntryValidator object is
 * looked for to extract the type correctly.  If no validator is attached to
 * the entry, then a new AnyNumberParameterEntryValidator object will be
 * created that that will allow the conversion from any supported type.
 *
 * The parameter must exist or an <tt>Exceptions::InvalidParameterName</tt>
 * object will be thrown.  The parameters type must be acceptable, or an
 * <tt>Exceptions::InvalidParameterType</tt> object will be thown.
 *
 * \relates ParameterList
 */
TEUCHOS_LIB_DLL_EXPORT int getIntParameter(
  ParameterList const& paramList, std::string const& paramName
  );


/** \brief Get double integer parameter.
 *
 * If the underlying parameter type is already a double, then all is good.
 * However, if it is not, then a AnyNumberParameterEntryValidator object is
 * looked for to extract the type correctly.  If no validator is attached to
 * the entry, then a new AnyNumberParameterEntryValidator object will be
 * created that that will allow the conversion from any supported type.
 *
 * The parameter must exist or an <tt>Exceptions::InvalidParameterName</tt>
 * object will be thrown.  The parameters type must be acceptable, or an
 * <tt>Exceptions::InvalidParameterType</tt> object will be thown.
 *
 * \relates ParameterList
 */
TEUCHOS_LIB_DLL_EXPORT double getDoubleParameter(
  ParameterList const& paramList,
  std::string const& paramName
  );


/** \brief Get std::string numeric parameter.
 *
 * If the underlying parameter type is already a std::string, then all is
 * good.  However, if it is not, then a AnyNumberParameterEntryValidator
 * object is looked for to extract the type correctly.  If no validator is
 * attached to the entry, then a new AnyNumberParameterEntryValidator object
 * will be created that that will allow the conversion from any supported
 * type.
 *
 * The parameter must exist or an <tt>Exceptions::InvalidParameterName</tt>
 * object will be thrown.  The parameters type must be acceptable, or an
 * <tt>Exceptions::InvalidParameterType</tt> object will be thown.
 *
 * \relates ParameterList
 */
TEUCHOS_LIB_DLL_EXPORT std::string getNumericStringParameter(
  ParameterList const& paramList,
  std::string const& paramName
  );

/** \brief Speicialized class for retrieving a dummy object of type
 * AnyNumberParameterEntryValidator.
 *
 * \relates AnyNumberParameterEntryValidator
 */
template<>
class TEUCHOS_LIB_DLL_EXPORT DummyObjectGetter<AnyNumberParameterEntryValidator>{

public:

  /** \name Getter Functions */
  //@{

  /** \brief Retrieves a dummy object of type
  * AnyNumberParameterEntryValidator.
  */
  static RCP<AnyNumberParameterEntryValidator > getDummyObject();
  
  //@}
  
};


/** \brief Default structure used by EnhancedNumberTraits<T> to produce a
 * compile time error when the specialization does not exist for type
 * <tt>T</tt>.
 */
template <class T>
struct UndefinedEnhancedNumberTraits{
  //! This function should not compile if there is an attempt to instantiate!
  static inline T notDefined() {
    return T::this_type_is_missing_a_specialization(); 
  }
};


/** \brief Class defining the traits of the number type being used in an
 * EnhancedNumberValidator
 *
 * This class defines some of the traits of a number type being used by an
 * EnhancedNumberValidator.  The number has the following traits:
 *
 * \li \c min Defines the minimum possible value the number type can take on.
 *
 * \li \c max Defines the maximum possible value the number type can take on.
 *
 * \li \c defaultStep Defines the default amount a value of the number type
 * should be incremented by when being incremented in a UI.
 *
 * \li \c defaultPrecision Defines the default number of decimals with which
 * the number type should be displayed in a UI. This trait is useless for
 * non-floating point number types.
 *
 * Note that simply using this class will result in compile time errors. Only
 * specializations of this class will produce valid code.
 */
template <class T>
class EnhancedNumberTraits{
public:

  /** \brief Gets the minimum possible value the number type can take on. */
  static inline T min()
    { return UndefinedEnhancedNumberTraits<T>::notDefined(); }

  /** \brief Gets the maximum possible value the number type can take on. */
  static inline T max()
    { return UndefinedEnhancedNumberTraits<T>::notDefined(); }

  /** \brief gets default amount a value of the number type should be
   * incremented by when being utilizied in a UI. */
  static inline T defaultStep()
    { return UndefinedEnhancedNumberTraits<T>::notDefined(); }

  /** \brief Gets the default precision with which the number type should be
   * displayed. */
  static inline unsigned short defaultPrecision()
     { return UndefinedEnhancedNumberTraits<T>::notDefined(); }

};


template<>
class EnhancedNumberTraits<short int>{
public:
  static inline short int min() { return std::numeric_limits<short int>::min(); }
  static inline short int max() { return std::numeric_limits<short int>::max(); }
  static inline short int defaultStep() { return 1; }
  static inline unsigned short defaultPrecision() { return 0; }
};


template<>
class EnhancedNumberTraits<short unsigned int>{
public:
  static inline short unsigned int min() { return std::numeric_limits<short unsigned int>::min(); }
  static inline short unsigned int max() { return std::numeric_limits<short unsigned int>::max(); }
  static inline short unsigned int defaultStep() { return 1; }
  static inline unsigned short defaultPrecision() { return 0; }
};


template<>
class EnhancedNumberTraits<int>{
public:
  static inline int min() { return std::numeric_limits<int>::min(); }
  static inline int max() { return std::numeric_limits<int>::max(); }
  static inline int defaultStep() { return 1; }
  static inline unsigned short defaultPrecision() { return 0; }
};


template<>
class EnhancedNumberTraits<unsigned int>{
public:
  static inline unsigned int min() { return std::numeric_limits<unsigned int>::min(); }
  static inline unsigned int max() { return std::numeric_limits<unsigned int>::max(); }
  static inline unsigned int defaultStep() { return 1; }
  static inline unsigned short defaultPrecision() { return 0; }
};


template<>
class EnhancedNumberTraits<long int>{
public:
  static inline long int min() { return std::numeric_limits<long int>::min(); }
  static inline long int max() { return std::numeric_limits<long int>::max(); }
  static inline long int defaultStep() { return 1; }
  static inline unsigned short defaultPrecision() { return 0; }
};


template<>
class EnhancedNumberTraits<long unsigned int>{
public:
  static inline long unsigned int min() { return std::numeric_limits<long unsigned int>::min(); }
  static inline long unsigned int max() { return std::numeric_limits<long unsigned int>::max(); }
  static inline long unsigned int defaultStep() { return 1; }
  static inline unsigned short defaultPrecision() { return 0; }
};


#ifdef HAVE_TEUCHOS_LONG_LONG_INT


template<>
class EnhancedNumberTraits<long long int>{
public:
  static inline long long int min() { return std::numeric_limits<long long int>::min(); }
  static inline long long int max() { return std::numeric_limits<long long int>::max(); }
  static inline long long int defaultStep() { return 1; }
  static inline unsigned short defaultPrecision() { return 0; }
};


template<>
class EnhancedNumberTraits<long long unsigned int>{
public:
  static inline long long unsigned int min() { return std::numeric_limits<long long unsigned int>::min(); }
  static inline long long unsigned int max() { return std::numeric_limits<long long unsigned int>::max(); }
  static inline long long unsigned int defaultStep() { return 1; }
  static inline unsigned short defaultPrecision() { return 0; }
};


#endif // HAVE_TEUCHOS_LONG_LONG_INT

  
template<>
class EnhancedNumberTraits<double>{
public:
  static inline double min() { return -std::numeric_limits<double>::max(); }
  static inline double max() { return std::numeric_limits<double>::max(); }
  static inline double defaultStep() { return 1; }
  static inline unsigned short defaultPrecision() { return 100; }
};


template<>
class EnhancedNumberTraits<float>{
public:
  static inline float min() { return -std::numeric_limits<float>::max(); }
  static inline float max() { return std::numeric_limits<float>::max(); }
  static inline float defaultStep() { return 1; }
  static inline unsigned short defaultPrecision() { return 100; }
};

 
/** \brief Class uesd to validate a particular type of number.
 *
 * Please see <tt>EnhancedNumberValidatorXMLConverter</tt> for documenation 
 * regarding the XML representation of this validator.
 */
template <class T>
class EnhancedNumberValidator : public ParameterEntryValidator{

public:

  /** \name Constructors/Destructor */
  //@{

  /** \brief Constructs a EnhancedNumberValidator.
   *
   * @param min The minimum acceptable value for this validator.
   *
   * @param max The maximum acceptable value for this validator.
   *
   * @param step The increments at which the values being validated should be
   * changed when incremented in a UI.
   *
   * @param precision The number of decimials places to which the values
   * validated shold be compared to the min and max and the number of decimals
   * which are displayed in a UI. This parameter is pretty much meamingless
   * for non-floating point types.
   */
  EnhancedNumberValidator(
    T min, 
    T max, 
    T step=EnhancedNumberTraits<T>::defaultStep(), 
    unsigned short precision=EnhancedNumberTraits<T>::defaultPrecision()):
    ParameterEntryValidator(),
    minVal(min), maxVal(max), step_(step), precision_(precision), 
    containsMin(true), containsMax(true){}

  /** \brief Constructs a EnhancedNumberValidator without an explicit minimum
   * or maximum.
   */
  EnhancedNumberValidator():
    ParameterEntryValidator(),
    minVal(EnhancedNumberTraits<T>::min()),
    maxVal(EnhancedNumberTraits<T>::max()),
    step_(EnhancedNumberTraits<T>::defaultStep()),
    precision_(EnhancedNumberTraits<T>::defaultPrecision()),
    containsMin(false),
    containsMax(false){}
  
  //@}
    
  //! \name Setter Functions
  //@{
  
  /** \brief Sets the minimum acceptable value for the validator.
   * 
   * @param min The desired minimum acceptable value for the validator.
   */
  void setMin(T min){
    minVal = min;
    containsMin = true;
  }

  /** \brief Sets the maximum acceptable value for the validator.
   * 
   * @param min The desired maximum acceptable value for the validator.
   */
  void setMax(T max){
    maxVal = max;
    containsMax = true;
  }

  /** \brief Sets the step being used for the validator.
   *
   * @param The step to be used for the validator.
   */
  void setStep(T step){
    step_ = step;
  }

  /** \brief Sets the precision specified for the validator.
   *
   * @param The precision specific for the validator.
   */
  void setPrecision(unsigned short precision){
    precision_ = precision;
  }
  
  //@}

  /** \name Getter Functions */
  //@{

  /** \brief Gets the minimum acceptable value for the validator.
   *
   *@return The minimum acceptable value for the validator.
   */
  T getMin() const{
    return minVal;
  }

  /** \brief Gets the maximum acceptable value for the validator.
   *
   *@return The maximum acceptable value for the validator.
   */
  T getMax() const{
    return maxVal;
  }

  /** \brief Gets the step being used for the validator.
   *
   * @return The step being used for the validator.
   */
  T getStep() const{
    return step_;
  }

  /** \brief  Gets the precision specified for the validator.
   *
   * @return The precision specific for the validator.
   */
  unsigned short getPrecision() const{
    return precision_;
  }
  
  //@}

  //! \name Attribute/Query Methods 
  //@{

  /** \brief Determines whether or not the validator has a minimum value.
   *
   * @return True if the validator has a minimum value, false otherwise.
   */
  bool hasMin() const{
    return containsMin;
  }

  /** \brief  Determines whether or not the validator has a maximum value.
   *
   * @return True if the validator has a maximum value, false otherwise.
   */ 
  bool hasMax() const{
    return containsMax;
  }
  
  //@}

  /** \name Overridden from ParameterEntryValidator */
  //@{

  /** \brief . */
  ValidStringsList validStringValues() const{
    return null;
  }

  /** \brief . */
  void validate(ParameterEntry const &entry, std::string const &paramName,
    std::string const &sublistName) const;

  /** \brief . */
  const std::string getXMLTypeName() const{
    return  "EnhancedNumberValidator(" + TypeNameTraits<T>::name()+ ")"; 
  }

  /** \brief . */
  void printDoc(std::string const &docString, std::ostream &out) const{
    StrUtils::printLines(out,"# ",docString);
    out << "#\tValidator Used: " << std::endl;
    out << "#\t\tNumber Validator" << std::endl;
    out << "#\t\tType: " << Teuchos::TypeNameTraits<T>::name() << 
      std::endl;
    out << "#\t\tMin (inclusive): " << minVal << std::endl;
    out << "#\t\tMax (inclusive): " << maxVal << std::endl;
  }
  
  //@}

private:

  /** \name Private Members */
  //@{
  
  /** \brief The minimum value accepted by the validator.
   */
  T minVal;

  /** \brief The maximum value accepted by the validator.
   */
  T maxVal;

  /** \brief The increment to use when increaseing or decreaseing the value the validator is validating.
   */
  T step_;

  /** \brief The number of decimal places with which the nubmer will be displayed in a UI. This value
   * is meaningless for non-floating point number types.
   */
  unsigned short precision_;

  /** \brief Whether or not a minimum value has been specified for this validator.
   */
  bool containsMin;

  /** \brief Whetehr or not a maximum value has been specified for this validator.
   */
  bool containsMax;
  
  //@}
  
};

template<class T>
void EnhancedNumberValidator<T>::validate(ParameterEntry const &entry, std::string const &paramName,
  std::string const &sublistName) const
{
  any anyValue = entry.getAny(true);
  const std::string &entryName = entry.getAny(false).typeName();

  TEST_FOR_EXCEPTION(anyValue.type() != typeid(T),
    Exceptions::InvalidParameterType,
    "Aww shoot! Sorry bud, but it looks like the \"" << 
    paramName << "\"" <<
    " parameter in the \"" << sublistName 
    << "\" sublist didn't quite work out." << std::endl <<
    "No need to fret though. I'm sure it's just a small mistake. "
    "Maybe the information below "<<
    "can help you figure out what went wrong." << std::endl << std::endl <<
    "Error: The value that you entered was the wrong type." << std::endl <<
    "Parameter: " << paramName << std::endl <<
    "Type specified: " << entryName << std::endl <<
    "Type accepted: " << Teuchos::TypeNameTraits<T>::name() << std::endl);

  bool isValueInRange; 
  any_cast<T>(anyValue) >= minVal && any_cast<T>(anyValue) <= maxVal
    ? isValueInRange = true : isValueInRange=false;
  TEST_FOR_EXCEPTION(!(isValueInRange),
    Exceptions::InvalidParameterValue,
    "Aww shoot! Sorry bud, but it looks like the \"" << 
    paramName << "\"" <<
    " parameter in the \"" << sublistName << 
    "\" sublist didn't quite work "
    "out." << std::endl <<
    "No need to fret though. I'm sure it's just a small mistake. "
    "Maybe the information below "<<
    "can help you figure out what went wrong." << std::endl << 
    std::endl <<
    "Error: The value that was entered doesn't fall with in " <<
    "the range set by the validator" << std::endl <<
    "Parameter: " << paramName  << std::endl <<
    "Min: " << minVal << std::endl <<
    "Max: " << maxVal << std::endl <<
    "Value entered: " << 
    (any_cast<T>(anyValue)) << std::endl << std::endl);
}

/** \brief Speicialized class for retrieving a dummy object of type
 * EnhancedNumberValidator<T>.
 *
 * \relates EnhancedNumberValidator<T>
 */
template<class T>
class DummyObjectGetter<EnhancedNumberValidator<T> >{

public:

  /** \name Getter Functions */
  //@{

  /** \brief Retrieves a dummy object of type
  * EnhancedNumberValidator<T>.
  */
  static RCP<EnhancedNumberValidator<T> > getDummyObject();
  
  //@}
};

template<class T>
RCP<EnhancedNumberValidator<T> >
  DummyObjectGetter<EnhancedNumberValidator<T> >::getDummyObject()
{
  return rcp(new EnhancedNumberValidator<T>);
}

/** \brief Validate a file name entry.
 *
 * Simply indicates that the parameter entry with this validator should
 * contain a filename.
 *
 * Please see <tt>FileNameValidatorXMLConverter</tt> for documenation 
 * regarding the XML representation of this validator.
 */
class TEUCHOS_LIB_DLL_EXPORT FileNameValidator : public ParameterEntryValidator {

public:

  /** \name Public types */
  //@{

  /** \brief The default value of the mustAlreadyExist parameter in the
   * constructor. */
  static bool mustAlreadyExistDefault() { return false; }
  
  //@}

  /** \name Constructors/Destructor */
  //@{

  /** \brief Constructs a FileNameValidator.
   *
   * @param mustAlreadyExist True if the file the user specifies should
   * already exists, false otherwise.
   */
  FileNameValidator(bool mustAlreadyExist = mustAlreadyExistDefault());
  
  //@}

  //! \name Attribute/Query Functions 
  //@{
  
  /** \brief Gets the variable describing whether or not this validator wants
   * the file that is specified to already exist.
   *
   * @return Whether or not the validator requires the file to already exist
   */
  bool fileMustExist() const;
  
  //@}

  //! \name Setter Functions
  //@{
  
  /** \brief Sets whether or not the validator requires the file to already
   * exist.
   *
   * @param shouldFileExist True if the file should already exist, false
   * otherwise.
   *
   * @return The new value of the shouldFileExist variable.
   */
  bool setFileMustExist(bool shouldFileExist);
  
  //@}

  /** \name Overridden from ParameterEntryValidator */
  //@{

  /** \brief . */
  ValidStringsList validStringValues() const;

  /** \brief . */
  void validate(
    ParameterEntry const &entry,
    std::string const &paramName, 
    std::string const &sublistName) const;

  /** \brief . */
  const std::string getXMLTypeName() const;

  /** \brief . */
  void printDoc(std::string const &docString, std::ostream &out) const;
  
  //@}

private:

  /** \name Private Members */
  //@{
  
  /** \brief Whether or not the file specified in the parameter should 
   * already exist.
   */
  bool mustAlreadyExist_;
  
  //@}

};

/** \brief Speicialized class for retrieving a dummy object of type
 * FileNameValidator.
 *
 * \relates FileNameValidator
 */
template<>
class DummyObjectGetter<FileNameValidator>{

public:

  /** \name Getter Functions */
  //@{

  /** \brief Retrieves a dummy object of type
  * FileNameValidator.
  */
  static RCP<FileNameValidator> getDummyObject();
  
  //@}
  
};

/** \brief A simple validator that only allows certain string values to be
 * choosen or simply enforces that a particular parameter have a std::string
 * for a value.
 *
 * Please see <tt>StringValidatorXMLConverter</tt> for documenation 
 * regarding the XML representation of this validator.
 */
class TEUCHOS_LIB_DLL_EXPORT StringValidator : public ParameterEntryValidator {

public:

  /** \name Constructors/Destructor */
  //@{

  /** \brief Constructs a StringValidator.
   */
  StringValidator();
  
  /** \brief Constructs a StringValidator.
   *
   * @param validStrings A list of valid string values for this validator.
   */
  StringValidator(const Teuchos::Array<std::string> &validStrings);

  //@}

  //! \name Setter Functions
  //@{
 
  /** \brief Sets the Array of valid strings and returns what the current
   * array of valid string now is.
   *
   * @param validStrings What the array for the valid strings should contain.
   *
   * @return What the arry for the valid strings now conatians.
   */
  ValidStringsList setValidStrings(
    const Teuchos::Array<std::string> &validStrings);
  
  //@}

  /** \name Overridden from ParameterEntryValidator */
  //@{

  /** \brief . */
  ValidStringsList validStringValues() const;

  /** \brief . */
  void validate(ParameterEntry const &entry, std::string const &paramName,
    std::string const &sublistName) const;

  /** \brief . */
  const std::string getXMLTypeName() const;

  /** \brief . */
  void printDoc(std::string const &docString, std::ostream &out) const;
  
  //@}

private:

  /** \name Private Members */
  //@{
  
  /** \brief An array containing a list of all the valid string values.
   */
  ValidStringsList validStrings_;
  
  //@}

};

/** \brief Speicialized class for retrieving a dummy object of type
 * StringValidator.
 *
 * \relates StringValidator
 */
template<>
class DummyObjectGetter<StringValidator>{

public:

  /** \name Getter Functions */
  //@{

  /** \brief Retrieves a dummy object of type
  * StringValidator.
  */
  static RCP<StringValidator> getDummyObject();
  
  //@}
  
};


/**
 * \brief An abstract base class for all ArrayValidators.
 */
template<class ValidatorType, class EntryType>
class AbstractArrayValidator : public ParameterEntryValidator {

public:

  /** @name Constructors */
  //@{

  /**
   * \brief Constructs an AbstractArrayValidator.
   *
   * @param prototypeValidator The prototype validator to be applied
   * to each entry in the array.
   */
  AbstractArrayValidator(RCP<const ValidatorType> prototypeValidator):
    ParameterEntryValidator(),
    prototypeValidator_(prototypeValidator){}

  //@}

  /** \name Getter Functions */
  //@{

  /** \brief Returns the prototype validator for this Array Validator */
  RCP<const ValidatorType> getPrototype() const{
    return prototypeValidator_;
  }

  //@}

  /** \name Overridden from ParameterEntryValidator */
  //@{

  /** \brief . */
  ValidStringsList validStringValues() const {
    return prototypeValidator_->validStringValues();
  }

  //@}

private:

  /** \name Private Members */
  //@{

  /** \brief The prototype validator to be applied to each entry in the Array.
   */
  RCP<const ValidatorType> prototypeValidator_;

  /** \brief Hidden default constructor. */
  AbstractArrayValidator<ValidatorType, EntryType>();
  
  //@}

};

/**
 * \brief Takes a validator, wraps it, and applies it to a TwoDArray.
 *
 * This class is a wrapper, allowing you to apply a normal validator to a
 * TwoDArray of values.  It is templated on both the validator type and the type
 * of the entries contained within the array.
 *
 * Please see <tt>TwoDArrayValidatorXMLConverter</tt> for documenation 
 * regarding the XML representation of this validator.
 *
 * \relates TwoDArray
 */
template<class ValidatorType, class EntryType>
class TwoDArrayValidator : public AbstractArrayValidator<ValidatorType, EntryType>{
public: 
  /** @name Constructor */
  //@{

  /** \brief Constructs a ArrayValidator.
   *
   * @param prototypeValidator The validator to be used on each
   * entry in the array.
   */

  TwoDArrayValidator(RCP<const ValidatorType> prototypeValidator):
    AbstractArrayValidator<ValidatorType, EntryType>(prototypeValidator){}
  
  //@}

  /** \name Overridden from ParameterEntryValidator */
  //@{

  /** \brief . */
  virtual void validate(ParameterEntry const &entry, std::string const &paramName,
    std::string const &sublistName) const;

  /** \brief . */
  const std::string getXMLTypeName() const{
    return "TwoDArrayValidator(" + 
      this->getPrototype()->getXMLTypeName() + ", " +
      TypeNameTraits<EntryType>::name() + ")";
  }

  /** \brief . */
  virtual void printDoc(std::string const &docString, std::ostream &out) const
  {
    StrUtils::printLines(out,"# ",docString);
    std::string toPrint;
    toPrint += "TwoDArrayValidator:\n";
    toPrint += "Prototype Validator:\n";
    this->getPrototype()->printDoc(toPrint, out);
  }
  
  //@}

};

template<class ValidatorType, class EntryType>
void TwoDArrayValidator<ValidatorType, EntryType>::validate(ParameterEntry const &entry, std::string const &paramName,
  std::string const &sublistName) const
{
  any anyValue = entry.getAny(true);
  const std::string &entryName = entry.getAny(false).typeName();
  TEST_FOR_EXCEPTION(anyValue.type() != typeid(TwoDArray<EntryType>),
    Exceptions::InvalidParameterType,
    "Aww shoot! Sorry bud, but it looks like the \"" <<
    paramName << "\"" <<
    " parameter in the \"" << sublistName << 
    "\" sublist didn't quite work out." << std::endl <<
    "No need to fret though. I'm sure it's just a small mistake. "
    "Maybe the information below "<<
    "can help you figure out what went wrong." << std::endl << 
    std::endl <<
    "Error: The value you entered was the wrong type." << std::endl <<
    "Parameter: " << paramName << std::endl <<
    "Type specified: " << entryName << std::endl <<
    "Type accepted: " << TypeNameTraits<TwoDArray<EntryType> >::name() <<
    std::endl << std::endl);

  TwoDArray<EntryType> extracted = 
    getValue<Teuchos::TwoDArray<EntryType> >(entry);
  RCP<const ParameterEntryValidator> prototype = this->getPrototype();
  for(int i = 0; i<extracted.getNumRows(); ++i){
    for(int j = 0; j<extracted.getNumCols(); ++j){
      ParameterEntry dummyParameter;
      dummyParameter.setValue(extracted(i,j));
      try{
        prototype->validate(
          dummyParameter, paramName, sublistName);
      }
      catch(Exceptions::InvalidParameterValue& e){
        std::stringstream oss;
        oss << "TwoDArray Validator Exception:" << std::endl <<
        "Bad Index: (" << i << "," << j << ")" << std::endl << e.what();
        throw Exceptions::InvalidParameterValue(oss.str());
      }
    }
  }
}


/** \brief Speicialized class for retrieving a dummy object of type
 * TwoDArrayValidator.
 *
 * \relates TwoDArrayValidator
 */
template<class ValidatorType, class EntryType>
class DummyObjectGetter<TwoDArrayValidator<ValidatorType, EntryType> >{

public:

  /** \name Getter Functions */
  //@{

  /** \brief Retrieves a dummy object of type
  * TwoDArrayValidator<ValidatorType, EntryType>.
  */
  static RCP<TwoDArrayValidator<ValidatorType, EntryType> > getDummyObject();
  
  //@}
  
};

template<class ValidatorType, class EntryType>
RCP<TwoDArrayValidator<ValidatorType, EntryType> >
  DummyObjectGetter<TwoDArrayValidator<ValidatorType, EntryType> >::getDummyObject()
{
  return rcp(new TwoDArrayValidator<ValidatorType, EntryType>(
    DummyObjectGetter<ValidatorType>::getDummyObject()));
}

/** \brief Convience class for StringValidators that are to be applied to
 * TwoDArrays.
 */
class TEUCHOS_LIB_DLL_EXPORT TwoDArrayStringValidator : 
  public TwoDArrayValidator<StringValidator, std::string>{

public:

  /** \name Constructors/Destructor */
  //@{

  /** \brief . */
  TwoDArrayStringValidator(RCP<const StringValidator> prototypeValidator):
    TwoDArrayValidator<StringValidator, std::string>(prototypeValidator){}
  
  //@}

};


/** \brief Convience class for FileNameValidators that are to be applied to
 * TwoDArrays.
 *
 */
class TEUCHOS_LIB_DLL_EXPORT TwoDArrayFileNameValidator : 
  public TwoDArrayValidator<FileNameValidator, std::string>{

public:

  /** \name Constructors/Destructor */
  //@{

  /** \brief . */
  TwoDArrayFileNameValidator(RCP<const FileNameValidator> prototypeValidator):
    TwoDArrayValidator<FileNameValidator, std::string>(prototypeValidator){}
  
  //@}

};


/** \brief Convience class for EnhancedNumberValidators that are to be applied
 * to TwoDArray.
 */
template<class T>
class TwoDArrayNumberValidator : public TwoDArrayValidator<EnhancedNumberValidator<T>, T>{
public:
  /** \name Constructors/Destructor */
  //@{

  /** \brief . */
  TwoDArrayNumberValidator(
    RCP<const EnhancedNumberValidator<T> > prototypeValidator):
    TwoDArrayValidator<EnhancedNumberValidator<T>, T>(prototypeValidator){}
  
  //@}
  
};


/** \brief Takes a validator, wraps it, and applies it to an array.
 *
 * This class is a wrapper, allowing you to apply a normal validator to an
 * array of values.  It is templated on both the validator type and the type
 * of the entries contained within the array.
 *
 * Please see <tt>ArrayValidatorXMLConverter</tt> for documenation 
 * regarding the XML representation of this validator.
 *
 * \relates Array
 */
template<class ValidatorType, class EntryType>
class ArrayValidator : public AbstractArrayValidator<ValidatorType, EntryType>{

public:

  /** \name Constructors/Destructor */
  //@{

  /** \brief Constructs a ArrayValidator.
   *
   * @param prototypeValidator The validator to be used on each
   * entry in the array.
   */
  ArrayValidator(RCP<const ValidatorType> prototypeValidator):
    AbstractArrayValidator<ValidatorType, EntryType>(prototypeValidator){}
  
  //@}

  /** \name Overridden from ParameterEntryValidator */
  //@{

  /** \brief . */
  virtual void validate(ParameterEntry const &entry, std::string const &paramName,
    std::string const &sublistName) const;

  /** \brief . */
  const std::string getXMLTypeName() const{
    return "ArrayValidator(" + 
      this->getPrototype()->getXMLTypeName() + ", " +
      TypeNameTraits<EntryType>::name() + ")";
  }

  /** \brief . */
  virtual void printDoc(std::string const &docString, std::ostream &out) const
  {
    StrUtils::printLines(out,"# ",docString);
    std::string toPrint;
    toPrint += "ArrayValidator:\n";
    toPrint += "Prototype Validator:\n";
    this->getPrototype()->printDoc(toPrint, out);
  }
  
  //@}

};

template<class ValidatorType, class EntryType>
void ArrayValidator<ValidatorType, EntryType>::validate(ParameterEntry const &entry, std::string const &paramName,
  std::string const &sublistName) const
{
  any anyValue = entry.getAny(true);
  const std::string &entryName = entry.getAny(false).typeName();
  TEST_FOR_EXCEPTION(anyValue.type() != typeid(Array<EntryType>),
    Exceptions::InvalidParameterType,
    "Aww shoot! Sorry bud, but it looks like the \"" <<
    paramName << "\"" <<
    " parameter in the \"" << sublistName << 
    "\" sublist didn't quite work out." << std::endl <<
    "No need to fret though. I'm sure it's just a small mistake. "
    "Maybe the information below "<<
    "can help you figure out what went wrong." << std::endl << 
    std::endl <<
    "Error: The value you entered was the wrong type." << std::endl <<
    "Parameter: " << paramName << std::endl <<
    "Type specified: " << entryName << std::endl <<
    "Type accepted: " << TypeNameTraits<Array<EntryType> >::name() <<
    std::endl << std::endl);

  Array<EntryType> extracted = 
    getValue<Teuchos::Array<EntryType> >(entry);
  RCP<const ParameterEntryValidator> prototype = this->getPrototype();
  for(int i = 0; i<extracted.size(); ++i){
    ParameterEntry dummyParameter;
    dummyParameter.setValue(extracted[i]);
    try{
      prototype->validate(
        dummyParameter, paramName, sublistName);
    }
    catch(Exceptions::InvalidParameterValue& e){
      std::stringstream oss;
      oss << "Array Validator Exception:" << std::endl <<
      "Bad Index: " << i << std::endl << e.what();
      throw Exceptions::InvalidParameterValue(oss.str());
    }
  }
}

/** \brief Speicialized class for retrieving a dummy object of type
 * ArrayValidator.
 *
 * \relates ArrayValidator
 */
template<class ValidatorType, class EntryType>
class DummyObjectGetter<ArrayValidator<ValidatorType, EntryType> >{

public:

  /** \name Getter Functions */
  //@{

  /** \brief Retrieves a dummy object of type
  * ArrayValidator<ValidatorType, EntryType>.
  */
  static RCP<ArrayValidator<ValidatorType, EntryType> > getDummyObject();
  
  //@}
  
};

template<class ValidatorType, class EntryType>
RCP<ArrayValidator<ValidatorType, EntryType> >
  DummyObjectGetter<ArrayValidator<ValidatorType, EntryType> >::getDummyObject()
{
  return rcp(new ArrayValidator<ValidatorType, EntryType>(
    DummyObjectGetter<ValidatorType>::getDummyObject()));
}


/** \brief Convience class for StringValidators that are to be applied to
 * arrays.
 *
 * Also needed for maintaining backwards compatiblitiy with the earliest
 * versions of the Optika package.  This class would be a simple typedef,
 * however I wanted to maintain consistency with the ArrayNumberValidator
 * class which cannot be typedef'd.
 */
class TEUCHOS_LIB_DLL_EXPORT ArrayStringValidator : 
  public ArrayValidator<StringValidator, std::string>{

public:

  /** \name Constructors/Destructor */
  //@{

  /** \brief . */
  ArrayStringValidator(RCP<const StringValidator> prototypeValidator):
    ArrayValidator<StringValidator, std::string>(prototypeValidator){}
  
  //@}

};


/** \brief Convience class for FileNameValidators that are to be applied to
 * arrays.
 *
 * Also needed for maintaining backwards compatiblitiy with the earliest
 * versions of the Optika package.  This class would be a simple typedef,
 * however I wanted to maintain consistency with the ArrayNumberValidator
 * class which cannot be typedef'd.
 */
class TEUCHOS_LIB_DLL_EXPORT ArrayFileNameValidator : public ArrayValidator<FileNameValidator, std::string>{

public:

  /** \name Constructors/Destructor */
  //@{

  /** \brief . */
  ArrayFileNameValidator(RCP<const FileNameValidator> prototypeValidator):
    ArrayValidator<FileNameValidator, std::string>(prototypeValidator){}
  
  //@}

};


/** \brief Convience class for EnhancedNumberValidators that are to be applied
 * to arrays.
 *
 * Also needed for maintaining backwards compatiblitiy with the earliest
 * versions of the Optika package.  This class would be a simple typedef,
 * however the current c++ compilers do not support templated typedefs
 */
template<class T>
class ArrayNumberValidator : public ArrayValidator<EnhancedNumberValidator<T>, T>{
public:
  /** \name Constructors/Destructor */
  //@{

  /** \brief . */
  ArrayNumberValidator(
    RCP<const EnhancedNumberValidator<T> > prototypeValidator):
    ArrayValidator<EnhancedNumberValidator<T>, T>(prototypeValidator){}
  
  //@}
  
};



// ///////////////////////////
// Implementations


//
// StringToIntegralParameterEntryValidator
//


// Constructors


template<class IntegralType>
StringToIntegralParameterEntryValidator<IntegralType>::StringToIntegralParameterEntryValidator(
  ArrayView<const std::string> const& strings, std::string const& defaultParameterName
  ):
  ParameterEntryValidator(),
  defaultParameterName_(defaultParameterName)
{
  typedef typename map_t::value_type val_t;
  for( int i = 0; i < static_cast<int>(strings.size()); ++i ) {
    const bool unique = map_.insert( val_t( strings[i], (IntegralType)i ) ).second;
    TEST_FOR_EXCEPTION(
      !unique, std::logic_error
      ,"Error, the std::string \"" << strings[i] << "\" is a duplicate for parameter \""
      << defaultParameterName_ << "\"."
      );
  }
  setValidValues(strings);
}


template<class IntegralType>
StringToIntegralParameterEntryValidator<IntegralType>::StringToIntegralParameterEntryValidator(
  ArrayView<const std::string> const& strings, ArrayView<const IntegralType> const& integralValues 
  ,std::string const& defaultParameterName
  ):
  ParameterEntryValidator(),
  defaultParameterName_(defaultParameterName)
{
#ifdef TEUCHOS_DEBUG
  TEUCHOS_ASSERT_EQUALITY( strings.size(), integralValues.size() );
#endif
  TEST_FOR_EXCEPTION(
    strings.size() != integralValues.size(),
  std::logic_error,
  "Error, strings and integraValues must be of the same length."
  );
  typedef typename map_t::value_type val_t;
  for( int i = 0; i < static_cast<int>(strings.size()); ++i ) {
    const bool unique = map_.insert( val_t( strings[i], integralValues[i] ) ).second;
    TEST_FOR_EXCEPTION(
      !unique, std::logic_error
      ,"Error, the std::string \"" << strings[i] << "\" is a duplicate for parameter \""
      << defaultParameterName_ << "\""
      );
  }
  setValidValues(strings);
}

template<class IntegralType>
StringToIntegralParameterEntryValidator<IntegralType>::StringToIntegralParameterEntryValidator(
  ArrayView<const std::string>    const& strings
  ,ArrayView<const std::string>   const& stringsDocs
  ,ArrayView<const IntegralType>  const& integralValues 
  ,std::string          const& defaultParameterName
  ):
  ParameterEntryValidator(),
  defaultParameterName_(defaultParameterName)
{
#ifdef TEUCHOS_DEBUG
  TEUCHOS_ASSERT_EQUALITY( strings.size(), stringsDocs.size() );
  TEUCHOS_ASSERT_EQUALITY( strings.size(), integralValues.size() );
#endif
  TEST_FOR_EXCEPTION(
    strings.size() != integralValues.size(),
  std::logic_error,
  "Error, strings and integraValues must be of the same length."
  );
  TEST_FOR_EXCEPTION(
    strings.size() != stringsDocs.size(),
  std::logic_error,
  "Error, strings and stringsDocs must be of the same length."
  );
  typedef typename map_t::value_type val_t;
  for( int i = 0; i < static_cast<int>(strings.size()); ++i ) {
    const bool unique = map_.insert( val_t( strings[i], integralValues[i] ) ).second;
    TEST_FOR_EXCEPTION(
      !unique, std::logic_error
      ,"Error, the std::string \"" << strings[i] << "\" is a duplicate for parameter \""
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
std::string
StringToIntegralParameterEntryValidator<IntegralType>::getStringValue(
  const ParameterEntry &entry, const std::string &paramName
  ,const std::string &sublistName, const bool activeQuery
  ) const
{
  // Validate the parameter's type and value
  this->getIntegralValue(entry,paramName,sublistName,activeQuery);
  // Return the std::string value which is now validated!
  return any_cast<std::string>(entry.getAny(activeQuery)); // This cast should not fail!
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
ParameterEntryValidator::ValidStringsList
StringToIntegralParameterEntryValidator<IntegralType>::getStringDocs() const
{
  return validStringValuesDocs_; 
}

template<class IntegralType>
const std::string&
StringToIntegralParameterEntryValidator<IntegralType>::getDefaultParameterName() const
{
  return defaultParameterName_; 
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
const std::string 
StringToIntegralParameterEntryValidator<IntegralType>::getXMLTypeName() const{
  return "StringIntegralValidator(" +
    TypeNameTraits<IntegralType>::name() +
    ")";
}

template<class IntegralType>
void StringToIntegralParameterEntryValidator<IntegralType>::printDoc(
  std::string         const& docString
  ,std::ostream            & out
  ) const
{
  StrUtils::printLines(out,"# ",docString);
  out << "#   Valid std::string values:\n";
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
ParameterEntryValidator::ValidStringsList
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
  ArrayView<const std::string>   const& strings
  ,ArrayView<const std::string>  const* stringsDocs
  )
{
  validStringValues_ = rcp(new Array<std::string>(strings));
  if(stringsDocs)
    validStringValuesDocs_ = rcp(new Array<std::string>(*stringsDocs));
  // Here I build the list of valid values in the same order as passed in by
  // the client!
  std::ostringstream oss;
  for( int i = 0; i < static_cast<int>(strings.size()); ++i ) {
    oss << "    \""<<strings[i]<<"\"\n";
  }
  // Note: Above four spaces is designed for the error output above.
  validValues_ = oss.str();
}


} // namespace Teuchos


//
// Nonmember function implementations for StringToIntegralParameterEntryValidator
//


template<class IntegralType>
inline
Teuchos::RCP<Teuchos::StringToIntegralParameterEntryValidator<IntegralType> >
Teuchos::stringToIntegralParameterEntryValidator(
  ArrayView<const std::string> const& strings,
  std::string const& defaultParameterName
  )
{
  return rcp(
    new StringToIntegralParameterEntryValidator<IntegralType>(
      strings, defaultParameterName
      )
    );
}


template<class IntegralType>
inline
Teuchos::RCP<Teuchos::StringToIntegralParameterEntryValidator<IntegralType> >
Teuchos::stringToIntegralParameterEntryValidator(
  ArrayView<const std::string> const& strings,
  ArrayView<const IntegralType> const& integralValues, 
  std::string const& defaultParameterName
  )
{
  return rcp(
    new StringToIntegralParameterEntryValidator<IntegralType>(
      strings, integralValues, defaultParameterName
      )
    );
}


template<class IntegralType>
inline
Teuchos::RCP< Teuchos::StringToIntegralParameterEntryValidator<IntegralType> >
Teuchos::stringToIntegralParameterEntryValidator(
  ArrayView<const std::string> const& strings,
  ArrayView<const std::string> const& stringsDocs,
  ArrayView<const IntegralType> const& integralValues, 
  std::string const& defaultParameterName
  )
{
  return rcp(
    new StringToIntegralParameterEntryValidator<IntegralType>(
      strings, stringsDocs, integralValues, defaultParameterName
      )
    );
}

template<class IntegralType>
void Teuchos::setStringToIntegralParameter(
  std::string const& paramName,
  std::string const& defaultValue,
  std::string const& docString,
  ArrayView<const std::string> const& strings,
  ParameterList * paramList
  )
{
  typedef ParameterEntryValidator PEV;
  TEST_FOR_EXCEPT(0==paramList);
  paramList->set(
    paramName, defaultValue, docString,
    rcp_implicit_cast<const PEV>(
      stringToIntegralParameterEntryValidator<IntegralType>(
        strings, paramName
        )
      )
    );
}


template<class IntegralType>
void Teuchos::setStringToIntegralParameter(
  std::string const& paramName,
  std::string const& defaultValue,
  std::string const& docString,
  ArrayView<const std::string> const& strings,
  ArrayView<const IntegralType> const& integralValues, 
  ParameterList * paramList
  )
{
  typedef ParameterEntryValidator PEV;
  TEST_FOR_EXCEPT(0==paramList);
  paramList->set(
    paramName, defaultValue, docString,
    rcp_implicit_cast<const PEV>(
      stringToIntegralParameterEntryValidator<IntegralType>(
        strings, integralValues, paramName
        )
      )
    );
}


template<class IntegralType>
void Teuchos::setStringToIntegralParameter(
  std::string const& paramName,
  std::string const& defaultValue,
  std::string const& docString,
  ArrayView<const std::string> const& strings,
  ArrayView<const std::string> const& stringsDocs,
  ArrayView<const IntegralType> const& integralValues, 
  ParameterList * paramList
  )

{
  typedef ParameterEntryValidator PEV;
  TEST_FOR_EXCEPT(0==paramList);
  paramList->set(
    paramName, defaultValue, docString,
    rcp_implicit_cast<const PEV>(
      stringToIntegralParameterEntryValidator<IntegralType>(
        strings, stringsDocs, integralValues, paramName
        )
      )
    );
}


template<class IntegralType>
IntegralType Teuchos::getIntegralValue(
  ParameterList const& paramList, std::string const& paramName
  )
{
  const ParameterEntry &entry = paramList.getEntry(paramName);
  RCP<const StringToIntegralParameterEntryValidator<IntegralType> >
    integralValidator = getStringToIntegralParameterEntryValidator<IntegralType>(
      entry, paramList, paramName
      );
  return integralValidator->getIntegralValue(
    entry, paramName, paramList.name(), true
    );
}


template<class IntegralType>
std::string Teuchos::getStringValue(
  ParameterList const& paramList, std::string const& paramName
  )
{
  const ParameterEntry &entry = paramList.getEntry(paramName);
  RCP<const StringToIntegralParameterEntryValidator<IntegralType> >
    integralValidator = getStringToIntegralParameterEntryValidator<IntegralType>(
      entry, paramList, paramName
      );
  return integralValidator->getStringValue(
    entry, paramName, paramList.name(), true
    );
}


template<class IntegralType>
Teuchos::RCP<const Teuchos::StringToIntegralParameterEntryValidator<IntegralType> >
Teuchos::getStringToIntegralParameterEntryValidator(
  ParameterEntry const& entry, ParameterList const& paramList,
  std::string const& paramName
  )
{
  RCP<const ParameterEntryValidator>
    validator = entry.validator();
  TEST_FOR_EXCEPTION_PURE_MSG(
    is_null(validator), Exceptions::InvalidParameterType,
    "Error!  The parameter \""<<paramName<<"\" exists\n"
    "in the parameter (sub)list \""<<paramList.name()<<"\"\n"
    "but it does not contain any validator needed to extract\n"
    "an integral value of type \""<<TypeNameTraits<IntegralType>::name()<<"\"!"
    );
  RCP<const StringToIntegralParameterEntryValidator<IntegralType> >
    integralValidator
    =
    rcp_dynamic_cast<const StringToIntegralParameterEntryValidator<IntegralType> >(
      validator
      );
  TEST_FOR_EXCEPTION_PURE_MSG(
    is_null(integralValidator), Exceptions::InvalidParameterType,
    "Error!  The parameter \""<<paramName<<"\" exists\n"
    "in the parameter (sub)list \""<<paramList.name()<<"\"\n"
    "but it contains the wrong type of validator.  The expected validator type\n"
    "is \""<<TypeNameTraits<StringToIntegralParameterEntryValidator<IntegralType> >::name()<<"\"\n"
    "but the contained validator type is \""<<typeName(*validator)<<"\"!"
    );
  return integralValidator;
}


#endif // TEUCHOS_STANDARD_PARAMETER_ENTRY_VALIDATORS_H
