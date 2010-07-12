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
#include "Teuchos_Assert.hpp"
#include "Teuchos_StrUtils.hpp"
#include "Teuchos_TypeNameTraits.hpp"
#include <sys/stat.h>


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
   *             [in] Array of unique std::string names.
   * \param  defaultParameterName
   *             [in] The default name of the parameter (used in error messages)
   */
  StringToIntegralParameterEntryValidator(
    ArrayView<const std::string> const& strings,
    std::string const& defaultParameterName
    );

  /** \brief Construct with a mapping from strings to aribitrary typed
   * integral values.
   *
   * \param  strings
   *             [in] Array of unique std::string names.
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
    ArrayView<const std::string> const& strings,
    ArrayView<const IntegralType> const& integralValues, 
    std::string const& defaultParameterName
    );

  /** \brief Construct with a mapping from strings (with documentation) to
   * aribitrary typed integral values.
   *
   * \param  strings
   *             [in] Array of unique std::string names.
   * \param  stringsDocs
   *             [in] Array of documentation strings for each std::string value.
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
   * \param  str  [in] String that is being used to lookup the corresponding
   *              integral value.
   * \param  paramName
   *              [in] Optional name that will be used to generate error messages.
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
   * \param  entry
   *              [in] The std::string entry.
   * \param  paramName
   *              [in] Optional name that will be used to generate error messages.
   * \param  sublistName
   *              [in] The name of the sublist.
   * \param  activeQuery
   *              [in] If true, then this lookup will be recored as an active query
   *              which will turn the <tt>isUsed</tt> bool to <tt>true</tt>.
   */
  IntegralType getIntegralValue(
    const ParameterEntry &entry, const std::string &paramName = "",
    const std::string &sublistName = "", const bool activeQuery = true
    ) const;

  /** \brief Get and validate a std::string value embedded in a
   * <tt>ParameterEntry</tt> object.
   *
   *
   * \param  entry
   *              [in] The std::string entry.
   * \param  paramName
   *              [in] Optional name that will be used to generate error messages.
   * \param  sublistName
   *              [in] The name of the sublist.
   * \param  activeQuery
   *              [in] If true, then this lookup will be recored as an active query
   *              which will turn the <tt>isUsed</tt> bool to <tt>true</tt>.
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

  /** \brief Get a pointer to the array containing all the documentation strings.
   *
   * \return A point to the array containing all the documentation strings.
   */
  ValidStringsList getStringDocs() const;

  /** \brief Get the name of the default parameter for the validator.
   *
   * \return The name of the default parameter for the validator.
   */
  const std::string& getDefaultParameterName() const;

  /** \brief Validate the std::string and pass it on..
   *
   * \param  str  [in] String that is being used to lookup the corresponding
   *              integral value.
   * \param  name [in] Optional name that will be used to generate error messages.
   *
   * If the std::string name <tt>str</tt> does not exist, the an std::exception will be
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
  const std::string getXMLTagName() const;

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
  static const std::string& tagName(){
  	static const std::string tagName_ = TypeNameTraits<IntegralType>::name() + "stringtointegralvalidator";
	return tagName_;
  }
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


/** \brief Set up a std::string parameter that will use an embedded validator to
 * allow the extraction of an integral value.
 *
 * The function <tt>getIntegralValue()</tt> can then be used to extract the
 * integral value of the std::string parameter.  In this case, the integral value
 * return will just be the zero-based index of the std::string value in the list
 * <tt>strings</tt>.
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


/** \brief Set up a std::string parameter that will use an embedded validator to
 * allow the extraction of an integral value from a list of integral values.
 *
 * The function <tt>getIntegralValue()</tt> can then be used to extract the
 * integral value of the std::string parameter.  In this case, the integral value
 * return will just be the zero-based index of the std::string value in the list
 * <tt>strings</tt>.
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


/** \brief Set up a std::string parameter with documentation strings for each valid
 * value that will use an embedded validator to allow the extraction of an
 * integral value from a list of integral values.
 *
 * The function <tt>getIntegralValue()</tt> can then be used to extract the
 * integral value of the std::string parameter.  In this case, the integral value
 * return will just be the zero-based index of the std::string value in the list
 * <tt>strings</tt>.
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


/** \brief Get an integral value for a parameter that is assumed to already be set.
 *
 * This function does a dynamic cast to get the underlying valiator of type
 * StringToIntegralParameterEntryValidator<IntegralType>.  If this dynamic
 * cast failes then an <tt>Exceptions::InvalidParameterType</tt> std::exception is
 * thrown with an excellent error message.
 *
 * \relates ParameterList
 */
template<class IntegralType>
IntegralType getIntegralValue(
  ParameterList const& paramList, std::string const& paramName
  );


/** \brief Get a std::string value for a parameter that is assumed to already be set.
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


/** \brief Get a StringToIntegralParameterEntryValidator<IntegralType> object out of
 * a ParameterEntry object.
 *
 * This function with thrown of the validator does not exist.
 */
template<class IntegralType>
RCP<const StringToIntegralParameterEntryValidator<IntegralType> >
getStringToIntegralParameterEntryValidator(
  ParameterEntry const& entry, ParameterList const& paramList,
  std::string const& paramName
  );


/** \brief Return the std::string name of the verbosity level as it is accepted by the
 * verbosity level parameter.
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


/** \brief Standard implementation of a ParameterEntryValidator that accepts
 * numbers from a number of different formats and converts them to numbers in
 * another format.
 *
 * Objects of this type are meant to be used as both abstract objects passed
 * to <tt>Teuchos::ParameterList</tt> objects to be used to validate parameter
 * types and values, and to be used by the code that reads parameter values.
 * Having a single definition for the types of valids input and outputs for a
 * parameter value makes it easier to write error-free validated code.
 */
class TEUCHOS_LIB_DLL_EXPORT AnyNumberParameterEntryValidator : public ParameterEntryValidator {
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

  /** \name Constructors */
  //@{

  /** \brief Construct with a preferrded type of double and accept all
   * types.
   */
  AnyNumberParameterEntryValidator();

  /** \brief Construct with allowed input and output types and the preferred
   * type.
   *
   * \param preferredType
   *          [in] Determines the preferred type.  This enum value is used to 
   *          set the default value in the override <tt>validateAndModify()</tt>.
   * \param acceptedType
   *          [in] Determines the types that are allowed in the parameter list.
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

  bool allowDouble() const;

  bool allowInt() const;
  
  bool allowString() const;

  EPreferredType prefferedType() const;

  static const std::string& getPrefferedTypeString(EPreferredType enumValue){
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
			throw std::runtime_error("Cannot convert enumValue: " + toString(enumValue) + " to a string");
	}
	//Should never get here. This code is here so that a warning is not generated.
	static const std::string& emptyString("");
	return emptyString;
  }

  static EPreferredType getPrefferedTypeStringEnum(const std::string& enumString){
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
  const std::string getXMLTagName() const;

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

  static const std::string& getIntEnumString(){
  	static const std::string intEnumString_ = TypeNameTraits<int>::name();
	return intEnumString_;
  }

  static const std::string& getDoubleEnumString(){
  	static const std::string doubleEnumString_ = TypeNameTraits<double>::name();
	return doubleEnumString_;
  }

  static const std::string& getStringEnumString(){
  	static const std::string stringEnumString_ = TypeNameTraits<std::string>::name();
	return stringEnumString_;
  }

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
 * If the underlying parameter type is already a std::string, then all is good.
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
TEUCHOS_LIB_DLL_EXPORT std::string getNumericStringParameter(
  ParameterList const& paramList,
  std::string const& paramName
  );

template <class T>
struct UndefinedEnhancedNumberTraits{
	static inline T notDefined() { return T::this_type_is_missing_a_specialization(); }
};

template <class T>
class EnhancedNumberTraits{
public:
	static inline T min() { return UndefinedEnhancedNumberTraits<T>::notDefined(); }
	static inline T max() { return UndefinedEnhancedNumberTraits<T>::notDefined(); }
	static inline T defaultStep() { return UndefinedEnhancedNumberTraits<T>::notDefined(); }
	static inline unsigned short defaultPrecision() { return UndefinedEnhancedNumberTraits<T>::notDefined(); }
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
	static inline unsigned short defaultPrecision() { return 2; }
};

template<>
class EnhancedNumberTraits<float>{
public:
	static inline float min() { return -std::numeric_limits<float>::max(); }
	static inline float max() { return std::numeric_limits<float>::max(); }
	static inline float defaultStep() { return 1; }
	static inline unsigned short defaultPrecision() { return 2; }
};
 
/**
 * A Template base class for NumberValidators.
 */
template <class T>
class EnhancedNumberValidator : public ParameterEntryValidator{
public:
	/**
	 * Constructs a EnhancedNumberValidator.
	 *
	 * @param min The minimum acceptable value for this validator.
	 * @param max The maximum acceptable value for this validator.
	 * @param step The increments at which the values being validated should be changed.
	 * @param precision The number of decimials places to which the values validated shold be compared to the
	 * min and max. This parameter is pretty much meamingless for non-floating point types.
	 */
	EnhancedNumberValidator(T min, T max, T step=EnhancedNumberTraits<T>::defaultStep(), unsigned short precision=EnhancedNumberTraits<T>::defaultPrecision()):ParameterEntryValidator(),
	minVal(min), maxVal(max), step_(step), precision_(precision), containsMin(true), containsMax(true){}

	/**
	 * Constructs a EnhancedNumberValidator without an explicit minimum or maximum.
	 *
	 * @param step The increments at which the values being validated should be changed.
	 * @param precision The number of decimials places to which the values validated shold be compared to the
	 * min and max. This parameter is pretty much meamingless for non-floating point types.
	 */
	EnhancedNumberValidator():
		ParameterEntryValidator(),
		minVal(EnhancedNumberTraits<T>::min()),
		maxVal(EnhancedNumberTraits<T>::max()),
		step_(EnhancedNumberTraits<T>::defaultStep()),
		precision_(EnhancedNumberTraits<T>::defaultPrecision()),
		containsMin(false),
		containsMax(false){}
		
	/**
	 * Sets the minimum acceptable value for the validator.
	 * 
	 * @param min The desired minimum acceptable value for the validator.
	 */
	void setMin(T min){
		minVal = min;
		containsMin = true;
	}

	/**
	 * Sets the maximum acceptable value for the validator.
	 * 
	 * @param min The desired maximum acceptable value for the validator.
	 */
	void setMax(T max){
		maxVal = max;
		containsMax = true;
	}

	/**
	 * Gets the minimum acceptable value for the validator.
	 *
	 *@return The minimum acceptable value for the validator.
	 */
	T getMin() const{
		return minVal;
	}

	/**
	 * Gets the maximum acceptable value for the validator.
	 *
	 *@return The maximum acceptable value for the validator.
	 */
	T getMax() const{
		return maxVal;
	}

	/**
	 * Determines whether or not the validator has a minimum value.
	 *
	 * @return True if the validator has a minimum value, false otherwise.
	 */
	bool hasMin() const{
		return containsMin;
	}

	/**
	 * Determines whether or not the validator has a maximum value.
	 *
	 * @return True if the validator has a maximum value, false otherwise.
	 */ 
	bool hasMax() const{
		return containsMax;
	}

	/**
	 * Gets the step being used for the validator.
	 *
	 * @return The step being used for the validator.
	 */
	T getStep() const{
		return step_;
	}

	/**
	 * Sets the step being used for the validator.
	 *
	 * @param The step to be used for the validator.
	 */
	void setStep(T step){
		step_ = step;
	}

	/**
	 * Sets the precision specified for the validator.
	 *
	 * @param The precision specific for the validator.
	 */
	void setPrecision(unsigned short precision){
		precision_ = precision;
	}

	/**
	 * Gets the precision specified for the validator.
	 *
	 * @return The precision specific for the validator.
	 */
	unsigned short getPrecision() const{
		return precision_;
	}

	ValidStringsList validStringValues() const{
		return null;
	}

	void validate(ParameterEntry const &entry, std::string const &paramName, std::string const &sublistName) const{
		any anyValue = entry.getAny(true);
		if(anyValue.type() == typeid(T) ){
			bool isValueInRange = false;
			if(precision_ != 0){
				T precisionPadding = pow((T)10,(-((T)precision_)));
				any_cast<T>(anyValue) >= minVal-((T)precisionPadding) && any_cast<T>(anyValue) <= maxVal+((T)precisionPadding) ?
				isValueInRange = true : isValueInRange=false;
			}
			else{
				any_cast<T>(anyValue) >= minVal && any_cast<T>(anyValue) <= maxVal ?
				isValueInRange = true : isValueInRange=false;
			}
			if(!(isValueInRange)){
				std::stringstream oss;
				std::string msg;
				oss << "Aww shoot! Sorry bud, but it looks like the \"" << paramName << "\"" <<
				" parameter in the \"" << sublistName << "\" sublist didn't quite work out.\n" <<
				"No need to fret though. I'm sure it's just a small mistake. Maybe the information below "<<
				"can help you figure out what went wrong.\n\n"
				"Error: The value that was entered doesn't fall with in " <<
				"the range set by the validator.\n" <<
				"Parameter: " << paramName << "\n" <<
				"Min: " << minVal << "\n" <<
				"Max: " << maxVal << "\n" <<
				"Value entered: " << (any_cast<T>(anyValue)) << "\n";
				msg = oss.str();
				throw Exceptions::InvalidParameterValue(msg);
			}	
		}
		else{
			const std::string &entryName = entry.getAny(false).typeName();
			std::stringstream oss;
			std::string msg;
			oss << "Aww shoot! Sorry bud, but it looks like the \"" << paramName << "\"" <<
			" parameter in the \"" << sublistName << "\" sublist didn't quite work out.\n" <<
			"No need to fret though. I'm sure it's just a small mistake. Maybe the information below "<<
			"can help you figure out what went wrong.\n\n"
			"Error: The value that you entered was the wrong type.\n" <<
			"Parameter: " << paramName << "\n" <<
			"Type specified: " << entryName << "\n" <<
			"Type accepted: " << Teuchos::TypeNameTraits<T>::name() << "\n";
			msg = oss.str();
			throw Exceptions::InvalidParameterType(msg);
		}
	}

	/** \brief . */
	const std::string getXMLTagName() const{
		return TypeNameTraits<T>::name() + "enhancednumbervalidator";
	}

	void printDoc(std::string const &docString, std::ostream &out) const{
		StrUtils::printLines(out,"# ",docString);
		out << "#  Validator Used: \n";
		out << "#  	Number Validator\n";
		out << "#  	Type: " << Teuchos::TypeNameTraits<T>::name() << "\n";
		out << "#  	Min (inclusive): " << minVal << "\n";
		out << "#  	Max (inclusive): " << maxVal << "\n";
	}

private:
	/**
	 * The minimum value accepted by the validator.
	 */
	T minVal;

	/**
	 * The maximum value accepted by the validator.
	 */
	T maxVal;

	/**
	 * The increment to use when increaseing or decreaseing the value the validator is validating.
	 */
	T step_;

	/**
	 * The amount of decimal to which the value to be validated will be compared to the
	 * maximum and minimum. A precision of 0 means the the value to be validated
	 * must fall exactly within the mininmum and maximum. Notd this value is pretty much
	 * meaningingless for non-floating point value types.
	 */
	unsigned short precision_;

	/**
	 * Whether or not a minimum value has been specified for this validator.
	 */
	bool containsMin;

	/**
	 * Whetehr or not a maximum value has been specified for this validator.
	 */
	bool containsMax;
};


/**
 * Simply indicates that the parameter entry with this validator should
 * contain a filename.
 */
class FileNameValidator : public ParameterEntryValidator{
public:
	static const bool mustAlreadyExistDefault=false;
	/**
	 * Constructs a FileNameValidator.
	 *
	 * @param mustAlreadyExist True if the file the user specifies should already exists, false otherwise.
	 */
	FileNameValidator(bool mustAlreadyExist=mustAlreadyExistDefault);

	/**
	 * Gets the variable describing whether or not this validator wants the file that is specified to
	 * already exist.
	 *
	 * @return Whether or not the validator requires the file to already exist
	 */
	bool fileMustExist() const;

	/**
	 * Sets whether or not the validator requires the file to already exist.
	 *
	 * @param shouldFileExist True if the file should already exist, false otherwise.
	 * @return The new value of the shouldFileExist variable.
	 */
	bool setFileMustExist(bool shouldFileExist);

	ValidStringsList validStringValues() const;

	void validate(ParameterEntry const &entry, std::string const &paramName, std::string const &sublistName) const;

	/** \brief . */
	const std::string getXMLTagName() const;

	void printDoc(std::string const &docString, std::ostream &out) const;

private:
	/**
	 * Whether or not the file specified in the parameter should already exist.
	 */
	bool mustAlreadyExist_;
};

/**
 * A simple validator that only allows certain string values to be choosen or simply
 * enforces that a particular parameter have a std::string for a value.
 */
class StringValidator : public ParameterEntryValidator{
public:
	/**
	 * Constructs a StringValidator.
	 *
	 * @param validStrings A list of valid string values for this validator.
	 */
	StringValidator(const Teuchos::Array<std::string> &validStrings);

	/**
	 * Constructs a StringValidator.
	 */
	StringValidator();

	/**
	 * Sets the Array of valid strings and returns what the current array of valid
	 * string now is.
	 *
	 * @param validStrings What the array for the valid strings should contain.
	 * @return What the arry for the valid strings now conatians.
	 */
	ValidStringsList setValidStrings(const Teuchos::Array<std::string> &validStrings);

	ValidStringsList validStringValues() const;

	void validate(ParameterEntry const &entry, std::string const &paramName, std::string const &sublistName) const;

	/** \brief . */
	const std::string getXMLTagName() const;

	void printDoc(std::string const &docString, std::ostream &out) const;

	//virtual XMLObject getXML() const;

  /*static const std::string& getTagName(){
	static const std::string tagName = "stringvalidator";
	return tagName;
  }*/
private:
	/**
	 * An array containing a list of all the valid string values.
	 */
	ValidStringsList validStrings_;
};

/**
 * An Abstract base class for all ArrayValidators
 */
template<class ValidatorType, class EntryType>
class ArrayValidator : public ParameterEntryValidator{
public:
	/**
	 * Constructs a ArrayValidator.
	 *
	 * @param prototypeValidator The validator to be used on each
	 * entry in the array.
	 */
	ArrayValidator(RCP<ValidatorType> prototypeValidator):
		ParameterEntryValidator(),
	    prototypeValidator_(prototypeValidator){}

	ValidStringsList validStringValues() const{
		return prototypeValidator_->validStringValues();
	}

	virtual void validate(ParameterEntry const &entry, std::string const &paramName, std::string const &sublistName) const{
		any anyValue = entry.getAny(true);
		if(anyValue.type() == typeid(Array<EntryType>)){
			Array<EntryType> extracted = getValue<Teuchos::Array<EntryType> >(entry);
			for(int i = 0; i<extracted.size(); ++i){
				ParameterEntry dummyParameter;
				dummyParameter.setValue(extracted[i]);
				try{
					prototypeValidator_->validate(dummyParameter, paramName, sublistName);
				}
				catch(Exceptions::InvalidParameterValue& e){
					std::stringstream oss;
					oss << "Array Validator Exception: \n" << "Bad Index: " << i << "\n" << e.what();
					throw Exceptions::InvalidParameterValue(oss.str());
				}
			}
		}
		else{
			const std::string &entryName = entry.getAny(false).typeName();
			std::stringstream oss;
			std::string msg;
			oss << "Aww shoot! Sorry bud, but it looks like the \"" << paramName << "\"" <<
			" parameter in the \"" << sublistName << "\" sublist didn't quite work out.\n" <<
			"No need to fret though. I'm sure it's just a small mistake. Maybe the information below "<<
			"can help you figure out what went wrong.\n\n"
			"Error: The value you entered was the wrong type.\n" <<
			"Parameter: " << paramName << "\n" <<
			"Type specified: " << entryName << "\n" <<
			"Type accepted: " << TypeNameTraits<Array<std::string> >::name() << "\n";
			msg = oss.str();
			throw Exceptions::InvalidParameterType(msg);
		}
	}

	/** \brief . */
	const std::string getXMLTagName() const{
		return "arrayvalidator";
	}

	virtual void printDoc(std::string const &docString, std::ostream &out) const{
		StrUtils::printLines(out,"# ",docString);
		std::string toPrint;
		toPrint += "ArrayValidator:\n";
		toPrint += "Prototype Validator:\n";
		prototypeValidator_->printDoc(toPrint, out);
	}

	RCP<const ValidatorType> getPrototype() const{
		return prototypeValidator_;
	}

private:
	/**
	 * The prototype validator to be applied to each entry in the Array.
	 */
	RCP<ValidatorType> prototypeValidator_;
};

class ArrayStringValidator : public ArrayValidator<StringValidator, std::string>{
public:
	ArrayStringValidator(RCP<StringValidator> prototypeValidator):
		ArrayValidator<StringValidator, std::string>(prototypeValidator){}
};

class ArrayFileNameValidator : public ArrayValidator<FileNameValidator, std::string>{
public:
	ArrayFileNameValidator(RCP<FileNameValidator> prototypeValidator):
		ArrayValidator<FileNameValidator, std::string>(prototypeValidator){}
};

template<class T>
class ArrayNumberValidator : public ArrayValidator<EnhancedNumberValidator<T>, T>{
public:
	ArrayNumberValidator(RCP<EnhancedNumberValidator<T> > prototypeValidator):
		ArrayValidator<EnhancedNumberValidator<T>, T>(prototypeValidator){}
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
  )
  :defaultParameterName_(defaultParameterName)
{
  typedef typename map_t::value_type val_t;
  for( int i = 0; i < static_cast<int>(strings.size()); ++i ) {
    const bool unique = map_.insert( val_t( strings[i], i ) ).second;
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
  )
  :defaultParameterName_(defaultParameterName)
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
  )
  :defaultParameterName_(defaultParameterName)
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
const std::string StringToIntegralParameterEntryValidator<IntegralType>::getXMLTagName() const{
	return TypeNameTraits<IntegralType>::name() + "stringintegralvalidator";
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

/*template<class IntegralType>
XMLObject StringToIntegralParameterEntryValidator<IntegralType>::getXML() const
{
	XMLObject valiTag(getTagName());
	XMLObject stringsTag("strings");
	XMLObject stringsDocsTag("stringdocs");
	XMLObject integralValuesTag("integralvalues");
	for(typename map_t::const_iterator it = map_.begin(); it != map_.end(); ++it){
		XMLObject stringTag("string");
		stringsTag.addAttribute("value", it->first);
		stringsTag.addChild(stringTag);
		XMLObject integralValueTag("integralvalue");
		integralValueTag.addAttribute("value",toString(it->second));
		integralValuesTag.addChild(integralValueTag);
	}
	if(!validStringValuesDocs_.is_null()){
		for(Array<std::string>::const_iterator it = validStringValuesDocs_->begin(); it != validStringValuesDocs_->end(); ++it){
			XMLObject stringDocTag("stringdoc");
			stringDocTag.addAttribute("value", *it);
		}
	}
	XMLObject defaultParameterNameTag("defaultparametername");
	defaultParameterNameTag.addAttribute("value", defaultParameterName_);

	valiTag.addChild(stringsTag);
	valiTag.addChild(stringsDocsTag);
	valiTag.addChild(integralValuesTag);
	valiTag.addChild(defaultParameterNameTag);
	return valiTag;
}*/


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
  typename map_t::const_iterator itr = map_.begin();
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
