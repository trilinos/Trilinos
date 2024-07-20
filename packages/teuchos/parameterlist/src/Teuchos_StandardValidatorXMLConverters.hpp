// @HEADER
// *****************************************************************************
//                    Teuchos: Common Tools Package
//
// Copyright 2004 NTESS and the Teuchos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef TEUCHOS_STANDARDVALIDATORXMLCONVERTERS_HPP
#define TEUCHOS_STANDARDVALIDATORXMLCONVERTERS_HPP

/*! \file Teuchos_StandardValidatorXMLConverters.hpp
    \brief A collection of standard ValidatorXMLConverters.
*/

#include "Teuchos_ValidatorXMLConverter.hpp"
#include "Teuchos_StandardParameterEntryValidators.hpp"
#include "Teuchos_ValidatorXMLConverterDB.hpp"
#include "Teuchos_XMLParameterListReader.hpp"
#include "Teuchos_DummyObjectGetter.hpp"


namespace Teuchos {


/** \brief Convert a StringToIntegralParameterEntryValidator to and from XML.

This class knows how to write a
StringToIntegralParameterEntryValidator to XML, and create an
StringToIntegralParameterEntryValidator from its XML representation.

Here is the XML representation of a StringToIntegralValidator:

\code
 <Validator
   type="StringToIntegralValidator(NumberType)"
   defaultParameterName="Name of default parameter"
   caseSensitive="[true|false]"
   validatorId="Validator id"
 >
   <String
     stringValue="Value 1"
     integralValue="int value 1"
     stringDoc="Documentation for Value 1"
   />
   <String
     stringValue="Value 2"
     integralValue="int value 2"
     stringDoc="Documentation for Value 2"
   />
   ...More String Values...
 </Validator>
\endcode

Here, "[true|false]" means the XML string representation of either
Boolean true, or Boolean false.

The "integralValue", "stringDoc", and "caseSensitive" XML attributes
are optional. However, if one of the "String" tags includes an
"integralValue" and/or a "stringDoc" XML attribute, all other "String"
tags must do so as well.

 */
template<class IntegralType>
class StringToIntegralValidatorXMLConverter :
  public ValidatorXMLConverter
{

public:

  /** \name Overridden from ValidatorXMLConverter */
  //@{

  /** \brief . */
  RCP<ParameterEntryValidator> convertXML(
    const XMLObject& xmlObj,
    const IDtoValidatorMap& validatorIDsMap) const;

  /** \brief . */
  void convertValidator(
    const RCP<const ParameterEntryValidator> validator,
    XMLObject& xmlObj,
    const ValidatortoIDMap& validatorIDsMap) const;

  #ifdef HAVE_TEUCHOS_DEBUG
  /** \brief . */
  RCP<const ParameterEntryValidator >
  getDummyValidator() const{
    return DummyObjectGetter<
    StringToIntegralParameterEntryValidator<IntegralType> >::getDummyObject();
  }
  #endif

  //@}

private:

  /** \name Private Members */
  //@{

  /** \brief . */
  static const std::string& getIntegralValueAttributeName() {
    static const std::string integralValueAttributeName_ = "integralValue";
    return integralValueAttributeName_;
  }

  /** \brief . */
  static const std::string& getStringTagName() {
    static const std::string stringTagName_ = "String";
    return stringTagName_;
  }

  /** \brief . */
  static const std::string& getStringValueAttributeName() {
    static const std::string stringValueAttributeName_ = "stringValue";
    return stringValueAttributeName_;
  }

  /** \brief . */
  static const std::string& getStringDocAttributeName() {
    static const std::string stringDocAttributeName_ = "stringDoc";
    return stringDocAttributeName_;
  }

  //! Name (tag) of the default parameter attribute.
  static const std::string& getDefaultParameterAttributeName() {
    static const std::string defaultParameterAttributeName_ =
      "defaultParameterName";
    return defaultParameterAttributeName_;
  }

  //! Name (tag) of the caseSensitive attribute.
  static const std::string& getCaseSensitiveAttributeName() {
    static const std::string caseSensitiveAttributeName_ =
      "caseSensitive";
    return caseSensitiveAttributeName_;
  }

  //@}

};


//
// Implementations
//


template<class IntegralType>
RCP<ParameterEntryValidator>
StringToIntegralValidatorXMLConverter<IntegralType>::convertXML(
  const XMLObject& xmlObj,
  const IDtoValidatorMap& /*validatorIDsMap*/) const
{
  Array<std::string> strings;
  Array<std::string> stringDocs;
  Array<IntegralType> integralValues;
  for (int i=0; i<xmlObj.numChildren(); ++i) {
    XMLObject currentChild = xmlObj.getChild(i);
    TEUCHOS_TEST_FOR_EXCEPTION(currentChild.getTag() != getStringTagName(),
      BadTagException,
      "Error converting xmlObject to "
      "StringToIntegralParameterEntryValidator." << std::endl <<
      "Unrecognized tag: " << currentChild.getTag());
    strings.append(currentChild.getRequired(getStringValueAttributeName()));
    if (currentChild.hasAttribute(getIntegralValueAttributeName())) {
      integralValues.append(
        currentChild.getRequired<IntegralType>(
          getIntegralValueAttributeName()));
    }
    if (currentChild.hasAttribute(getStringDocAttributeName())) {
      stringDocs.append(
        currentChild.getRequired<std::string>(getStringDocAttributeName()));
    }
  }
  std::string defaultParameterName =
    xmlObj.getRequired(getDefaultParameterAttributeName());

  // The "caseSensitive" attribute is not required.  It is true by default.
  const bool caseSensitive =
    xmlObj.getWithDefault<bool> (getCaseSensitiveAttributeName (), true);

  typedef StringToIntegralParameterEntryValidator<IntegralType> ret_type;
  if (stringDocs.size() != 0 && integralValues.size() != 0) {
    return rcp (new ret_type (strings, stringDocs, integralValues (), defaultParameterName, caseSensitive));
  }
  else if (integralValues.size() != 0) {
    return rcp (new ret_type (strings, integralValues(), defaultParameterName, caseSensitive));
  }
  else {
    return rcp (new ret_type (strings, defaultParameterName, caseSensitive));
  }
}


template<class IntegralType>
void StringToIntegralValidatorXMLConverter<IntegralType>::convertValidator(
  const RCP<const ParameterEntryValidator> validator,
  XMLObject& xmlObj,
  const ValidatortoIDMap& /*validatorIDsMap*/) const
{
  RCP<const StringToIntegralParameterEntryValidator<IntegralType> >
    castedValidator =
    rcp_dynamic_cast<
      const StringToIntegralParameterEntryValidator<IntegralType> >(
        validator, true);

  RCP<const Array<std::string> > stringValues =
    castedValidator->validStringValues();
  RCP<const Array<std::string> > stringDocValues =
    castedValidator->getStringDocs();

  bool hasStringDocs =
    !(stringDocValues.is_null()) && (stringDocValues->size() != 0);
  for (int i =0; i<stringValues->size(); ++i) {
    XMLObject stringTag(getStringTagName());
    stringTag.addAttribute(getStringValueAttributeName(), (*stringValues)[i]);
    stringTag.addAttribute(getIntegralValueAttributeName(),
      castedValidator->getIntegralValue((*stringValues)[i]));
    if (hasStringDocs) {
      stringTag.addAttribute(
        getStringDocAttributeName(), (*stringDocValues)[i]);
    }
    xmlObj.addChild(stringTag);
  }
  xmlObj.addAttribute(getDefaultParameterAttributeName(),
    castedValidator->getDefaultParameterName());

  // Add "caseSensitive" bool attribute here.
  const bool caseSensitive = castedValidator->isCaseSensitive ();
  xmlObj.addBool (getCaseSensitiveAttributeName (), caseSensitive);

  xmlObj.addAttribute(getIntegralValueAttributeName(),
    TypeNameTraits<IntegralType>::name());
}

/**
 * \brief Converts BoolParameterEntryValidators to and from XML.
 *
 * The valid XML representation for an BoolParameterEntryValidators is:
 * \code
  <Validator type="BoolValidator"
  />
  \endcode
 */
class TEUCHOSPARAMETERLIST_LIB_DLL_EXPORT BoolValidatorXMLConverter : public ValidatorXMLConverter
{

public:

  /** \name Overridden from ValidatorXMLConverter */
  //@{

  /** \brief . */
  RCP<ParameterEntryValidator> convertXML(
    const XMLObject& xmlObj,
    const IDtoValidatorMap& validatorIDsMap) const;

  /** \brief . */
  void convertValidator(
    const RCP<const ParameterEntryValidator> validator,
    XMLObject& xmlObj,
    const ValidatortoIDMap& validatorIDsMap) const;

  #ifdef HAVE_TEUCHOS_DEBUG
  /** \brief . */
  RCP<const ParameterEntryValidator> getDummyValidator() const;
  #endif

  //@}

private:

  /** \name Private Members */
  //@{

  // currently empty

  //@}

};

/**
 * \brief Converts AnyNumberParameterEntryValidators to and from XML.
 *
 * The valid XML representation for an AnyNumberParameterEntryValidator is:
 * \code
  <Validator type="AnyNumberValidator"
   allowInt="True or False"
   allowDouble="True or False"
   allowString="True or False"
   prefferedType="Prefered type"
  />
  \endcode
 */
class TEUCHOSPARAMETERLIST_LIB_DLL_EXPORT AnyNumberValidatorXMLConverter : public ValidatorXMLConverter
{

public:

  /** \name Overridden from ValidatorXMLConverter */
  //@{

  /** \brief . */
  RCP<ParameterEntryValidator> convertXML(
    const XMLObject& xmlObj,
    const IDtoValidatorMap& validatorIDsMap) const;

  /** \brief . */
  void convertValidator(
    const RCP<const ParameterEntryValidator> validator,
    XMLObject& xmlObj,
    const ValidatortoIDMap& validatorIDsMap) const;

  #ifdef HAVE_TEUCHOS_DEBUG
  /** \brief . */
  RCP<const ParameterEntryValidator> getDummyValidator() const;
  #endif

  //@}

private:

  /** \name Private Members */
  //@{

  /** \brief . */
  static const std::string& getAllowIntAttributeName() {
    static const std::string allowIntAttributeName_ = "allowInt";
    return allowIntAttributeName_;
  }

  /** \brief . */
  static const std::string& getAllowDoubleAttributeName() {
    static const std::string allowDoubleAttributeName_ = "allowDouble";
    return allowDoubleAttributeName_;
  }

  /** \brief . */
  static const std::string& getAllowStringAttributeName() {
    static const std::string allowStringAttributeName_ = "allowString";
    return allowStringAttributeName_;
  }

  /** \brief . */
  static const std::string& getPrefferedTypeAttributeName() {
    static const std::string prefferedTypeAttributeName_ = "prefferedType";
    return prefferedTypeAttributeName_;
  }

  //@}

};


/** \brief Converts EnhancedNumberValidators to and from XML.
 *
 * The valid XML representation of an EnhancedNumberValidator is:
 * \code
  <Validator type="EnhancedNumberValidator(numbertype)"
   min="Minimum Value"
   max="Maximum Value"
   step="Step Value"
   precision="Precision Value"
   validatorId="Validator Id"
  />
  \endcode
 * The "min", "max", "step", and "precision" XML attributes are all optional.
 */
template<class T>
class EnhancedNumberValidatorXMLConverter : public ValidatorXMLConverter
{

public:

  /** \name Overridden from ValidatorXMLConverter */
  //@{

  /** \brief . */
  RCP<ParameterEntryValidator> convertXML(
    const XMLObject& xmlObj,
    const IDtoValidatorMap& validatorIDsMap) const;

  /** \brief . */
  void convertValidator(
    const RCP<const ParameterEntryValidator> validator,
    XMLObject& xmlObj,
    const ValidatortoIDMap& validatorIDsMap) const;

#ifdef HAVE_TEUCHOS_DEBUG
  /** \brief . */
  RCP<const ParameterEntryValidator> getDummyValidator() const{
    return DummyObjectGetter<EnhancedNumberValidator<T> >::getDummyObject();
  }
#endif

  //@}

private:

  /** \name Private Members */
  //@{

  /** \brief . */
  static const std::string& getMinAttributeName() {
    static const std::string minAttributeName = "min";
    return minAttributeName;
  }

  /** \brief . */
  static const std::string& getMaxAttributeName() {
    static const std::string maxAttributeName = "max";
    return maxAttributeName;
  }

  /** \brief . */
  static const std::string& getStepAttributeName() {
    static const std::string stepAttributeName = "step";
    return stepAttributeName;
  }

  /** \brief . */
  static const std::string& getPrecisionAttributeName() {
    static const std::string precisionAttributeName = "precision";
    return precisionAttributeName;
  }

  //@}

};


template<class T>
RCP<ParameterEntryValidator>
EnhancedNumberValidatorXMLConverter<T>::convertXML(
  const XMLObject& xmlObj,
  const IDtoValidatorMap& /*validatorIDsMap*/) const
{
  RCP<EnhancedNumberValidator<T> > toReturn =
    rcp(new EnhancedNumberValidator<T>);
  T step = xmlObj.getWithDefault(
    getStepAttributeName(), EnhancedNumberTraits<T>::defaultStep());
  toReturn->setStep(step);
  unsigned short int precision = xmlObj.getWithDefault(
   getPrecisionAttributeName(),
   EnhancedNumberTraits<T>::defaultPrecision());
  toReturn->setPrecision(precision);
  if (xmlObj.hasAttribute(getMinAttributeName())) {
    toReturn->setMin(xmlObj.getRequired<T>(getMinAttributeName()));
  }
  if (xmlObj.hasAttribute(getMaxAttributeName())) {
    toReturn->setMax(xmlObj.getRequired<T>(getMaxAttributeName()));
  }
  return toReturn;
}


template<class T>
void EnhancedNumberValidatorXMLConverter<T>::convertValidator(
  const RCP<const ParameterEntryValidator > validator,
  XMLObject& xmlObj,
  const ValidatortoIDMap& /*validatorIDsMap*/) const
{
  RCP<const EnhancedNumberValidator<T> > castedValidator =
    rcp_dynamic_cast<const EnhancedNumberValidator<T> >(validator, true);
  if (castedValidator->hasMin()) {
    xmlObj.addAttribute<T>(getMinAttributeName(), castedValidator->getMin());
  }
  if (castedValidator->hasMax()) {
    xmlObj.addAttribute<T>(getMaxAttributeName(), castedValidator->getMax());
  }
  xmlObj.addAttribute<T>(getStepAttributeName(), castedValidator->getStep());
  xmlObj.addAttribute<short unsigned int>(
    getPrecisionAttributeName(), castedValidator->getPrecision());
}


/**
 * \brief Converts FileNameValidators to and from XML.
 *
 * The valid XML representation of a FileNameValidator is:
 *
 * \code
  <Validator type="FilenameValidator"
   fileMustExist="Bool Value"
   validatorId="Validator Id"
  />
  \endcode
 *
 * The "fileMustExist" XML attribute is optional.
 */
class TEUCHOSPARAMETERLIST_LIB_DLL_EXPORT FileNameValidatorXMLConverter : public ValidatorXMLConverter
{

public:

  /** \name Overridden from ValidatorXMLConverter */
  //@{

  /** \brief . */
  RCP<ParameterEntryValidator> convertXML(
    const XMLObject& xmlObj,
    const IDtoValidatorMap& validatorIDsMap) const;

  /** \brief . */
  void convertValidator(
    const RCP<const ParameterEntryValidator> validator,
    XMLObject& xmlObj,
    const ValidatortoIDMap& validatorIDsMap) const;

  #ifdef HAVE_TEUCHOS_DEBUG
  /** \brief . */
  RCP<const ParameterEntryValidator> getDummyValidator() const;
  #endif

  //@}

private:

  /** \name Private Members */
  //@{

  /** \brief . */
  static const std::string& getFileMustExistAttributeName() {
    static const std::string fileMustExistAttributeName = "fileMustExist";
    return fileMustExistAttributeName;
  }

  //@}

};


/**
 * \brief Converts StringValidators to and from XML.
 *
 * The valid XML represenation of a StringValidator is:
 * \code
  <Validator type="StringValidator"
   validatorId="Validator id"
  >
    <String value="Value 1"/>
    <String value="Value 2"/>
    ...Other String Values...
  </Validator>
 \endcode
 */
class TEUCHOSPARAMETERLIST_LIB_DLL_EXPORT StringValidatorXMLConverter : public ValidatorXMLConverter
{

public:

  /** \name Overridden from ValidatorXMLConverter */
  //@{

  /** \brief . */
  RCP<ParameterEntryValidator> convertXML(
    const XMLObject& xmlObj,
    const IDtoValidatorMap& validatorIDsMap) const;

  /** \brief . */
  void convertValidator(
    const RCP<const ParameterEntryValidator> validator,
    XMLObject& xmlObj,
    const ValidatortoIDMap& validatorIDsMap) const;

  #ifdef HAVE_TEUCHOS_DEBUG
  /** \brief . */
  RCP<const ParameterEntryValidator> getDummyValidator() const;
  #endif

  //@}

private:

  /** \name Private Members */
  //@{

  /** \brief . */
  static const std::string& getStringTagName() {
    static const std::string stringTagName = "String";
    return stringTagName;
  }

  /** \brief . */
  static const std::string& getStringValueAttributeName() {
    static const std::string stringValueAttributeName = "value";
    return stringValueAttributeName;
  }

  //@}

};

template<class ValidatorType, class EntryType>
class AbstractArrayValidatorXMLConverter : public ValidatorXMLConverter{
public:

  /** \name Overridden from ValidatorXMLConverter */
  //@{

  /** \brief . */
  RCP<ParameterEntryValidator> convertXML(
    const XMLObject& xmlObj,
    const IDtoValidatorMap& validatorIDsMap) const;

  /** \brief . */
  void convertValidator(
    const RCP<const ParameterEntryValidator> validator,
    XMLObject& xmlObj,
    const ValidatortoIDMap& validatorIDsMap) const;

  //@}

  /** \name Pure Virtual Fuctions */
  //@{

  /** \brief Returns a concrete validator that has
   * AbstractArrayValidator as it's parent class.
   */
  virtual RCP<AbstractArrayValidator<ValidatorType, EntryType> >
    getConcreteValidator(RCP<ValidatorType> prototypeValidator) const = 0;

  //@}
};


template<class ValidatorType, class EntryType>
RCP<ParameterEntryValidator>
AbstractArrayValidatorXMLConverter<ValidatorType, EntryType>::convertXML(
    const XMLObject& xmlObj,
    const IDtoValidatorMap& validatorIDsMap) const
{
  RCP<ValidatorType> prototypeValidator;
  if(xmlObj.hasAttribute(
    ValidatorXMLConverter::getPrototypeIdAttributeName()))
  {
    IDtoValidatorMap::const_iterator result =
      validatorIDsMap.find(
        xmlObj.getRequired<ParameterEntryValidator::ValidatorID>(
          getPrototypeIdAttributeName()));
    if (result != validatorIDsMap.end() ) {
      prototypeValidator =
        rcp_dynamic_cast<ValidatorType>(result->second, true);
    }
    else {
      TEUCHOS_TEST_FOR_EXCEPTION(true,
        MissingValidatorDefinitionException,
        "Could not find prototype validator with id: "
        << xmlObj.getRequired<ParameterEntryValidator::ValidatorID>(
          getPrototypeIdAttributeName()) << std::endl<< std::endl);
    }
  }
  else {
    prototypeValidator = rcp_dynamic_cast<ValidatorType>(
      ValidatorXMLConverterDB::convertXML(
        xmlObj.getChild(0), validatorIDsMap), true);
  }
  return getConcreteValidator(prototypeValidator);
}

template<class ValidatorType, class EntryType>
void
AbstractArrayValidatorXMLConverter<ValidatorType, EntryType>::convertValidator(
  const RCP<const ParameterEntryValidator> validator,
  XMLObject& xmlObj,
  const ValidatortoIDMap& validatorIDsMap) const
{
  RCP<const AbstractArrayValidator<ValidatorType, EntryType> > castedValidator =
    rcp_dynamic_cast<const AbstractArrayValidator<ValidatorType, EntryType> >(
      validator, true);
  if(validatorIDsMap.find(castedValidator->getPrototype())
    == validatorIDsMap.end())
  {
    xmlObj.addChild(ValidatorXMLConverterDB::convertValidator(
      castedValidator->getPrototype(), validatorIDsMap, false));
  }
  else{
    ParameterEntryValidator::ValidatorID prototypeID =
      validatorIDsMap.find(castedValidator->getPrototype())->second;

    xmlObj.addAttribute<ParameterEntryValidator::ValidatorID>(
      getPrototypeIdAttributeName(), prototypeID);
  }
}

/**
 * \brief Converts ArrayValidators to and from XML.
 *
 * ArrayValidators can be represented in XML one of two ways.
 * The first just creates the prototype validator as a child of
 * the ArrayValidator. In this case, the prototype validator does
 * NOT use a validatorId.
 *\code
  <Validator
   type="ArrayValidator(PrototypeValidatorType,ParameterArrayType)"
   validatorId="Validator id"
  >
     ...Prototype Validator Goes Here...
  </Validator>
 \endcode
 *
 * The second way to define an ArrayValidator in XML is to just use
 * the "prototypeId" attribute to specify the prototype validator as
 * some other validator you've already defined.
 * \code
   <Validator
     type="ArrayValidator(PrototypeValidatorType,ParameterArrayType)"
     validatorId="Validator id"
     prototypeId="Prototype Validator Id"
   />
 * \endcode
 */
template<class ValidatorType, class EntryType>
class ArrayValidatorXMLConverter :
  public AbstractArrayValidatorXMLConverter<ValidatorType, EntryType>
{
  /** @name Overridden from AbstractArrayValidatorXMLConverter */
  //@{

  virtual RCP<AbstractArrayValidator<ValidatorType, EntryType> > getConcreteValidator(
    RCP<ValidatorType> prototypeValidator) const
  {
    return rcp(new ArrayValidator<ValidatorType, EntryType>(prototypeValidator));
  }

#ifdef HAVE_TEUCHOS_DEBUG
  /** @name Overridden ValidatorXMLConverter*/
  //@{
  /** \brief . */
  RCP<const ParameterEntryValidator> getDummyValidator() const{
    return DummyObjectGetter<ArrayValidator<ValidatorType, EntryType> >::
      getDummyObject();
  }
  //@}
#endif
};

/**
 * \brief Converts TwoDArrayValidators to and from XML.
 *
 * TwoDArrayValidators can be represented in XML one of two ways.
 * The first just creates the prototype validator as a child of
 * the ArrayValidator. In this case, the prototype validator does
 * NOT use a validatorId.
 *\code
  <Validator
   type="TwoDArrayValidator(PrototypeValidatorType,ParameterArrayType)"
   validatorId="Validator id"
  >
     ...Prototype Validator Goes Here...
  </Validator>
 \endcode
 *
 * The second way to define an TwoDArrayValidator in XML is to just use
 * the "prototypeId" attribute to specify the prototype validator as
 * some other validator you've already defined.
 * \code
   <Validator
     type="TwoDArrayValidator(PrototypeValidatorType,ParameterArrayType)"
     validatorId="Validator id"
     prototypeId="Prototype Validator Id"
   />
 * \endcode
 */
template<class ValidatorType, class EntryType>
class TwoDArrayValidatorXMLConverter :
  public AbstractArrayValidatorXMLConverter<ValidatorType, EntryType>
{
  /** @name Overridden from AbstractArrayValidatorXMLConverter */
  //@{

  virtual RCP<AbstractArrayValidator<ValidatorType, EntryType> > getConcreteValidator(
    RCP<ValidatorType> prototypeValidator) const
  {
    return rcp(new TwoDArrayValidator<ValidatorType, EntryType>(prototypeValidator));
  }

  //@}

#ifdef HAVE_TEUCHOS_DEBUG
  /** @name Overridden ValidatorXMLConverter*/
  //@{
  /** \brief . */
  RCP<const ParameterEntryValidator> getDummyValidator() const{
    return DummyObjectGetter<TwoDArrayValidator<ValidatorType, EntryType> >::
      getDummyObject();
  }
  //@}
#endif

};



} // namespace Teuchos


#endif  // TEUCHOS_STANDARDVALIDATORXMLCONVERTERS_HPP

