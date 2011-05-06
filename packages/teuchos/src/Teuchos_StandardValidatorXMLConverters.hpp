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


/** \brief Converts StringToIntegralParameterEntryValidators to and from XML.
 *
 * The XML Representation for a StringToIntegralValidator is:
 * \code
  <Validator 
    type="StringToIntegralValidator(NumberType)"
    defaultParameterName="Name of default parameter"
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
 *
 * The "integralValue" and "stringDoc" XML attributes are optional. However,
 * if one of the "String" tags includes an "integralValue" and/or a "stringDoc"
 * XML attribute, all other "String" tags must do so as well.
 */
template<class IntegralType>
class TEUCHOS_LIB_DLL_EXPORT StringToIntegralValidatorXMLConverter : 
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

  /** \brief . */
  static const std::string& getDefaultParameterAttributeName() {
    static const std::string defaultParameterAttributeName_ =
      "defaultParameterName";
    return defaultParameterAttributeName_;
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
    TEST_FOR_EXCEPTION(currentChild.getTag() != getStringTagName(), 
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

  if(stringDocs.size() != 0 && integralValues.size() != 0){
    return stringToIntegralParameterEntryValidator<IntegralType>(
      strings, stringDocs, integralValues(), defaultParameterName);
  }
  else if(integralValues.size() != 0){
    return stringToIntegralParameterEntryValidator<IntegralType>(
      strings, integralValues(), defaultParameterName);
  }
  else{
    return stringToIntegralParameterEntryValidator<IntegralType>(
      strings, defaultParameterName);
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
  xmlObj.addAttribute(getIntegralValueAttributeName(),
    TypeNameTraits<IntegralType>::name());
}

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
class TEUCHOS_LIB_DLL_EXPORT AnyNumberValidatorXMLConverter : public ValidatorXMLConverter
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
class TEUCHOS_LIB_DLL_EXPORT EnhancedNumberValidatorXMLConverter : public ValidatorXMLConverter
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
class TEUCHOS_LIB_DLL_EXPORT FileNameValidatorXMLConverter : public ValidatorXMLConverter
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
class TEUCHOS_LIB_DLL_EXPORT StringValidatorXMLConverter : public ValidatorXMLConverter
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
class TEUCHOS_LIB_DLL_EXPORT ArrayValidatorXMLConverter : public ValidatorXMLConverter
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
    return DummyObjectGetter<ArrayValidator<ValidatorType, EntryType> >::
      getDummyObject();
  }
#endif
  
  //@}

};


template<class ValidatorType, class EntryType>
RCP<ParameterEntryValidator>
ArrayValidatorXMLConverter<ValidatorType, EntryType>::convertXML(
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
      TEST_FOR_EXCEPTION(true,
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
  return rcp(new ArrayValidator<ValidatorType, EntryType>(prototypeValidator));
}


template<class ValidatorType, class EntryType>
void 
ArrayValidatorXMLConverter<ValidatorType, EntryType>::convertValidator(
  const RCP<const ParameterEntryValidator> validator,
  XMLObject& xmlObj,
  const ValidatortoIDMap& validatorIDsMap) const
{
  RCP<const ArrayValidator<ValidatorType, EntryType> > castedValidator = 
    rcp_dynamic_cast<const ArrayValidator<ValidatorType, EntryType> >(
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

} // namespace Teuchos


#endif  // TEUCHOS_STANDARDVALIDATORXMLCONVERTERS_HPP

