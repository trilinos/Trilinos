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
#define TEUCHOS_STANDARDVALIDATORXMLCONVERTERSL_HPP

#include "Teuchos_ValidatorXMLConverter.hpp"
#include "Teuchos_StandardParameterEntryValidators.hpp"
#include "Teuchos_ValidatorXMLConverterDB.hpp"
#include "Teuchos_XMLParameterListReader.hpp"
#include "Teuchos_DummyObjectGetter.hpp"


namespace Teuchos {


/** \brief Converts StringToIntegralParameterEntryValidators to and from XML.
 */
template<class IntegralType>
class StringToIntegralValidatorXMLConverter : 
  public ValidatorXMLConverter
{

public:

  /** \brief . */
  RCP<ParameterEntryValidator>
  convertXML(const XMLObject& xmlObj, IDtoValidatorMap& validatorMap) const;

  /** \brief . */
  XMLObject convertValidator(const RCP<const ParameterEntryValidator> validator, 
    const ValidatortoIDMap& validatorMap) const;

  #ifdef HAVE_TEUCHOS_DEBUG
  /** \brief . */
  RCP<const ParameterEntryValidator > 
  getDummyValidator() const{
    return DummyObjectGetter<
	  StringToIntegralParameterEntryValidator<IntegralType> >::getDummyObject();
  }
  #endif

private:

  static const std::string& getIntegralValueAttributeName() {
    static const std::string integralValueAttributeName_ = "integralvalue";
    return integralValueAttributeName_;
  }

  static const std::string& getStringTagName() {
    static const std::string stringTagName_ = "string";
    return stringTagName_;
  }

  static const std::string& getStringValueAttributeName() {
    static const std::string stringValueAttributeName_ = "stringvalue";
    return stringValueAttributeName_;
  }

  static const std::string& getStringDocAttributeName() {
    static const std::string stringDocAttributeName_ = "stringdoc";
    return stringDocAttributeName_;
  }

  static const std::string& getDefaultParameterAttributeName() {
    static const std::string defaultParameterAttributeName_ = "defaultparametername";
    return defaultParameterAttributeName_;
  }

};


//
// Implementations
//


template<class IntegralType>
RCP<ParameterEntryValidator>
StringToIntegralValidatorXMLConverter<IntegralType>::convertXML(
  const XMLObject& xmlObj, IDtoValidatorMap& validatorMap) const
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
        currentChild.getRequired<IntegralType>(getIntegralValueAttributeName()));
    }
    if (currentChild.hasAttribute(getStringDocAttributeName())) {
      stringDocs.append(currentChild.getRequired<std::string>(getStringDocAttributeName()));
    }
  }
  std::string defaultParameterName = xmlObj.getRequired(getDefaultParameterAttributeName());
  if (stringDocs.size() != 0 && integralValues.size() != 0) {
    return stringToIntegralParameterEntryValidator<IntegralType>(strings, stringDocs,
      integralValues, defaultParameterName);
  }
  else if (stringDocs.size() == 0 && integralValues.size() !=0) {
    return stringToIntegralParameterEntryValidator<IntegralType>(strings,
      integralValues, defaultParameterName);
  }
  else {
    return stringToIntegralParameterEntryValidator<IntegralType>(strings,
      defaultParameterName);
  }
}


template<class IntegralType>
XMLObject StringToIntegralValidatorXMLConverter<IntegralType>::convertValidator(
  const RCP<const ParameterEntryValidator> validator, 
  const ValidatortoIDMap& validatorMap) const
{
  RCP<const StringToIntegralParameterEntryValidator<IntegralType> > 
    castedValidator =
    rcp_dynamic_cast<const StringToIntegralParameterEntryValidator<IntegralType> >(
      validator, true);

  XMLObject toReturn(validator->getXMLTagName());

  RCP<const Array<std::string> > stringValues =
    castedValidator->validStringValues();
  RCP<const Array<std::string> > stringDocValues = 
    castedValidator->getStringDocs();

  bool hasStringDocs = !(stringDocValues.is_null()) && (stringDocValues->size() != 0);
  for (int i =0; i<stringValues->size(); ++i) {
    XMLObject stringTag(getStringTagName());
    stringTag.addAttribute(getStringValueAttributeName(), (*stringValues)[i]);
    stringTag.addAttribute(getIntegralValueAttributeName(),
      castedValidator->getIntegralValue((*stringValues)[i]));
    if (hasStringDocs) {
      stringTag.addAttribute(getStringDocAttributeName(), (*stringDocValues)[i]);
    }
    toReturn.addChild(stringTag);
  }
  toReturn.addAttribute(getDefaultParameterAttributeName(),
    castedValidator->getDefaultParameterName());
  toReturn.addAttribute(getIntegralValueAttributeName(),
    TypeNameTraits<IntegralType>::name());
  return toReturn;
}

/*
 * \brief Converts AnyNumberParameterEntryValidators to and from XML.
 */
class AnyNumberValidatorXMLConverter : public ValidatorXMLConverter
{

public:

  /** \brief . */
  RCP<ParameterEntryValidator> convertXML(const XMLObject& xmlObj, 
    IDtoValidatorMap& validatorMap) const;

  /** \brief . */
  XMLObject convertValidator(const RCP<const ParameterEntryValidator> validator, 
    const ValidatortoIDMap& validatorMap) const;

  #ifdef HAVE_TEUCHOS_DEBUG
  /** \brief . */
  RCP<const ParameterEntryValidator> getDummyValidator() const;
  #endif

private:

  static const std::string& getAllowIntAttributeName() {
    static const std::string allowIntAttributeName_ = "allowInt";
    return allowIntAttributeName_;
  }

  static const std::string& getAllowDoubleAttributeName() {
    static const std::string allowDoubleAttributeName_ = "allowDouble";
    return allowDoubleAttributeName_;
  }

  static const std::string& getAllowStringAttributeName() {
    static const std::string allowStringAttributeName_ = "allowString";
    return allowStringAttributeName_;
  }
  
  static const std::string& getPrefferedTypeAttributeName() {
    static const std::string prefferedTypeAttributeName_ = "prefferedType";
    return prefferedTypeAttributeName_;
  }
 
};


/** \brief Converts EnhancedNumberValidators to and from XML.
 */
template<class T>
class EnhancedNumberValidatorXMLConverter : public ValidatorXMLConverter
{

public:

  /** \brief . */
  RCP<ParameterEntryValidator> convertXML(const XMLObject& xmlObj,
    IDtoValidatorMap& validatorMap) const;

  /** \brief . */
  XMLObject convertValidator(const RCP<const ParameterEntryValidator> validator,
    const ValidatortoIDMap& validatorMap) const;

#ifdef HAVE_TEUCHOS_DEBUG
  /** \brief . */
  RCP<const ParameterEntryValidator> getDummyValidator() const{
    return DummyObjectGetter<EnhancedNumberValidator<T> >::getDummyObject();
  }
#endif

private:

  static const std::string& getMinAttributeName() {
    static const std::string minAttributeName = "min";
    return minAttributeName;
  }

  static const std::string& getMaxAttributeName() {
    static const std::string maxAttributeName = "max";
    return maxAttributeName;
  }

  static const std::string& getStepAttributeName() {
    static const std::string stepAttributeName = "step";
    return stepAttributeName;
  }

  static const std::string& getPrecisionAttributeName() {
    static const std::string precisionAttributeName = "precision";
    return precisionAttributeName;
  }

};


template<class T>
RCP<ParameterEntryValidator> EnhancedNumberValidatorXMLConverter<T>::convertXML(
  const XMLObject& xmlObj, IDtoValidatorMap& validatorMap) const
{
  RCP<EnhancedNumberValidator<T> > toReturn = rcp(new EnhancedNumberValidator<T>());
  xmlObj.getWithDefault(getStepAttributeName(), EnhancedNumberTraits<T>::defaultStep()),
  xmlObj.getWithDefault(getPrecisionAttributeName(),
    EnhancedNumberTraits<T>::defaultPrecision());
  if (xmlObj.hasAttribute(getMinAttributeName())) {
    toReturn->setMin(xmlObj.getRequired<T>(getMinAttributeName()));
  }
  if (xmlObj.hasAttribute(getMaxAttributeName())) {
    toReturn->setMax(xmlObj.getRequired<T>(getMaxAttributeName()));
  }
  return toReturn;
}


template<class T>
XMLObject EnhancedNumberValidatorXMLConverter<T>::convertValidator(
  const RCP<const ParameterEntryValidator > validator,
  const ValidatortoIDMap& validatorMap) const
{
  RCP<const EnhancedNumberValidator<T> > castedValidator =
    rcp_dynamic_cast<const EnhancedNumberValidator<T> >(validator, true);
  XMLObject toReturn(castedValidator->getXMLTagName());
  if (castedValidator->hasMin()) {
    toReturn.addAttribute<T>(getMinAttributeName(), castedValidator->getMin());
  }
  if (castedValidator->hasMax()) {
    toReturn.addAttribute<T>(getMaxAttributeName(), castedValidator->getMax());
  }
  toReturn.addAttribute<T>(getStepAttributeName(), castedValidator->getStep());
  toReturn.addAttribute<T>(getPrecisionAttributeName(), castedValidator->getPrecision());
  return toReturn;
}


/*
 * \brief Converts FileNameValidators to and from XML.
 */
class FileNameValidatorXMLConverter : public ValidatorXMLConverter
{

public:

  /** \brief . */
  RCP<ParameterEntryValidator> convertXML(const XMLObject& xmlObj, 
    IDtoValidatorMap& validatorMap) const;

  /** \brief . */
  XMLObject convertValidator(const RCP<const ParameterEntryValidator> validator,
    const ValidatortoIDMap& validatorMap) const;

  #ifdef HAVE_TEUCHOS_DEBUG
  /** \brief . */
  RCP<const ParameterEntryValidator> getDummyValidator() const;
  #endif

private:

  static const std::string& getFileMustExistAttributeName() {
    static const std::string fileMustExistAttributeName = "filemustexist";
    return fileMustExistAttributeName;
  }
  
};


/*
 * \brief Converts StringValidators to and from XML.
 */
class StringValidatorXMLConverter : public ValidatorXMLConverter
{

public:

  /** \brief . */
  RCP<ParameterEntryValidator> convertXML(const XMLObject& xmlObj,
    IDtoValidatorMap& validatorMap) const;

  /** \brief . */
  XMLObject convertValidator(const RCP<const ParameterEntryValidator> validator,
    const ValidatortoIDMap& validatorMap) const;

  #ifdef HAVE_TEUCHOS_DEBUG
  /** \brief . */
  RCP<const ParameterEntryValidator> getDummyValidator() const;
  #endif

private:
  
  static const std::string& getStringTagName() {
    static const std::string stringTagName = "string";
    return stringTagName;
  }

  static const std::string& getStringValueAttributeName() {
    static const std::string stringValueAttributeName = "stringvalue";
    return stringValueAttributeName;
  }

};


/*
 * \brief Converts ArrayValidators to and from XML.
 */
template<class ValidatorType, class EntryType>
class ArrayValidatorXMLConverter : public ValidatorXMLConverter
{
public:

  /** \brief . */
  RCP<ParameterEntryValidator> convertXML(
    const XMLObject& xmlObj, 
	IDtoValidatorMap& validatorMap) const;

  /** \brief . */
  XMLObject convertValidator(
    const RCP<const ParameterEntryValidator> validator, 
	const ValidatortoIDMap& validatorMap) const;

#ifdef HAVE_TEUCHOS_DEBUG
  /** \brief . */
  RCP<const ParameterEntryValidator> getDummyValidator() const{
    return DummyObjectGetter<ArrayValidator<ValidatorType, EntryType> >::getDummyObject();
  }
#endif

};


template<class ValidatorType, class EntryType>
RCP<ParameterEntryValidator>
ArrayValidatorXMLConverter<ValidatorType, EntryType>::convertXML(
  const XMLObject& xmlObj, IDtoValidatorMap& validatorMap) const
{
  RCP<ValidatorType> prototypeValidator;
  if (xmlObj.hasAttribute(getPrototypeIdAttributeName())) {
    IDtoValidatorMap::const_iterator result =
      validatorMap.getValidator(xmlObj.getRequiredInt(getPrototypeIdAttributeName()));
    if (result != validatorMap.end()) {
      prototypeValidator = rcp_dynamic_cast<ValidatorType>(result->second, true);
    }
    else {
      TEST_FOR_EXCEPTION(true,
        std::runtime_error,
        "Could not find prototype validtor with id: "
        << xmlObj.getRequiredInt(getPrototypeIdAttributeName()) << "\n");
    }
  }
  else {
    prototypeValidator = rcp_dynamic_cast<ValidatorType>(
	  ValidatorXMLConverterDB::convertXML(xmlObj.getChild(0), validatorMap), true);
  }
  return rcp(new ArrayValidator<ValidatorType, EntryType>(prototypeValidator));
}


template<class ValidatorType, class EntryType>
XMLObject ArrayValidatorXMLConverter<ValidatorType, EntryType>::convertValidator(
  const RCP<const ParameterEntryValidator> validator, 
  const ValidatortoIDMap& validatorMap) const
{
  XMLObject toReturn(validator->getXMLTagName());
  RCP<const ArrayValidator<ValidatorType, EntryType> > castedValidator = 
    rcp_dynamic_cast<const ArrayValidator<ValidatorType, EntryType> >(validator, true);
  if (validatorMap.getID(castedValidator->getPrototype()) != validatorMap.end()) {
    toReturn.addInt(getPrototypeIdAttributeName(),
      validatorMap.getID(castedValidator->getPrototype())->second);
  }
  else {
    toReturn.addChild(
	  ValidatorXMLConverterDB::convertValidator(
	    castedValidator->getPrototype(), validatorMap));
  }
  return toReturn;
}

} // namespace Teuchos


#endif  // TEUCHOS_STANDARDVALIDATORXMLCONVERTERS_HPP

