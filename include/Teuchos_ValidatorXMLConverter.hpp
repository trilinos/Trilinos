// @HEADER
// *****************************************************************************
//                    Teuchos: Common Tools Package
//
// Copyright 2004 NTESS and the Teuchos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef TEUCHOS_VALIDATORXMLCONVERTER_HPP
#define TEUCHOS_VALIDATORXMLCONVERTER_HPP

/*! \file Teuchos_ValidatorXMLConverter.hpp
 * \brief Converts back and forth between XML
 * and ParameterEntryValidators.
*/

#include "Teuchos_XMLObject.hpp"
#include "Teuchos_Describable.hpp"
#include "Teuchos_XMLParameterListExceptions.hpp"
#include "Teuchos_ParameterEntryValidator.hpp"
#include "Teuchos_XMLParameterListReader.hpp"
#include "Teuchos_XMLParameterListWriter.hpp"


namespace Teuchos {


/** \brief An abstract base class for converting ParameterEntryValidators to
 * and from XML.
 */
class TEUCHOSPARAMETERLIST_LIB_DLL_EXPORT ValidatorXMLConverter : public Describable {

public:

  /** \name Converter Functions */
  //@{

  /** \brief Converts a given XMLObject to a ParameterEntryValidator.
   *
   * @param xmlObj The XMLObject to convert to a ParameterEntryValidator.
   * @param validatorIDsMap A map containing ParameterEntryValidators and their
   * associated IDs.
   * @return The converted ParameterEntryValidator.
   */
  RCP<ParameterEntryValidator>
  fromXMLtoValidator(
    const XMLObject& xmlObj,
    const IDtoValidatorMap& validatorIDsMap) const;

  /** \brief Preforms any and all special xml conversion that
   * is specific to a
   * particular ParameterEntryValidator.
   *
   * @param xmlObj The xml to be converted.
   * @param validatorIDsMap A map containing ParameterEntryValidators and their
   * associated IDs.
   * @return The converted ParameterEntryValidator.
   */
  virtual RCP<ParameterEntryValidator>
    convertXML(const XMLObject& xmlObj,
    const IDtoValidatorMap& validatorIDsMap) const=0;

  /** \brief Converters a given ParameterEntryValidator to XML.
   *
   * @param validator The ParameterEntryValidator to be converted to XML.
   * @param validatorIDsMap A map containing ParameterEntryValidators and their
   * associated IDs.
   * @param assignedID Whether or not the validator to be converted has been
   * assigned an ID and is therefore in the validatorIDsMap and should have a
   * ID attribute.
   * @return An XML representation of the given ParameterEntryValidator.
   */
  XMLObject fromValidatortoXML(
    const RCP<const ParameterEntryValidator> validator,
    const ValidatortoIDMap& validatorIDsMap,
    bool assignedID=true) const;

  /** \brief Preforms any and all special validator conversion that is
   * specific to a particlar ParameterEntryValidator
   *
   * @param validator The validator to be converted.
   * @param xmlObj The XMLObject to store all serialization in.
   * @param validatorIDsMap A map containing ParameterEntryValidators and their
   * associated IDs.
   * being converted.
   */
  virtual void convertValidator(
    const RCP<const ParameterEntryValidator> validator,
    XMLObject& xmlObj,
    const ValidatortoIDMap& validatorIDsMap) const = 0;

  //@}

  #ifdef HAVE_TEUCHOS_DEBUG
  /** \name Debug Functions */
  //@{

  /**
   * \brief Returns dummy validator of the type the converter is designed to
   * convert.
   *
   * @return A default dummy validator.
   */
  virtual Teuchos::RCP<const ParameterEntryValidator>
    getDummyValidator() const = 0;

  //@}
  #endif

  //! \name Attribute/Query Functions
  //@{

  /** \brief . */
  static const std::string& getIdAttributeName(){
    static const std::string idAttributeName = "validatorId";
    return idAttributeName;
  }

  /** \brief . */
  static const std::string& getPrototypeIdAttributeName(){
    static const std::string prototypeIdAttributeName = "prototypeId";
    return prototypeIdAttributeName;
  }

  /** \brief . */
  static const std::string& getTypeAttributeName(){
    static const std::string typeAttributeName = "type";
    return typeAttributeName;
  }

  /** \brief . */
  static const std::string& getValidatorTagName(){
    static const std::string validatorTagName = "Validator";
    return validatorTagName;
  }
  //@}

};


} // namespace Teuchos


#endif // TEUCHOS_VALIDATORXMLCONVERTER_HPP
