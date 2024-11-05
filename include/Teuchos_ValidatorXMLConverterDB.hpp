// @HEADER
// *****************************************************************************
//                    Teuchos: Common Tools Package
//
// Copyright 2004 NTESS and the Teuchos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef TEUCHOS_VALIDATORXMLCONVERTERDB_HPP
#define TEUCHOS_VALIDATORXMLCONVERTERDB_HPP

/*! \file Teuchos_ValidatorXMLConverterDB.hpp
 * \brief A database for ValidatorXMLConverters.
*/

#include "Teuchos_ValidatorXMLConverter.hpp"


namespace Teuchos {

class ParameterEntryValidator;

/** \brief Provides ability to lookup ValidatorXMLConverterDB
 */
class TEUCHOSPARAMETERLIST_LIB_DLL_EXPORT ValidatorXMLConverterDB {
public:

  /** \name Modifier Functions */
  //@{

  /** \brief Add a converter to the database.
   *
   * \param validator A dummy validator representing the type of validator the
   * converter is designed to convert.
   * \param convertToAdd The converter to add to the database.
   */
  static void addConverter(RCP<const ParameterEntryValidator> validator,
    RCP<ValidatorXMLConverter> converterToAdd);

  //@}

  /** \name Converter Functions */
  //@{

  /** \brief Get an appropriate ValidatorXMLConverter given a
   * Validator.
   *
   * \param validator The ParameterEntryValidator for which a converter is
   * desired.
   */
  static RCP<const ValidatorXMLConverter> getConverter(
    const ParameterEntryValidator& validator);

  /** \brief Get an appropriate ValidatorXMLConverter given a XMLObject.
   *
   * @param xmlObject The XMLObject for which a converter is desired.
   */
  static RCP<const ValidatorXMLConverter>
    getConverter(const XMLObject& xmlObject);

  /**
   * \brief Given a validator converts the
   * validator to XML.
   *
   * \param validator The validator to be converter.
   * \param validatorIDsMap A map containing ParameterEntryValidators and their
   * associated IDs.
   * \param assignedID Whether or not the validator to be converted has been
   * assigned an ID and is therefore in the validatorIDsMap and should have a
   * ID attribute.
   *
   * \return XML representation of the validator.
   */
  static XMLObject convertValidator(
    RCP<const ParameterEntryValidator> validator,
    const ValidatortoIDMap& validatorIDsMap,
    bool assignedID=true);

  /**
   * \brief Given an XMLObject converts the XMLObject
   * to a ParameterEntryValidator and inserts the validator into the map.
   *
   * \param xmlObject The XMLObject representing the validator to be converted.
   * \param validatorIDsMap A map containing ParameterEntryValidators and their
   * associated IDs.
   * \return A ParameterEntryValidator that was represented by the XML.
   */
  static RCP<ParameterEntryValidator>
    convertXML(
      const XMLObject& xmlObject,
      const IDtoValidatorMap& validatorIDsMap);

  //@}

  /** \name I/O Functions */
  //@{

  /**
   * \brief prints the xml tags associated with all known converters
   *
   * \param out Stream to which tags should be printed.
   */
  static void printKnownConverters(std::ostream& out);

  //@}

private:

  /** \name Private types. */
  //@{

  /** \brief convience class. */
  typedef std::map<std::string, RCP<ValidatorXMLConverter> > ConverterMap;

  /** \brief convience typedef. */
  typedef std::pair<std::string, RCP<ValidatorXMLConverter> > ConverterPair;

  //@}

  /** \name Setup functions */
  //@{

  /** \brief Gets the default converter to be used to convert
   * Validators.
   *
   * This map is used to enable outside code to set up new converter types.
   */
  static ConverterMap& getConverterMap();

  //@}


};


} // end namespace Teuchos


//
// Helper Macros
//

// Doing this include so that when people macros like
// TEUCHOS_ADD_STRINGTOINTEGRALCONVERTER below they don't have to bother
// including it themselves. More likely they might not even know they have
// to.
#include "Teuchos_StandardValidatorXMLConverters.hpp"

/** \brief Add StringToIntegralParameterEntryValidator<INTEGRAL_TYPE> to set of
 * supported parameter types.
 */
#define TEUCHOS_ADD_STRINGTOINTEGRALVALIDATOR_CONVERTER(INTEGRALTYPE) \
  \
  Teuchos::ValidatorXMLConverterDB::addConverter( \
    Teuchos::DummyObjectGetter< \
      Teuchos::StringToIntegralParameterEntryValidator< INTEGRALTYPE > >:: \
        getDummyObject(), \
    Teuchos::rcp(new Teuchos::StringToIntegralValidatorXMLConverter< INTEGRALTYPE >));



/** \brief Add EnhancedNumberValidator<T> to the set of supported parameter
 * types.
 */
#define TEUCHOS_ADD_ENHANCEDNUMBERVALIDATOR_CONVERTER(T) \
  \
  Teuchos::ValidatorXMLConverterDB::addConverter( \
    Teuchos::DummyObjectGetter< \
      Teuchos::EnhancedNumberValidator< T > >:: \
        getDummyObject(), \
    Teuchos::rcp(new Teuchos::EnhancedNumberValidatorXMLConverter< T >));


/** \brief Add ArrayValidator<VALIDATORTYPE, ENTRYTYPE> to set of supported
 * parameter types.
 */
#define TEUCHOS_ADD_ARRAYVALIDATOR_CONVERTER(VALIDATORTYPE, ENTRYTYPE) \
  \
  Teuchos::ValidatorXMLConverterDB::addConverter( \
    Teuchos::DummyObjectGetter< \
      Teuchos::ArrayValidator< VALIDATORTYPE, ENTRYTYPE > >:: \
        getDummyObject(), \
    Teuchos::rcp(new Teuchos::ArrayValidatorXMLConverter< VALIDATORTYPE, ENTRYTYPE >)); \
  Teuchos::ValidatorXMLConverterDB::addConverter( \
    Teuchos::DummyObjectGetter< \
      Teuchos::TwoDArrayValidator< VALIDATORTYPE, ENTRYTYPE > >:: \
        getDummyObject(), \
    Teuchos::rcp(new Teuchos::TwoDArrayValidatorXMLConverter< VALIDATORTYPE, ENTRYTYPE >));


/** \brief Add numeric parameter types for type T. */
#define TEUCHOS_ADD_NUMBERTYPE_VALIDATOR_CONVERTERS(T) \
  TEUCHOS_ADD_STRINGTOINTEGRALVALIDATOR_CONVERTER(T); \
  TEUCHOS_ADD_ENHANCEDNUMBERVALIDATOR_CONVERTER(T); \
  TEUCHOS_ADD_ARRAYVALIDATOR_CONVERTER(Teuchos::EnhancedNumberValidator< T >, T );

/** \brief Add a validator converter of type CONVERTER_TYPE which converts
 *  validators of VALIDATOR_TYPE to the map CONVERTER_MAP.
 */
#define TEUCHOS_ADD_VALIDATOR_CONVERTER(VALIDATOR_TYPE, CONVERTER_TYPE) \
  Teuchos::ValidatorXMLConverterDB::addConverter( \
      Teuchos::DummyObjectGetter< VALIDATOR_TYPE > \
      ::getDummyObject(), \
      Teuchos::rcp(new CONVERTER_TYPE ));

#endif // TEUCHOS_VALIDATORXMLCONVERTERDB_HPP
