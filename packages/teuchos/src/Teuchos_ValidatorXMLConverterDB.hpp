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
class ValidatorXMLConverterDB {
public:

  /** \name Public types. */
  //@{
  
  /** \brief convience class. */
  typedef std::map<std::string, RCP<ValidatorXMLConverter> > ConverterMap;

  /** \brief convience typedef. */
  typedef std::pair<std::string, RCP<ValidatorXMLConverter> > ConverterPair;

  //@}

  /** \name Modifier Functions */
  //@{
  
  /** \brief Add a converter to the database.
   *
   * \param validator A dummy validator representing the type of validator the
   * converter is designed to convert.
   * \param convertToAdd The converter to add to the database.
   */
  static void addConverter(RCP<ParameterEntryValidator> validator,
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


/** \brief Add StringToIntegralParameterEntryValidator<INTEGRAL_TYPE> to set of
 * supported parameter types.
 */
#define TEUCHOS_ADD_STRINGTOINTEGRALCONVERTER(CONVERTER_MAP, INTEGRALTYPE) \
  \
  (CONVERTER_MAP).insert(Teuchos::ValidatorXMLConverterDB::ConverterPair( \
    Teuchos::DummyObjectGetter< \
      Teuchos::StringToIntegralParameterEntryValidator< INTEGRALTYPE > >:: \
        getDummyObject()->getXMLTypeName(), \
    Teuchos::rcp(new Teuchos::StringToIntegralValidatorXMLConverter< INTEGRALTYPE >)));


/** \brief Add EnhancedNumberValidator<T> to the set of supported parameter 
 * types.
 */
#define TEUCHOS_ADD_ENHANCEDNUMBERCONVERTER(CONVERTER_MAP, T) \
  \
  (CONVERTER_MAP).insert(Teuchos::ValidatorXMLConverterDB::ConverterPair( \
    Teuchos::DummyObjectGetter<Teuchos::EnhancedNumberValidator< T > >:: \
      getDummyObject()->getXMLTypeName(), \
    Teuchos::rcp(new Teuchos::EnhancedNumberValidatorXMLConverter< T >)));


/** \brief Add ArrayValidator<VALIDATORTYPE, ENTRYTYPE> to set of supported
 * parameter types.
 */
#define TEUCHOS_ADD_ARRAYCONVERTER(CONVERTER_MAP, VALIDATORTYPE, ENTRYTYPE) \
  \
  (CONVERTER_MAP).insert(Teuchos::ValidatorXMLConverterDB::ConverterPair( \
    Teuchos::DummyObjectGetter<Teuchos::ArrayValidator< VALIDATORTYPE , ENTRYTYPE > >:: \
      getDummyObject()->getXMLTypeName(), \
    Teuchos::rcp(new Teuchos::ArrayValidatorXMLConverter< VALIDATORTYPE, ENTRYTYPE >)));


/** \brief Add numeric parameter types for type T. */
#define TEUCHOS_ADD_NUMBERTYPECONVERTERS(CONVERTER_MAP, T) \
  TEUCHOS_ADD_STRINGTOINTEGRALCONVERTER(CONVERTER_MAP, T ); \
  TEUCHOS_ADD_ENHANCEDNUMBERCONVERTER(CONVERTER_MAP, T ); \
  TEUCHOS_ADD_ARRAYCONVERTER(CONVERTER_MAP, Teuchos::EnhancedNumberValidator< T >, T );


#endif // TEUCHOS_VALIDATORXMLCONVERTERDB_HPP
