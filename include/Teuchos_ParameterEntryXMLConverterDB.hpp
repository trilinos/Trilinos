// @HEADER
// *****************************************************************************
//                    Teuchos: Common Tools Package
//
// Copyright 2004 NTESS and the Teuchos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef TEUCHOS_PARAMETERENTRYXMLCONVERTERDB_HPP
#define TEUCHOS_PARAMETERENTRYXMLCONVERTERDB_HPP

#include "Teuchos_StandardParameterEntryXMLConverters.hpp"
#include "Teuchos_XMLParameterListExceptions.hpp"
#include "Teuchos_XMLParameterListWriter.hpp"


/*! \file Teuchos_ParameterEntryXMLCoverterDB.hpp
 * \brief A database for ParameterEntryXMLConverters.
*/


namespace Teuchos {

/** \brief Provides ability to lookup ParameterEntryXMLConverters
 */
class TEUCHOSPARAMETERLIST_LIB_DLL_EXPORT ParameterEntryXMLConverterDB {
public:

  /** \name Modifier Functions */
  //@{

  /** \brief Add a converter to the database.
   *
   * \param convertToAdd The converter to add to the database.
   */
  static void addConverter(RCP<ParameterEntryXMLConverter> converterToAdd){
    getConverterMap().insert(
      ConverterPair(converterToAdd->getTypeAttributeValue(), converterToAdd));
  }

  //@}

  /** \name Getter Functions */
  //@{


  /** \brief Get an appropriate ParameterEntryXMLConverter given a ParameterEntry.
   *
   * \param entry The ParameterEntry for which a converter is desired.
   */
  static RCP<const ParameterEntryXMLConverter>
    getConverter(RCP<const ParameterEntry> entry);

  /** \brief Get an appropriate ParameterEntryXMLConverter given a XMLObject.
   *
   * \param xmlObject The XMLObject for which a converter is desired.
   */
  static RCP<const ParameterEntryXMLConverter>
    getConverter(const XMLObject& xmlObject);

  /** \brief Gets the default converter to be used on Parameter Entries */
  static RCP<const ParameterEntryXMLConverter> getDefaultConverter();

  //@}

  // 2010/07/30: rabarlt: The above two functions should be moved into
  // Teuchos_ParameterEntryXMLConvergerDB.cpp.  These functions don't need to
  // be inlined and it will be easier to set breakpoints in the debugger if
  // they are in a *.cpp file.

  /** \name Converter Functions */
  //@{

  /**
   * \brief Converts the given ParameterEntry to XML.
   */
  static XMLObject convertEntry(
    RCP<const ParameterEntry> entry,
    const std::string& name,
    const ParameterEntry::ParameterEntryID& id,
    const ValidatortoIDMap& validatorIDsMap)
  {
    return getConverter(entry)->fromParameterEntrytoXML(
      entry, name, id, validatorIDsMap);
  }

  /**
   * \brief Converts XML to a ParameterEntry.
   */
  static ParameterEntry convertXML(const XMLObject& xmlObj)
  {
    return getConverter(xmlObj)->fromXMLtoParameterEntry(xmlObj);
  }

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

  /** \brief convience typedef */
  typedef std::map<std::string, RCP<ParameterEntryXMLConverter> > ConverterMap;

  /** \brief convience typedef */
  typedef std::pair<std::string, RCP<ParameterEntryXMLConverter> > ConverterPair;

  //@}

  /** \brief Gets the map containing all the ParameterEntry converters. */
  static ConverterMap& getConverterMap();


};


} // namespace Teuchos

//
// Helper Macros
//


/**
 * Add a converter of type T to map CONVERTER_MAP
 */
#define TEUCHOS_ADD_TYPE_CONVERTER(T) \
  \
  Teuchos::ParameterEntryXMLConverterDB::addConverter( \
    Teuchos::rcp(new Teuchos::StandardTemplatedParameterConverter< T >));

/**
 * Add a converter of type Array<T> to map CONVERTER_MAP
 */
#define TEUCHOS_ADD_ARRAYTYPE_CONVERTER(T) \
  Teuchos::ParameterEntryXMLConverterDB::addConverter( \
    Teuchos::rcp(new Teuchos::StandardTemplatedParameterConverter< Teuchos::Array< T > >)); \
  Teuchos::ParameterEntryXMLConverterDB::addConverter( \
    Teuchos::rcp(new Teuchos::StandardTemplatedParameterConverter< Teuchos::TwoDArray< T > >));

/**
 * Add both a converter for type T and Array<T> to CONVERTER_MAP
 */
#define TEUCHOS_ADD_TYPE_AND_ARRAYTYPE_CONVERTER(T) \
  \
  TEUCHOS_ADD_TYPE_CONVERTER(T); \
  TEUCHOS_ADD_ARRAYTYPE_CONVERTER(T);


#endif // TEUCHOS_PARAMETERENTRYXMLCONVERTERDB_HPP
