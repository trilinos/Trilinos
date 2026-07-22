// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ROL_PARAMETERENTRYXMLCONVERTERDB_HPP
#define ROL_PARAMETERENTRYXMLCONVERTERDB_HPP

#include "ROL_StandardParameterEntryXMLConverters.hpp"


/*! \file ROL_ParameterEntryXMLCoverterDB.hpp
 * \brief A database for ParameterEntryXMLConverters.
*/


namespace ROL {

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


} // namespace ROL

//
// Helper Macros
//


/**
 * Add a converter of type T to map CONVERTER_MAP
 */
#define TEUCHOS_ADD_TYPE_CONVERTER(T) \
  \
  ROL::ParameterEntryXMLConverterDB::addConverter( \
    Teuchos::rcp(new ROL::StandardTemplatedParameterConverter< T >));

#endif // ROL_PARAMETERENTRYXMLCONVERTERDB_HPP
