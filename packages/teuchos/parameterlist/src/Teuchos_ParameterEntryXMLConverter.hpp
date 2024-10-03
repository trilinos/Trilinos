// @HEADER
// *****************************************************************************
//                    Teuchos: Common Tools Package
//
// Copyright 2004 NTESS and the Teuchos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef TEUCHOS_PARAMETERENTRYXMLCONVERTER_HPP
#define TEUCHOS_PARAMETERENTRYXMLCONVERTER_HPP


/*! \file Teuchos_ParameterEntryXMLCoverter.hpp
 *  \brief The base class for all ParameterEntryXMLConverters.
*/


#include "Teuchos_ParameterEntry.hpp"
#include "Teuchos_XMLObject.hpp"
#include "Teuchos_Describable.hpp"
#include "Teuchos_XMLParameterListWriter.hpp"


namespace Teuchos {

/** \brief A class used to convert parameter entries to xml and vice versa.
 */
class TEUCHOSPARAMETERLIST_LIB_DLL_EXPORT ParameterEntryXMLConverter : public Describable {

public:

  /** \name Converter Functions */
  //@{

  /** \brief Converts the given xml into a parameter entry.
   *
   * \param xmlObj The xml to be converted to a parameter entry.
   * \returns A ParameterEntry with the aspects specified by the xml.
   */
  ParameterEntry fromXMLtoParameterEntry(const XMLObject &xmlObj) const;

  /** \brief Converts the given parameter entry to xml.
   *
   * \param entry The parameter entry to convert to xml.
   * \param name The name associated with the parameter entry.
   * \returns An XMLObject representing the parameter entry.
   */
  XMLObject fromParameterEntrytoXML(
    RCP<const ParameterEntry> entry,
    const std::string &name,
    const ParameterEntry::ParameterEntryID& id,
    const ValidatortoIDMap& validatorIDsMap) const;

  virtual any getAny(const XMLObject& xmlObj) const=0;

  //@}

  //! \name Attribute/Query Methods
  //@{

  /** \brief Gets a string representing the value that should be assigned to
   * the "type" attribute when converting a parameter entry to xml.
   *
   * \returns The value to be assigned to the "type" attribute when converting
   * a parameter entry to xml.
   */
  virtual const std::string getTypeAttributeValue() const=0;

  /** \brief Gets the value to be assigned to the "value" attribute when
   * converting the paramter entry to xml.
   *
   * \param entry The entry being converted.
   *
   * \returns The value to be assigned to the "value" attribute when
   * converting the parameter entry to xml.
   */
  virtual const std::string getValueAttributeValue(
    RCP<const ParameterEntry > entry) const=0;

  /** \brief . */
  static const std::string& getTypeAttributeName() {
    static const std::string typeAttributeName_ = "type";
    return typeAttributeName_;
  }

  /** \brief . */
  static const std::string& getIdAttributeName() {
    static const std::string idAttributeName_ = "id";
    return idAttributeName_;
  }

  /** \brief . */
  static const std::string& getValueAttributeName() {
    static const std::string valueAttributeName_ = "value";
    return valueAttributeName_;
  }

  //@}

private:

  /** \name Private Members */
  //@{

  /** \brief . */
  static const std::string& getDefaultAttributeName() {
    static const std::string defaultAttributeName_ = "isDefault";
    return defaultAttributeName_;
  }

  /** \brief . */
  static const std::string& getUsedAttributeName() {
    static const std::string usedAttributeName_ = "isUsed";
    return usedAttributeName_;
  }

  /** \brief . */
  static const std::string& getDocStringAttributeName() {
    static const std::string docStringAttributeName_ = "docString";
    return docStringAttributeName_;
  }

  //@}

};


} // namespace Teuchos


#endif // TEUCHOS_PARAMETERENTRYXMLCONVERTER_HPP
