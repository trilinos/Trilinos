// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ROL_PARAMETERENTRYXMLCONVERTER_HPP
#define ROL_PARAMETERENTRYXMLCONVERTER_HPP


/*! \file ROL_ParameterEntryXMLCoverter.hpp
 *  \brief The base class for all ParameterEntryXMLConverters.
*/

#include "ROL_XMLObject.hpp"
#include "Teuchos_ParameterEntry.hpp"


namespace ROL {

/** \brief A class used to convert parameter entries to xml and vice versa.
 */
class TEUCHOSPARAMETERLIST_LIB_DLL_EXPORT ParameterEntryXMLConverter {

public:

  /* Destructor */
  virtual ~ParameterEntryXMLConverter() {}

  /** \name Converter Functions */
  //@{

  /** \brief Converts the given xml into a parameter entry.
   *
   * \param xmlObj The xml to be converted to a parameter entry.
   * \returns A ParameterEntry with the aspects specified by the xml.
   */
  ParameterEntry fromXMLtoParameterEntry(const XMLObject &xmlObj) const;

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
    Teuchos::RCP<const ParameterEntry > entry) const=0;

  /** \brief . */
  static const std::string& getTypeAttributeName() {
    static const std::string typeAttributeName_ = "type";
    return typeAttributeName_;
  }

  /** \brief . */
  static const std::string& getNameAttributeName(){
    static const std::string nameAttributeName = "name";
    return nameAttributeName;
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


} // namespace ROL


#endif // TEUCHOS_PARAMETERENTRYXMLCONVERTER_HPP
