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
class TEUCHOS_LIB_DLL_EXPORT ParameterEntryXMLConverter : public Describable {

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
