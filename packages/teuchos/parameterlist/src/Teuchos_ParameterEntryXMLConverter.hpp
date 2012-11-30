// @HEADER
// ***********************************************************************
//
//                    Teuchos: Common Tools Package
//                 Copyright (2004) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
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
