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


#ifndef TEUCHOS_CONDITIONXMLCONVERTERDB_HPP
#define TEUCHOS_CONDITIONXMLCONVERTERDB_HPP

/*! \file Teuchos_ConditionXMLConverterDB.hpp
 * \brief A database for ConditionXMLConverters.
*/

// Both includes needed for convience macros below 
#include "Teuchos_StandardConditionXMLConverters.hpp"
#include "Teuchos_StandardConditions.hpp"


namespace Teuchos {

class Condition;

/** \brief Provides ability to lookup ConditionXMLConverters
 */
class TEUCHOS_LIB_DLL_EXPORT ConditionXMLConverterDB {

public:

  /** \name Modifier Functions */
  //@{
  
  /** \brief Add a converter to the database.
   *
   * \param condition A dummy condition representing the type of contidion
   * this converter will convert.
   * \param convertToAdd The converter to add to the database.
   */
  static void addConverter(RCP<const Condition> condition,
    RCP<ConditionXMLConverter> converterToAdd);
  
  //@}

  /** \name Converter Functions */
  //@{
  
  /** \brief Get an appropriate ConditionXMLConverter given a 
   *  Condition.
   *
   * \param condition The Condition for which a converter is
   * desired.
   *
   * \return A converter for the condition.
   */
  static RCP<const ConditionXMLConverter> 
    getConverter(const Condition& condition);

  /** \brief Get an appropriate ConditionXMLConverter given a XMLObject.
   *
   * @param xmlObject The XMLObject for which a converter is desired.
   *
   * @return A converter for the XMLObject.
   */
  static RCP<const ConditionXMLConverter> 
    getConverter(const XMLObject& xmlObject);

  /**
   * \brief Given a condition and ConditiontoIDMap, converts the
   * condition to XML.
   *
   * \param The condition to convert.
   * \param entryIDsMap A map containing ParameterEntrys and their associated
   * IDs.
   *
   * \return XML representation of the condition.
   */
  static XMLObject convertCondition(
    RCP<const Condition> condition,
    const XMLParameterListWriter::EntryIDsMap& entryIDsMap);

  /**
   * \brief Given an XMLObject and IDtoConditionMap, converts the XMLObject 
   * to a Condition.
   *
   * \param xmlObject The XMLObject to be converted.
   * \param entryIDsMap A map containing ParameterEntrys and their associated
   * IDs.
   *
   * \return A Condition that was represented by the XML.
   */
  static RCP<Condition> convertXML(
    const XMLObject& xmlObject,
    const XMLParameterListReader::EntryIDsMap& entryIDsMap);
  
  //@}

  /** \name I/O Functions */
  //@{

  /**
   * \brief prints the xml tags associated with all known converters
   *
   * \param out Stream to which tags should be printed.
   */
  static void printKnownConverters(std::ostream& out){
    out << "Known ConditionXMLConverters: " << std::endl;
    for(
      ConverterMap::const_iterator it = getConverterMap().begin();
      it != getConverterMap().end();
      ++it)
    {
      out << "\t" << it->first <<std::endl;
    }
  }
  
  //@}

private:

  /** \name Private Members */
  //@{
  
  /** \brief convience class. */
  typedef std::map<std::string, RCP<ConditionXMLConverter> > ConverterMap;

  /** \brief convience typedef. */
  typedef std::pair<std::string, RCP<ConditionXMLConverter> > ConverterPair;

  /** \brief Gets the default converter to be used to convert
   * Conditions.
   */
  static ConverterMap& getConverterMap();
  
  //@}

};


} // end namespace Teuchos

//
// Helper Macros
//

/** \brief Adds a NumberCondition of type T */
#define TEUCHOS_ADD_NUMBERCONVERTER(T) \
  Teuchos::ConditionXMLConverterDB::addConverter( \
    Teuchos::DummyObjectGetter<Teuchos::NumberCondition< T > >:: \
      getDummyObject(), \
    Teuchos::rcp(new Teuchos::NumberConditionConverter< T >));

#endif // TEUCHOS_CONDITIONXMLCONVERTERDB_HPP
