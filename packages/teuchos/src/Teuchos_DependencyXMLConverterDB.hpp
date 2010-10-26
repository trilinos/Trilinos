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


#ifndef TEUCHOS_DEPENDENCYXMLCONVERTERDB_HPP
#define TEUCHOS_DEPENDENCYXMLCONVERTERDB_HPP

/*! \file Teuchos_DependencyXMLConverterDB.hpp
 * \brief A database for DependencyXMLConverters.
*/

#include "Teuchos_DependencyXMLConverter.hpp"
#include "Teuchos_XMLParameterListReader.hpp"


namespace Teuchos {

class Dependency;

/** \brief Provides ability to lookup DependencyXMLConverterDB
 */
class DependencyXMLConverterDB {
public:

  /** \name Modifier Functions */
  //@{
  
  /** \brief Add a converter to the database.
   *
   * \param convertToAdd The converter to add to the database.
   */
  static void addConverter(Dependency& dependency,
    RCP<DependencyXMLConverter> converterToAdd);
  
  //@}

  /** \name Converter Functions */
  //@{
  
  /** \brief Get an appropriate DependencyXMLConverter given a 
   *  ParameterEntry.
   *
   * \param dependency The ParameterEntryDependency for which a 
   * converter is desired.
   */
  static RCP<const DependencyXMLConverter> getConverter(
    const Dependency& dependency);

  /** \brief Get an appropriate DependencyXMLConverter given a XMLObject.
   *
   * @param xmlObject The XMLObject for which a converter is desired.
   */
  static RCP<const DependencyXMLConverter> 
    getConverter(const XMLObject& xmlObject);

  /**
   * \brief Given a dependency converts the
   * dependency to XML.
   *
   * \return XML representation of the dependency.
   */
  static XMLObject convertDependency(
    RCP<const Dependency> dependency,
    const XMLParameterListWriter::EntryIDsMap& entryIDsMap,
    ValidatortoIDMap& validatorIDsMap); 

  /**
   * \brief Given an XMLObject converts the XMLObject 
   * to a Dependency.
   *
   * \return A Dependency that was represented by the XML.
   */
  static RCP<Dependency> convertXML(
    const XMLObject& xmlObject, 
    const XMLParameterListReader::EntryIDsMap& entryIDsMap,
    const IDtoValidatorMap& validatorIDsMap); 
  
  //@}

  /** \name I/O Functions */
  //@{

  /**
   * \brief prints the xml tags associated with all known converters
   *
   * \param out Stream to which tags should be printed.
   */
  static void printKnownConverters(std::ostream& out){
    out << "Known DependencyXMLConverters: " << std::endl;
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
  typedef std::map<std::string, RCP<DependencyXMLConverter> > ConverterMap;

  /** \brief convience typedef. */
  typedef std::pair<std::string, RCP<DependencyXMLConverter> > 
    ConverterPair;

  /** \brief Gets the default converter to be used to convert
   * Dependencies.
   */
  static ConverterMap& getConverterMap();
  
  //@}

};


} // end namespace Teuchos


#endif // TEUCHOS_DEPENDENCYXMLCONVERTERDB_HPP
