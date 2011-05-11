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

#ifndef TEUCHOS_DEPENDENCYXMLCONVERTER_HPP
#define TEUCHOS_DEPENDENCYXMLCONVERTER_HPP

/*! \file Teuchos_DependencyXMLConverter.hpp
 * \brief Converts back and forth between XML
 * and Dependencies.
*/

#include "Teuchos_Dependency.hpp"
#include "Teuchos_XMLObject.hpp"
#include "Teuchos_XMLParameterListWriter.hpp"
#include "Teuchos_XMLParameterListReader.hpp"
#include "Teuchos_Describable.hpp"


namespace Teuchos {


/** \brief An abstract base class for converting Dependencies to
 * and from XML.
 */
class TEUCHOS_LIB_DLL_EXPORT DependencyXMLConverter : public Describable {

public:

  /** \name Converter Functions */
  //@{
  
  /** \brief Converts a given XMLObject to a Dependency.
   *
   * @param xmlObj The XMLObject to convert to a Dependency.
   * @param entryIDsMap A map containing ParameterEntrys and their
   * associated IDs.
   * @param validtorIDsMap A map containing ParameterEntryValidators and their
   * associated IDs.
   * @return The converted Dependency.
   */
  RCP<Dependency> fromXMLtoDependency(
    const XMLObject& xmlObj,
    const XMLParameterListReader::EntryIDsMap& entryIDsMap,
    const IDtoValidatorMap& validatorIDsMap) const;

  /** \brief Preforms any and all special xml conversion that is specific to a
   * particular Dependency.
   *
   * @param xmlObj The xml to be converted.
   * in which this resulting dependency will be inserted.
   * @param dependees The dependees of the dependency.
   * @param dependents The dependents of the dependency.
   * @param entryIDsMap A map containing ParameterEntrys and their
   * associated IDs.
   * @param validtorIDsMap A map containing ParameterEntryValidators and their
   * associated IDs.
   * @return The converted Dependency.
   */
  virtual RCP<Dependency> convertXML(
    const XMLObject& xmlObj, 
    const Dependency::ConstParameterEntryList dependees,
    const Dependency::ParameterEntryList dependets,
    const XMLParameterListReader::EntryIDsMap& entryIDsMap,
    const IDtoValidatorMap& validatorIDsMap) const = 0;

  /** \brief Converters a given ParameterEntryValidator to XML.
   *
   * @param dependency The Dependency to be converted to XML.
   * @param entryIDsMap A map containing ParameterEntrys and their
   * associated IDs.
   * @param validtorIDsMap A map containing ParameterEntryValidators and their
   * associated IDs.
   * @return An XML representation of the given Dependency.
   */
  XMLObject fromDependencytoXML(
    const RCP<const Dependency> dependency,
    const XMLParameterListWriter::EntryIDsMap& entryIDsMap,
    ValidatortoIDMap& validatorIDsMap) const;
  
  /** \brief Preforms any and all special dependency conversion that is
   * specific to a particlar Dependency.
   *
   * @param dependency The validator to be converted.
   * @param xmlObj The XMLObject on which all child tags should be appended
   * and attributes added.
   * @param entryIDsMap A map containing ParameterEntrys and their
   * associated IDs.
   * @param validtorIDsMap A map containing ParameterEntryValidators and their
   * associated IDs.
   * @return An XML representation of the given Dependency.
   */
  virtual void convertDependency(
    const RCP<const Dependency> dependency, 
    XMLObject& xmlObj,
    const XMLParameterListWriter::EntryIDsMap& entryIDsMap,
    ValidatortoIDMap& validatorIDsMap) const = 0;
  
  //@}

  //! \name Attribute/Query Functions
  //@{

  /**
   * \brief Returns the string to be used for the dependee tag.
   */
  static const std::string& getDependeeTagName(){
    static const std::string dependeeTagName = "Dependee";
    return dependeeTagName;
  }

  /**
   * \brief Returns the string to be used for the dependent tag.
   */
  static const std::string& getDependentTagName(){
    static const std::string dependentTagName = "Dependent";
    return dependentTagName;
  }

  /**
   * \brief Returns the string to be used for the ParameterID attribute.
   */
  static const std::string& getParameterIdAttributeName(){
    static const std::string parameterIdAtrributeName = "parameterId";
    return parameterIdAtrributeName;
  }
 
  /**
   * \brief Returns the string to be used for the type attribute.
   */
  static const std::string& getTypeAttributeName(){
    static const std::string typeAttributeName = "type";
    return typeAttributeName;
  }
 
  //@}

};


} // namespace Teuchos


#endif // TEUCHOS_DEPENDENCYXMLCONVERTER_HPP
