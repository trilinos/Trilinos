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

#ifndef Teuchos_XMLPARAMETERLISTWRITER_H
#define Teuchos_XMLPARAMETERLISTWRITER_H


/*! \file Teuchos_XMLParameterListWriter.hpp
 *
 *   \brief Writes a ParameterList to an XML object
 */

#include "Teuchos_ParameterList.hpp"
#include "Teuchos_XMLObject.hpp"
#include "Teuchos_Utils.hpp"
#include "Teuchos_DependencySheet.hpp"
#include "Teuchos_ValidatorMaps.hpp"


namespace Teuchos {


/** \ingroup XML 
 * \brief Writes a ParameterList to an XML object
 */
class TEUCHOS_LIB_DLL_EXPORT XMLParameterListWriter {

public:

  /** \name Public Types */
  //@{

  /** \brief . */
  typedef std::map<RCP<const ParameterEntry>,
    ParameterEntry::ParameterEntryID, RCPConstComp> EntryIDsMap;
  
  //@}

  //! @name Constructors 
  //@{
  /** Construct a writer */
  XMLParameterListWriter();
  //@}

  /** Write the given list to an XML object */
  XMLObject toXML(
    const ParameterList& p, 
    RCP<const DependencySheet> depSheet = null) const;

  /** \brief . */
  static const std::string& getParameterListTagName(){
    static const std::string parameterListTagName = "ParameterList";
    return parameterListTagName;
  }

  /** \brief . */
  static const std::string& getNameAttributeName(){
    static const std::string nameAttributeName = "name";
    return nameAttributeName;
  } 

  /** \brief . */
  static const std::string& getValidatorsTagName(){
    static const std::string validatorsTagName = "Validators";
    return validatorsTagName;
  }

  /** \brief . */
  static const std::string& getDependenciesTagName(){
    static const std::string dependenciesTagName = "Dependencies";
    return dependenciesTagName;
  }

private:

  /** \brief Write the given list to an XML object.  */
  XMLObject convertParameterList(
      const ParameterList& p, 
      ParameterEntry::ParameterEntryID& idCounter, 
      EntryIDsMap& entryIDsMap,
      const ValidatortoIDMap& validatorIDsMap) const;

  /** \brief Convert all the validators. */
  XMLObject convertValidators(
    const ParameterList& p, 
    ValidatortoIDMap& validatorIDsMap) const;

  /** \brief Convert all the dependencies. */
  XMLObject convertDependencies(
    RCP<const DependencySheet> depSheet,
    const EntryIDsMap& entryIDsMap, 
    ValidatortoIDMap& validatorIDsMap) const;

  /** \brief Builds up the list of validators to be converted */
  void buildInitialValidatorMap(
    const ParameterList& p,
    ValidatortoIDMap& validatorIDsMap) const;
};


} // namespace teuchos


#endif
