// @HEADER
// *****************************************************************************
//                    Teuchos: Common Tools Package
//
// Copyright 2004 NTESS and the Teuchos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
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
class TEUCHOSPARAMETERLIST_LIB_DLL_EXPORT XMLParameterListWriter {

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
