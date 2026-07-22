// @HEADER
// *****************************************************************************
//                    Teuchos: Common Tools Package
//
// Copyright 2004 NTESS and the Teuchos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef Teuchos_XMLPARAMETERLISTREADER_H
#define Teuchos_XMLPARAMETERLISTREADER_H

/*! \file Teuchos_XMLParameterListReader.hpp
    \brief Writes an XML object to a parameter list
*/

#include "Teuchos_ParameterList.hpp"
#include "Teuchos_XMLObject.hpp"
#include "Teuchos_Utils.hpp"
#include "Teuchos_DependencySheet.hpp"
#include "Teuchos_ValidatorMaps.hpp"


namespace Teuchos {


/** \brief Writes an XML object to a parameter list.
 *
 * \ingroup XML
 */
class TEUCHOSPARAMETERLIST_LIB_DLL_EXPORT XMLParameterListReader{

public:

  /** \name Public Types */
  //@{

  /** \brief Convenience typedef */
  typedef std::map<ParameterEntry::ParameterEntryID,
    RCP<ParameterEntry> > EntryIDsMap;

  //@}

  //! @name Constructors
  //@{
  /** \brief . */
  XMLParameterListReader();
  //@}

  /** Write the given XML object to a parameter list */
  RCP<ParameterList> toParameterList(
    const XMLObject& xml, RCP<DependencySheet> depSheet) const;

  /** Write the given XML object to a parameter list */
  ParameterList toParameterList(const XMLObject& xml) const;

  /** \brief Set policy regarding duplicated sublists
    *
    * The default behavior of this class is to allow duplicated sublists,
    * although the resulting
    * ParameterList is undefined for the duplicated sublists (in most
    * cases, they will be merged in the order they are encountered in the
    * XML character stream).
    *
    * If set \c false, then duplicated sublists in the XML tree
    * will result in the Teuchos::DuplicateParameterSublist
    * exception being thrown.
    *
    * If set \c true, the default behavior is restored.
    */
  void setAllowsDuplicateSublists(bool policy);

  /** \brief Specifies the current policy regarding duplicated sublists.
      See setAllowsDuplicateSublists() for more details.
  */
  bool getAllowsDuplicateSublists() const;

private:

  bool _allowDuplicateSublists;

  /** \brief Write the given XML object to a parameter list along with the
   * validators located in the given map.
   */
  void convertParameterList(const XMLObject& xml,
    RCP<ParameterList> parentList,
    EntryIDsMap& entryIDsMap, const IDtoValidatorMap& validatorIDsMap) const;

  /** \brief Write the given XML object to appropriate validators. */
  void convertValidators(
    const XMLObject& xml, IDtoValidatorMap& validatorIDsMap) const;

  /** \brief Write the given XML object to appropriate dependencies. */
  void convertDependencies(
    RCP<DependencySheet> depSheet,
    const XMLObject& xml,
    const EntryIDsMap& entryIDsMap,
    const IDtoValidatorMap& validatorIDsMap) const;

  /** \brief Tests to see if there are duplicate validator IDs */
  void testForDuplicateValidatorIDs(
    ParameterEntryValidator::ValidatorID potentialNewID,
    const IDtoValidatorMap& currentMap) const;

  /** \brief . */
  void insertEntryIntoMap(
    const XMLObject& xmlObj,
    RCP<ParameterEntry> entryToInsert,
    EntryIDsMap& entryIDsMap) const;

};



} // namespace Teuchos


#endif
