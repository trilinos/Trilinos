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
class TEUCHOS_LIB_DLL_EXPORT XMLParameterListReader{

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


private:

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
