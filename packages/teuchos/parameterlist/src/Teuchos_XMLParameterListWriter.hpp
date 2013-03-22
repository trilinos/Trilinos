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
