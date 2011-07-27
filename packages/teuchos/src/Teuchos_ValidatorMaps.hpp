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

#ifndef Teuchos_VALIDATORMAPS_HPP
#define Teuchos_VALIDATORMAPS_HPP

/*! \file Teuchos_ValidatorMaps.hpp"
 */


#include "Teuchos_ParameterEntryValidator.hpp"
#include "Teuchos_ValidatorMaps.hpp"


namespace Teuchos {


/** \brief Maps Validators to integers. */
class TEUCHOS_LIB_DLL_EXPORT IDtoValidatorMap {
public:

  /** \brief . */
  typedef std::map<ParameterEntryValidator::ValidatorID,
    RCP<ParameterEntryValidator> > ValidatorMap;

  /** \brief . */
  typedef std::pair<ParameterEntryValidator::ValidatorID, 
    RCP<ParameterEntryValidator> > IDValidatorPair;

  /** \brief . */
  typedef ValidatorMap::iterator iterator;

  /** \brief . */
  typedef ValidatorMap::const_iterator const_iterator;

  /** \brief inserts an IDValidatorPair into the map. */
  void insert(IDValidatorPair toInsert);

  /** \brief Retrieves and iterator to a validator and id based on the id given.
   *
   * If no validator is found that has been mappend to the given id, a
   * reference to the end of the map is returned.
   */
  const_iterator find(int id) const;

  /** \brief Returns a const_reference to the beginning of the map. */
  const_iterator begin() const;

  /** \brief Returns a const_reference to the end of the map. */
  const_iterator end() const;

  /** \brief removes the specified validator from the map. */
  inline
  size_t erase(const ParameterEntryValidator::ValidatorID& x){
    return validatorMap.erase(x);
  }

private:

  ValidatorMap validatorMap;

};


/** \brief A class for mapping validators to integers. */
class TEUCHOS_LIB_DLL_EXPORT ValidatortoIDMap {
public:
  
  /** \brief . */
  typedef std::map<RCP<const ParameterEntryValidator>, int, RCPComp> 
    ValidatorMap;

  /** \brief . */
  typedef std::pair<RCP<const ParameterEntryValidator>, int> ValidatorIDPair;

  /** \brief . */
  typedef ValidatorMap::iterator iterator;

  /** \brief . */
  typedef ValidatorMap::const_iterator const_iterator;

  /** \brief . */
  ValidatortoIDMap();

  /** \brief inserts an IDValidatorPair into the map. */
  void insert(RCP<const ParameterEntryValidator> toInsert);
  
  /** \brief Returns an iterator to the validator and id specified by the validator.
   *
   * If no id is found with the associated validator, a reference to the end
   * of the map is returned.
   */
  const_iterator find(
    const RCP<const ParameterEntryValidator> validator) const;
  
  /** \brief Returns a const_reference to the beginning of the map. */
  const_iterator begin() const;

  /** \brief Returns a const_reference to the end of the map. */
  const_iterator end() const;

private:

  ValidatorMap validatorMap;

  int counter;

};


} // namespace Teuchos


#endif //Teuchos_VALIDATORMAPS_HPP
