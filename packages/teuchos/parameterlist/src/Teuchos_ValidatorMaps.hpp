// @HEADER
// *****************************************************************************
//                    Teuchos: Common Tools Package
//
// Copyright 2004 NTESS and the Teuchos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef Teuchos_VALIDATORMAPS_HPP
#define Teuchos_VALIDATORMAPS_HPP

/*! \file Teuchos_ValidatorMaps.hpp"
 */


#include "Teuchos_ParameterEntryValidator.hpp"


namespace Teuchos {


/** \brief Maps Validators to integers. */
class TEUCHOSPARAMETERLIST_LIB_DLL_EXPORT IDtoValidatorMap {
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
class TEUCHOSPARAMETERLIST_LIB_DLL_EXPORT ValidatortoIDMap {
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
