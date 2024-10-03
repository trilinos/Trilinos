// @HEADER
// *****************************************************************************
//                    Teuchos: Common Tools Package
//
// Copyright 2004 NTESS and the Teuchos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef TEUCHOS_VALIDATORFACTORY_HPP_
#define TEUCHOS_VALIDATORFACTORY_HPP_


#include "Teuchos_StandardParameterEntryValidators.hpp"


namespace Teuchos {


/** \brief Factory for ParameterEntryValidator objects.
 */
class ValidatorFactory {
public:

  enum ValidatorType{Int, Short, Double, Float, FileName,
  IntArray, ShortArray, DoubleArray, FloatArray, FileNameArray};

  /** \brief Creates a validator of the given type.
   *
   * @param valiType The type of validator to be created.
   */
  static RCP<ParameterEntryValidator> createValidator(ValidatorType valiType);

  /** \brief Creates and returns a Enhanced Number Validator of type int.   */
  static RCP<EnhancedNumberValidator<int> > getIntValidator();

  /** \brief Creates and returns a Enhanced Number Validator of type short. */
  static RCP<EnhancedNumberValidator<short> > getShortValidator();

  /** \brief Creates and returns a Enhanced Number Validator of type double. */
  static RCP<EnhancedNumberValidator<double> > getDoubleValidator();

  /** \brief Creates and returns a Enhanced Number Validator of type float. */
  static RCP<EnhancedNumberValidator<float> > getFloatValidator();

  /** \brief Creates and returns FileNameValidator. */
  static RCP<FileNameValidator> getFileNameValidator();

  /** \brief Creates and returns an Array Number Validator of type int. */
  static RCP<ArrayNumberValidator<int> > getArrayIntValidator();

  /** \brief Creates and returns an Array Number Validator of type short. */
  static RCP<ArrayNumberValidator<short> > getArrayShortValidator();

  /** \brief Creates and returns an Array Number Validator of type double. */
  static RCP<ArrayNumberValidator<double> > getArrayDoubleValidator();

  /** \brief Creates and returns an Array Number Validator of type float. */
  static RCP<ArrayNumberValidator<float> > getArrayFloatValidator();

  /** \brief Creates and returns an Array File Name Validator. */
  static RCP<ArrayFileNameValidator> getArrayFileNameValidator();

};


} // namespace Teuchos


#endif /* TEUCHOS_VALIDATORFACTORY_HPP_ */
