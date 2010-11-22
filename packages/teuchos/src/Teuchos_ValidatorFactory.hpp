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
