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
