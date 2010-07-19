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


namespace Teuchos{

class ValidatorFactory{
public:
	enum ValidatorType{Int, Short, Double, Float, FileName,
	IntArray, ShortArray, DoubleArray, FloatArray, FileNameArray};

	/**
	 * Creates a validator of the given type.
	 * 
	 * @param valiType The type of validator to be created.
	 * @return A validator of the specified type.
	 */
	static RCP<ParameterEntryValidator> createValidator(ValidatorType valiType);

	/**
	 * Creates and returns a Enhanced Number Validator of type int.
	 *
	 * @return An Enhanced Number Validator of type int.
	 */
	static RCP<EnhancedNumberValidator<int> > getIntValidator();

	/**
	 * Creates and returns a Enhanced Number Validator of type short.
	 *
	 * @return An Enhanced Number Validator of type short.
	 */
	static RCP<EnhancedNumberValidator<short> > getShortValidator();

	/**
	 * Creates and returns a Enhanced Number Validator of type double.
	 *
	 * @return An Enhanced Number Validator of type double.
	 */
	static RCP<EnhancedNumberValidator<double> > getDoubleValidator();

	/**
	 * Creates and returns a Enhanced Number Validator of type float.
	 *
	 * @return An Enhanced Number Validator of type float.
	 */
	static RCP<EnhancedNumberValidator<float> > getFloatValidator();

	/**
	 * Creates and returns FileNameValidator.
	 *
	 * @return A FileNameValidator.
	 */
	static RCP<FileNameValidator> getFileNameValidator();

	/**
	 * Creates and returns an Array Number Validator of type int.
	 *
	 * @return An Enhanced Number Validator of type int.
	 */
	static RCP<ArrayNumberValidator<int> > getArrayIntValidator();

	/**
	 * Creates and returns an Array Number Validator of type short.
	 *
	 * @return An Enhanced Number Validator of type short.
	 */
	static RCP<ArrayNumberValidator<short> > getArrayShortValidator();

	/**
	 * Creates and returns an Array Number Validator of type double.
	 *
	 * @return An Enhanced Number Validator of type double.
	 */
	static RCP<ArrayNumberValidator<double> > getArrayDoubleValidator();

	/**
	 * Creates and returns an Array Number Validator of type float.
	 *
	 * @return An Enhanced Number Validator of type float.
	 */
	static RCP<ArrayNumberValidator<float> > getArrayFloatValidator();

	/**
	 * Creates and returns an Array File Name Validator.
	 *
	 * @return An Array File Name Validator.
	 */
	static RCP<ArrayFileNameValidator> getArrayFileNameValidator();
};

}

#endif /* TEUCHOS_VALIDATORFACTORY_HPP_ */
