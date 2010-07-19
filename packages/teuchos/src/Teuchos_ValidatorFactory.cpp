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

#include "Teuchos_ValidatorFactory.hpp"
namespace Teuchos{

RCP<ParameterEntryValidator> ValidatorFactory::createValidator(ValidatorType valiType){
	switch(valiType){
		case Int:
			return RCP<EnhancedNumberValidator<int> >(new EnhancedNumberValidator<int>());
			break;
		case Short:
			return RCP<EnhancedNumberValidator<short> >(new EnhancedNumberValidator<short>());
			break;
		case Double:
			return RCP<EnhancedNumberValidator<double> >(new EnhancedNumberValidator<double>());
			break;
		case Float:
			return RCP<EnhancedNumberValidator<float> >(new EnhancedNumberValidator<float>());
			break;
		case IntArray:
			return RCP<ArrayNumberValidator<int> >( new ArrayNumberValidator<int>( RCP<EnhancedNumberValidator<int> >( new EnhancedNumberValidator<int>())));
			break;
		case ShortArray:
			return RCP<ArrayNumberValidator<short> >( new ArrayNumberValidator<short>( RCP<EnhancedNumberValidator<short> >( new EnhancedNumberValidator<short>())));
			break;
		case DoubleArray:
			return RCP<ArrayNumberValidator<double> >( new ArrayNumberValidator<double>( RCP<EnhancedNumberValidator<double> >( new EnhancedNumberValidator<double>())));
			break;
		case FloatArray:
			return RCP<ArrayNumberValidator<float> >( new ArrayNumberValidator<float>(RCP<EnhancedNumberValidator<float> >( new EnhancedNumberValidator<float>())));
			break;
		case FileName:
			return RCP<FileNameValidator>(new FileNameValidator());
			break;
		case FileNameArray:
			return RCP<ArrayFileNameValidator>(new ArrayFileNameValidator(RCP<FileNameValidator>(new FileNameValidator())));
			break;
		default:
			RCP<ParameterEntryValidator> toReturn;
			return toReturn;
			break;
	}
	RCP<ParameterEntryValidator> toReturn;
	return toReturn;
}

RCP<EnhancedNumberValidator<int> > ValidatorFactory::getIntValidator(){
	return RCP<EnhancedNumberValidator<int> >(new EnhancedNumberValidator<int>());
}

RCP<EnhancedNumberValidator<short> > ValidatorFactory::getShortValidator(){
	return RCP<EnhancedNumberValidator<short> >(new EnhancedNumberValidator<short>());
}

RCP<EnhancedNumberValidator<double> > ValidatorFactory::getDoubleValidator(){
	return RCP<EnhancedNumberValidator<double> >(new EnhancedNumberValidator<double>());
}

RCP<EnhancedNumberValidator<float> > ValidatorFactory::getFloatValidator(){
	return RCP<EnhancedNumberValidator<float> >(new EnhancedNumberValidator<float>());
}

RCP<FileNameValidator> ValidatorFactory::getFileNameValidator(){
	return RCP<FileNameValidator>(new FileNameValidator());
}

RCP<ArrayNumberValidator<int> > ValidatorFactory::getArrayIntValidator(){
	return RCP<ArrayNumberValidator<int> >( new ArrayNumberValidator<int>( RCP<EnhancedNumberValidator<int> >( new EnhancedNumberValidator<int>())));
}

RCP<ArrayNumberValidator<short> > ValidatorFactory::getArrayShortValidator(){
	return RCP<ArrayNumberValidator<short> >( new ArrayNumberValidator<short>( RCP<EnhancedNumberValidator<short> >( new EnhancedNumberValidator<short>())));
}

RCP<ArrayNumberValidator<double> > ValidatorFactory::getArrayDoubleValidator(){
	return RCP<ArrayNumberValidator<double> >( new ArrayNumberValidator<double>( RCP<EnhancedNumberValidator<double> >( new EnhancedNumberValidator<double>())));
}

RCP<ArrayNumberValidator<float> > ValidatorFactory::getArrayFloatValidator(){
	return RCP<ArrayNumberValidator<float> >( new ArrayNumberValidator<float>( RCP<EnhancedNumberValidator<float> >( new EnhancedNumberValidator<float>())));
}

RCP<ArrayFileNameValidator> ValidatorFactory::getArrayFileNameValidator(){
	return RCP<ArrayFileNameValidator>(new ArrayFileNameValidator(RCP<FileNameValidator>(new FileNameValidator())));
}

}
