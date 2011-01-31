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


namespace Teuchos {


RCP<ParameterEntryValidator>
ValidatorFactory::createValidator(ValidatorType valiType)
{
  switch(valiType){
    case Int:
      return rcp(new EnhancedNumberValidator<int>());
      break;
    case Short:
      return rcp(new EnhancedNumberValidator<short>());
      break;
    case Double:
      return rcp(new EnhancedNumberValidator<double>());
      break;
    case Float:
      return rcp(new EnhancedNumberValidator<float>());
      break;
    case IntArray:
      return rcp(new ArrayNumberValidator<int>(rcp(new EnhancedNumberValidator<int>())));
      break;
    case ShortArray:
      return rcp(new ArrayNumberValidator<short>(rcp(new EnhancedNumberValidator<short>())));
      break;
    case DoubleArray:
      return rcp(new ArrayNumberValidator<double>(rcp(new EnhancedNumberValidator<double>())));
      break;
    case FloatArray:
      return rcp( new ArrayNumberValidator<float>(rcp(new EnhancedNumberValidator<float>())));
      break;
    case FileName:
      return rcp(new FileNameValidator());
      break;
    case FileNameArray:
      return rcp(new ArrayFileNameValidator(rcp(new FileNameValidator())));
      break;
      // 2010/07/30: Don't put a default clause in that does this!  This will hide defects!
  }
  // Default return!
  RCP<ParameterEntryValidator> toReturn;
  return toReturn;

}


RCP<EnhancedNumberValidator<int> > ValidatorFactory::getIntValidator(){
  return rcp(new EnhancedNumberValidator<int>());
}


RCP<EnhancedNumberValidator<short> > ValidatorFactory::getShortValidator(){
  return rcp(new EnhancedNumberValidator<short>());
}


RCP<EnhancedNumberValidator<double> > ValidatorFactory::getDoubleValidator(){
  return rcp(new EnhancedNumberValidator<double>());
}


RCP<EnhancedNumberValidator<float> > ValidatorFactory::getFloatValidator(){
  return rcp(new EnhancedNumberValidator<float>());
}


RCP<FileNameValidator> ValidatorFactory::getFileNameValidator(){
  return rcp(new FileNameValidator());
}


RCP<ArrayNumberValidator<int> > ValidatorFactory::getArrayIntValidator(){
  return rcp(new ArrayNumberValidator<int>(rcp( new EnhancedNumberValidator<int>())));
}


RCP<ArrayNumberValidator<short> > ValidatorFactory::getArrayShortValidator(){
  return rcp(new ArrayNumberValidator<short>(rcp(new EnhancedNumberValidator<short>())));
}


RCP<ArrayNumberValidator<double> > ValidatorFactory::getArrayDoubleValidator(){
  return rcp(new ArrayNumberValidator<double>(rcp(new EnhancedNumberValidator<double>())));
}


RCP<ArrayNumberValidator<float> > ValidatorFactory::getArrayFloatValidator(){
  return rcp(new ArrayNumberValidator<float>(rcp(new EnhancedNumberValidator<float>())));
}


RCP<ArrayFileNameValidator> ValidatorFactory::getArrayFileNameValidator(){
  return rcp(new ArrayFileNameValidator(rcp(new FileNameValidator())));
}


} // namespace Teuchos
