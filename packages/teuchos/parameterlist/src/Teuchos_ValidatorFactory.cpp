// @HEADER
// *****************************************************************************
//                    Teuchos: Common Tools Package
//
// Copyright 2004 NTESS and the Teuchos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Teuchos_ValidatorFactory.hpp"


namespace Teuchos {


RCP<ParameterEntryValidator>
ValidatorFactory::createValidator(ValidatorType valiType)
{
  switch (valiType){
  case Int:
    return rcp(new EnhancedNumberValidator<int>());
  case Short:
    return rcp(new EnhancedNumberValidator<short>());
  case Double:
    return rcp(new EnhancedNumberValidator<double>());
  case Float:
    return rcp(new EnhancedNumberValidator<float>());
  case IntArray:
    return rcp(new ArrayNumberValidator<int>(rcp(new EnhancedNumberValidator<int>())));
  case ShortArray:
    return rcp(new ArrayNumberValidator<short>(rcp(new EnhancedNumberValidator<short>())));
  case DoubleArray:
    return rcp(new ArrayNumberValidator<double>(rcp(new EnhancedNumberValidator<double>())));
  case FloatArray:
    return rcp(new ArrayNumberValidator<float>(rcp(new EnhancedNumberValidator<float>())));
  case FileName:
    return rcp(new FileNameValidator());
  case FileNameArray:
    return rcp(new ArrayFileNameValidator(rcp(new FileNameValidator())));
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
