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
