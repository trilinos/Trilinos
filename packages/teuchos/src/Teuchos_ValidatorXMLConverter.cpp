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

#include "Teuchos_ValidatorXMLConverter.hpp"

namespace Teuchos{

RCP<ParameterEntryValidator>
ValidatorXMLConverter::fromXMLtoValidator(
  const XMLObject& xmlObj,
  const IDtoValidatorMap& validatorIDsMap) const
{
  #ifdef HAVE_TEUCHOS_DEBUG
  RCP<const ParameterEntryValidator> dummyValidator = getDummyValidator();
  TEST_FOR_EXCEPTION(
    xmlObj.getRequired(getTypeAttributeName())
    !=
    dummyValidator->getXMLTypeName(), 
    BadValidatorXMLConverterException, 
    "Cannot convert xmlObject " << 
    ". Expected a " << dummyValidator->getXMLTypeName() <<
    " tag but got a " << xmlObj.getRequired(getTypeAttributeName()) << "type");
  #endif
  RCP<ParameterEntryValidator> toReturn = 
    convertXML(xmlObj, validatorIDsMap);
  return toReturn;
}

XMLObject 
ValidatorXMLConverter::fromValidatortoXML(
  const RCP<const ParameterEntryValidator> validator,
  const ValidatortoIDMap& validatorIDsMap,
  bool assignedID) const
{
  #ifdef HAVE_TEUCHOS_DEBUG
  RCP<const ParameterEntryValidator> dummyValidator = getDummyValidator();
  TEST_FOR_EXCEPTION(
    validator->getXMLTypeName() 
    !=
    dummyValidator->getXMLTypeName(), 
    BadValidatorXMLConverterException, 
    "Cannot convert Validator " << 
    ". Expected a " << dummyValidator->getXMLTypeName() <<
    " validator but got a " << validator->getXMLTypeName() << "type");
  #endif
  XMLObject toReturn(getValidatorTagName());
  toReturn.addAttribute(getTypeAttributeName(), validator->getXMLTypeName());
  if(assignedID){
    TEST_FOR_EXCEPTION(validatorIDsMap.find(validator) == validatorIDsMap.end(),
      MissingValidatorDefinitionException,
      "Could not find an id associated with the validator in the "
      "given validatorIDsMap to use when " <<
      "writing it to XML!" << std::endl << std::endl);
    toReturn.addAttribute(getIdAttributeName(),
      validatorIDsMap.find(validator)->second);
  }
  convertValidator(validator, toReturn, validatorIDsMap);
  return toReturn;
}


}

