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

#include "Teuchos_StandardValidatorXMLConverters.hpp"
#include "Teuchos_RCP.hpp"

/*! \file Teuchos_StandardValidatorXMLConverters.hpp
*/


namespace Teuchos {


RCP<ParameterEntryValidator> AnyNumberValidatorXMLConverter::convertXML(
  const XMLObject& xmlObj,
  const IDtoValidatorMap& /*validatorIDsMap*/) const
{

  AnyNumberParameterEntryValidator::AcceptedTypes acceptedTypes;
  acceptedTypes.allowInt(xmlObj.getRequiredBool(getAllowIntAttributeName()));
  acceptedTypes.allowDouble(
    xmlObj.getRequiredBool(getAllowDoubleAttributeName()));
  acceptedTypes.allowString(
    xmlObj.getRequiredBool(getAllowStringAttributeName()));
  return anyNumberParameterEntryValidator(
      AnyNumberParameterEntryValidator::getPrefferedTypeStringEnum(
        xmlObj.getRequired(getPrefferedTypeAttributeName())), 
      acceptedTypes);
}


void AnyNumberValidatorXMLConverter::convertValidator(
  const RCP<const ParameterEntryValidator> validator,
  XMLObject& xmlObj,
  const ValidatortoIDMap& /*validatorIDsMap*/) const
{
  RCP<const AnyNumberParameterEntryValidator> castedValidator = 
    rcp_dynamic_cast<const AnyNumberParameterEntryValidator>(validator, true);
  xmlObj.addBool(
    getAllowIntAttributeName(), castedValidator->isIntAllowed());
  xmlObj.addBool(
    getAllowDoubleAttributeName(), castedValidator->isDoubleAllowed());
  xmlObj.addBool(
    getAllowStringAttributeName(), castedValidator->isStringAllowed());
  xmlObj.addAttribute(getPrefferedTypeAttributeName(),
    castedValidator->getPrefferedTypeString(
      castedValidator->getPreferredType()));
}

#ifdef HAVE_TEUCHOS_DEBUG
RCP<const ParameterEntryValidator> 
AnyNumberValidatorXMLConverter::getDummyValidator() const{
  return DummyObjectGetter<AnyNumberParameterEntryValidator>::getDummyObject();
}
#endif

RCP<ParameterEntryValidator> FileNameValidatorXMLConverter::convertXML(
  const XMLObject& xmlObj,
  const IDtoValidatorMap& /*validatorIDsMap*/) const
{
  return rcp(
    new FileNameValidator(
      xmlObj.getWithDefault<bool>(
        getFileMustExistAttributeName(),
        FileNameValidator::mustAlreadyExistDefault()
        )
      )
    );
}


void FileNameValidatorXMLConverter::convertValidator(
  const RCP<const ParameterEntryValidator> validator,
  XMLObject& xmlObj,
  const ValidatortoIDMap& /*validatorIDsMap*/) const
{
  RCP<const FileNameValidator> castedValidator =
    rcp_dynamic_cast<const FileNameValidator>(validator);
  xmlObj.addBool(
    getFileMustExistAttributeName(), castedValidator->fileMustExist());
}


#ifdef HAVE_TEUCHOS_DEBUG
RCP<const ParameterEntryValidator> 
FileNameValidatorXMLConverter::getDummyValidator() const{
  return DummyObjectGetter<FileNameValidator>::getDummyObject();
}
#endif


RCP<ParameterEntryValidator> StringValidatorXMLConverter::convertXML(
  const XMLObject& xmlObj,
  const IDtoValidatorMap& /*validatorIDsMap*/) const
{
  Array<std::string> strings(xmlObj.numChildren());
  if(xmlObj.numChildren()!=0){
    for(int i=0; i<xmlObj.numChildren(); ++i){
      XMLObject currentChild = xmlObj.getChild(i);
      TEST_FOR_EXCEPTION(currentChild.getTag() != getStringTagName(), 
        BadTagException,  
        "Error converting xmlObject to StringValidator." << std::endl << 
		    "Unrecognized tag: " << currentChild.getTag());
        strings[i] = (currentChild.getRequired(getStringValueAttributeName()));
    }
  }
  return rcp(new StringValidator(strings));
}


void StringValidatorXMLConverter::convertValidator(
  const RCP<const ParameterEntryValidator> validator,
  XMLObject& xmlObj,
  const ValidatortoIDMap& /*validatorIDsMap*/) const
{
  RCP<const StringValidator> castedValidator = 
    rcp_dynamic_cast<const StringValidator>(validator);

  if(!is_null(validator->validStringValues())){
    Array<std::string>::const_iterator it = 
     validator->validStringValues()->begin();
    for(; it != validator->validStringValues()->end(); ++it){
      XMLObject stringTag(getStringTagName());
      stringTag.addAttribute(getStringValueAttributeName(), *it);
      xmlObj.addChild(stringTag);
    }
  }
}


#ifdef HAVE_TEUCHOS_DEBUG
RCP<const ParameterEntryValidator> 
StringValidatorXMLConverter::getDummyValidator() const{
  return DummyObjectGetter<StringValidator>::getDummyObject();
}
#endif

} // namespace Teuchos

