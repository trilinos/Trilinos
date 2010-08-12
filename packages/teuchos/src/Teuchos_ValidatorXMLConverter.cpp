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

#include "Teuchos_ValidatorXMLConverter.hpp"

namespace Teuchos{

RCP<ParameterEntryValidator>
ValidatorXMLConverter::fromXMLtoValidator(const XMLObject& xmlObj) const {
  #ifdef HAVE_TEUCHOS_DEBUG
  RCP<const ParameterEntryValidator> dummyValidator = getDummyValidator();
  TEST_FOR_EXCEPTION(xmlObj.getTag() != dummyValidator->getXMLTagName(), 
    BadValidatorXMLConverterException, 
    "Cannot convert xmlObject " << 
    ". Expected a " << dummyValidator->getXMLTagName() <<
    " tag but got a " << xmlObj.getTag() << "tag");
  #endif
  ParameterEntryValidator::ValidatorID validatorID = 
    xmlObj.getRequired<ParameterEntryValidator::ValidatorID>(
      getIdAttributeName());
  RCP<ParameterEntryValidator> toReturn = convertXML(xmlObj, validatorID);
  return toReturn;
}

XMLObject 
ValidatorXMLConverter::fromValidatortoXML(
  const RCP<const ParameterEntryValidator> validator) const
{
  XMLObject toReturn = convertValidator(validator);
  toReturn.addAttribute(getIdAttributeName(),
    ParameterEntryValidator::getValidatorID(validator));
  return toReturn;
}


}

