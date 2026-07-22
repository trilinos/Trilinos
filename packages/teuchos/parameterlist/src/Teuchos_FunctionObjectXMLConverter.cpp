// @HEADER
// *****************************************************************************
//                    Teuchos: Common Tools Package
//
// Copyright 2004 NTESS and the Teuchos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Teuchos_FunctionObjectXMLConverter.hpp"

namespace Teuchos{


RCP<FunctionObject>
FunctionObjectXMLConverter::fromXMLtoFunctionObject(
  const XMLObject &xmlObj) const
{
  return convertXML(xmlObj);
}

XMLObject
FunctionObjectXMLConverter::fromFunctionObjecttoXML(
  const RCP<const FunctionObject> function) const
{
  XMLObject toReturn(FunctionObject::getXMLTagName());
  toReturn.addAttribute(
    getTypeAttributeName(), function->getTypeAttributeValue());
  convertFunctionObject(function, toReturn);
  return toReturn;
}


} //namespace teuchos
