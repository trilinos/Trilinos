// @HEADER
// *****************************************************************************
//                    Teuchos: Common Tools Package
//
// Copyright 2004 NTESS and the Teuchos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Teuchos_ParameterXMLFileReader.hpp"	
#include "Teuchos_XMLParameterListReader.hpp"
#include "Teuchos_Assert.hpp"


namespace Teuchos {


ParameterXMLFileReader::ParameterXMLFileReader(const std::string& filename)
  : fis_(filename)
{;}


ParameterList ParameterXMLFileReader::getParameters() const
{
  XMLParameterListReader paramReader;
  XMLObject xml = fis_.getObject();

  return paramReader.toParameterList(xml);
}


} // namespace Teuchos
