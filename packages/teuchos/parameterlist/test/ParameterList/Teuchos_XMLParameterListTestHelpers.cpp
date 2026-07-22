// @HEADER
// *****************************************************************************
//                    Teuchos: Common Tools Package
//
// Copyright 2004 NTESS and the Teuchos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Teuchos_XMLParameterListTestHelpers.hpp"
#include "Teuchos_XMLParameterListCoreHelpers.hpp"


Teuchos::RCP<Teuchos::ParameterList>
Teuchos::writeThenReadPL(ParameterList& myList)
{
  std::ostringstream xmlOut;
  writeParameterListToXmlOStream(myList, xmlOut);
  return getParametersFromXmlString(xmlOut.str());
}


Teuchos::RCP<Teuchos::ParameterList>
Teuchos::writeThenReadPL(
  ParameterList& myList,
  RCP<DependencySheet> depSheetIn,
  RCP<DependencySheet> depSheetOut)
{
  std::ostringstream xmlOut;
  writeParameterListToXmlOStream(myList, xmlOut, depSheetIn);
  return getParametersFromXmlString(xmlOut.str(), depSheetOut);
}
