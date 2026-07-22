// @HEADER
// *****************************************************************************
//                    Teuchos: Common Tools Package
//
// Copyright 2004 NTESS and the Teuchos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Teuchos_XMLParameterListHelpers.hpp"
#include "Teuchos_FileInputSource.hpp"
#include "Teuchos_XMLParameterListReader.hpp"
#include "Teuchos_CommHelpers.hpp"

void Teuchos::updateParametersFromXmlFileAndBroadcast(
  const std::string &xmlFileName,
  const Ptr<ParameterList> &paramList,
  const Comm<int> &comm,
  bool overwrite
  )
{
  if (comm.getSize()==1)
    updateParametersFromXmlFile(xmlFileName, paramList);
  else {
    if (comm.getRank()==0) {
      XMLParameterListReader xmlPLReader;
      xmlPLReader.setAllowsDuplicateSublists( false );
      FileInputSource xmlFile(xmlFileName);
      XMLObject xmlParams = xmlFile.getObject();
      std::string xmlString = toString(xmlParams);
      int strsize = static_cast<int>(xmlString.size());
      broadcast<int, int>(comm, 0, &strsize);
      broadcast<int, char>(comm, 0, strsize, &xmlString[0]);
      updateParametersFromXmlString(xmlString, paramList,overwrite);
    }
    else {
      int strsize;
      broadcast<int, int>(comm, 0, &strsize);
      std::string xmlString;
      xmlString.resize(strsize);
      broadcast<int, char>(comm, 0, strsize, &xmlString[0]);
      updateParametersFromXmlString(xmlString, paramList,overwrite);
    }
  }
}
