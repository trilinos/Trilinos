// @HEADER
// *****************************************************************************
//                    Teuchos: Common Tools Package
//
// Copyright 2004 NTESS and the Teuchos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Teuchos_RawParameterListHelpers.hpp"
#include "Teuchos_XMLParameterListWriter.hpp"
#include "Teuchos_XMLParameterListHelpers.hpp"
#include "Teuchos_CommHelpers.hpp"

TEUCHOSCOMM_LIB_DLL_EXPORT
void Teuchos::updateParametersAndBroadcast(
  const Ptr<ParameterList> &inParamList,                                     
  const Ptr<ParameterList> &paramList,
  const Comm<int> &comm,
  int root,
  bool overwrite)
{
  
  if (comm.getSize()==1) {
    if(overwrite) paramList->setParameters(*inParamList);
    else paramList->setParametersNotAlreadySet(*inParamList);
  }
  else {
    if (comm.getRank()==root) {
      XMLParameterListWriter w;
      std::string xmlString = toString(w.toXML(*inParamList));
      int strsize = static_cast<int>(xmlString.size());
      broadcast<int, int>(comm, root, &strsize);
      broadcast<int, char>(comm, root, strsize, &xmlString[0]);
      updateParametersFromXmlString(xmlString, paramList,overwrite);
    }
    else {
      int strsize;
      broadcast<int, int>(comm, root, &strsize);
      std::string xmlString;
      xmlString.resize(strsize);
      broadcast<int, char>(comm, root, strsize, &xmlString[0]);
      updateParametersFromXmlString(xmlString, paramList,overwrite);
    }
  }
}
