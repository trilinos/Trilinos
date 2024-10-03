// @HEADER
// *****************************************************************************
//                    Teuchos: Common Tools Package
//
// Copyright 2004 NTESS and the Teuchos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Teuchos_XMLParameterListCoreHelpers.hpp"
#include "Teuchos_FileInputSource.hpp"
#include "Teuchos_StringInputSource.hpp"
#include "Teuchos_XMLParameterListReader.hpp"
#include "Teuchos_XMLParameterListWriter.hpp"
#include <fstream>

void Teuchos::updateParametersFromXmlFile(
  const std::string &xmlFileName,
  const Ptr<ParameterList> &paramList
  )
{
  XMLParameterListReader xmlPLReader;
  xmlPLReader.setAllowsDuplicateSublists( false );
  FileInputSource xmlFile(xmlFileName);
  XMLObject xmlParams = xmlFile.getObject();
  paramList->setParameters(xmlPLReader.toParameterList(xmlParams));
}


Teuchos::RCP<Teuchos::ParameterList>
Teuchos::getParametersFromXmlFile(const std::string &xmlFileName)
{
  RCP<ParameterList> pl = parameterList();
  updateParametersFromXmlFile(xmlFileName, pl.ptr());
  return pl;
}


Teuchos::RCP<Teuchos::ParameterList>
Teuchos::getParametersFromXmlFile(
  const std::string &xmlFileName,
  RCP<DependencySheet> depSheet)
{
  XMLParameterListReader xmlPLReader;
  xmlPLReader.setAllowsDuplicateSublists( false );
  FileInputSource xmlFile(xmlFileName);
  XMLObject xmlParams = xmlFile.getObject();
  return xmlPLReader.toParameterList(xmlParams, depSheet);
}


void Teuchos::updateParametersFromXmlString(
  const std::string &xmlStr,
  const Ptr<ParameterList> &paramList,
  bool overwrite
  )
{
  XMLParameterListReader xmlPLReader;
  xmlPLReader.setAllowsDuplicateSublists( false );
  StringInputSource xmlStrSrc(xmlStr);
  XMLObject xmlParams = xmlStrSrc.getObject();
  if(overwrite) paramList->setParameters(xmlPLReader.toParameterList(xmlParams));
  else paramList->setParametersNotAlreadySet(xmlPLReader.toParameterList(xmlParams));
}


Teuchos::RCP<Teuchos::ParameterList>
Teuchos::getParametersFromXmlString(const std::string &xmlStr)
{
  RCP<ParameterList> pl = parameterList();
  updateParametersFromXmlString(xmlStr, pl.ptr());
  return pl;
}


Teuchos::RCP<Teuchos::ParameterList>
Teuchos::getParametersFromXmlString(const std::string &xmlStr,
  RCP<DependencySheet> depSheet)
{
  XMLParameterListReader xmlPLReader;
  xmlPLReader.setAllowsDuplicateSublists( false );
  StringInputSource xmlStrSrc(xmlStr);
  XMLObject xmlParams = xmlStrSrc.getObject();
  return xmlPLReader.toParameterList(xmlParams, depSheet);
}


void Teuchos::writeParameterListToXmlOStream(
  const ParameterList &paramList,
  std::ostream &xmlOut,
  RCP<const DependencySheet> depSheet
  )
{
  XMLParameterListWriter plWriter;
  XMLObject xml = plWriter.toXML(paramList, depSheet);
  xmlOut << xml << std::endl;
}


void Teuchos::writeParameterListToXmlFile(
  const ParameterList &paramList,
  const std::string &xmlFileName,
  RCP<const DependencySheet> depSheet
  )
{
  std::ofstream ofs(xmlFileName.c_str());
  writeParameterListToXmlOStream(paramList,ofs, depSheet);
}
