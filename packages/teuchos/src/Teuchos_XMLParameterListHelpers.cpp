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


#include "Teuchos_XMLParameterListHelpers.hpp"
#include "Teuchos_FileInputSource.hpp"
#include "Teuchos_StringInputSource.hpp"
#include "Teuchos_XMLParameterListReader.hpp"
#include "Teuchos_XMLParameterListWriter.hpp"
#include "Teuchos_CommHelpers.hpp"


void Teuchos::updateParametersFromXmlFile(
  const std::string &xmlFileName,
  ParameterList *paramList
  )
{
  TEST_FOR_EXCEPT(paramList==NULL);
  XMLParameterListReader xmlPLReader;
  FileInputSource xmlFile(xmlFileName);
  XMLObject xmlParams = xmlFile.getObject();
  paramList->setParameters(xmlPLReader.toParameterList(xmlParams));
}


void Teuchos::updateParametersFromXmlFileAndBroadcast(
  const std::string &xmlFileName,
  ParameterList *paramList,
  const Comm<int> &comm
  )
{
  TEST_FOR_EXCEPT(paramList==NULL);
  if (comm.getSize()==1)
    updateParametersFromXmlFile(xmlFileName,paramList);
  else {
    if (comm.getRank()==0) {
      XMLParameterListReader xmlPLReader;
      FileInputSource xmlFile(xmlFileName);
      XMLObject xmlParams = xmlFile.getObject();
      std::string xmlString = toString(xmlParams);
      int strsize = xmlString.size();
      broadcast<int, int>(comm, 0, &strsize);
      broadcast<int, char>(comm, 0, strsize, &xmlString[0]);
      updateParametersFromXmlString(xmlString, paramList);
    }
    else {
      int strsize;
      broadcast<int, int>(comm, 0, &strsize);
      std::string xmlString;
      xmlString.resize(strsize);
      broadcast<int, char>(comm, 0, strsize, &xmlString[0]);
      updateParametersFromXmlString(xmlString, paramList);
    }
  }
}


Teuchos::RCP<Teuchos::ParameterList>
Teuchos::getParametersFromXmlFile( const std::string &xmlFileName )
{
  RCP<ParameterList> pl = parameterList();
  updateParametersFromXmlFile( xmlFileName, &*pl );
  return pl;
}


Teuchos::RCP<Teuchos::ParameterList>
Teuchos::getParametersFromXmlFile(
  const std::string &xmlFileName,
  RCP<DependencySheet> depSheet)
{
  XMLParameterListReader xmlPLReader;
  FileInputSource xmlFile(xmlFileName);
  XMLObject xmlParams = xmlFile.getObject();
  return xmlPLReader.toParameterList(xmlParams, depSheet);
}


void Teuchos::updateParametersFromXmlString(
  const std::string &xmlStr,
  ParameterList *paramList
  )
{
  TEST_FOR_EXCEPT(paramList==NULL);
  XMLParameterListReader xmlPLReader;
  StringInputSource xmlStrSrc(xmlStr);
  XMLObject xmlParams = xmlStrSrc.getObject();
  paramList->setParameters(xmlPLReader.toParameterList(xmlParams));
}


Teuchos::RCP<Teuchos::ParameterList>
Teuchos::getParametersFromXmlString( const std::string &xmlStr )
{
  RCP<ParameterList> pl = parameterList();
  updateParametersFromXmlString( xmlStr, &*pl );
  return pl;
}


Teuchos::RCP<Teuchos::ParameterList>
Teuchos::getParametersFromXmlString( const std::string &xmlStr,
  RCP<DependencySheet> depSheet)
{
  XMLParameterListReader xmlPLReader;
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
