// @HEADER
// *****************************************************************************
//            LOCA: Library of Continuation Algorithms Package
//
// Copyright 2001-2005 NTESS and the LOCA contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "IOContFileUtils.H"

#include <fstream>
#include <iostream>
#include <sstream>

using namespace std;

bool WriteHeaderToContFile( const std::string & fileName,
    const Teuchos::ParameterList & fileParams )
{
  // The file to open
  std::ofstream oFile(fileName.c_str());

  // Writing the header
  oFile << "#" << setw(6) << "ID";

  // Looping on the parameters
  Teuchos::ParameterList::ConstIterator i;
  for (i = fileParams.begin(); i !=fileParams.end(); ++i)
    oFile << setw(15) << fileParams.name(i);
  oFile << std::endl;

  // Closing
  oFile.close();

  return true;
}

bool UpdateContFile( const std::string & fileName,
    const int & idStep,
    const Teuchos::ParameterList & fileParams )
{

  // The file to open
  std::ofstream oFile(fileName.c_str(), ios_base::app);

  // Writing the id
  oFile << scientific << setw(7)  << idStep;

  // Looping on the parameters
  Teuchos::ParameterList::ConstIterator i;
  for (i = fileParams.begin(); i !=fileParams.end(); ++i)
    oFile << scientific << setw(15) << fileParams.entry(i);
  oFile << std::endl;

  // Closing
  oFile.close();

  return true;
}

bool RestartContFile( const std::string & fileName, const int & idStep )
{
  // Copying the continuation in a backup file
  std::string fileNameBak = fileName + ".bak";
  std::string command = "cp " + fileName + " " + fileNameBak;
  TEUCHOS_ASSERT_EQUALITY(0, system(command.c_str()));

  // String of the line to cut from
  std::ostringstream os;
  os << idStep + 1;
  std::string lineNumber = os.str();

  // Cutting the file
  command =  "sed '" + lineNumber + ",$ d' " + fileNameBak + " > " + fileName;
  TEUCHOS_ASSERT_EQUALITY(0, system(command.c_str()));

  return true;
}

bool TouchContFileParameters( Teuchos::ParameterList & fileParams )
{
  // Either int or double type

  // Looping the list
  Teuchos::ParameterList::ConstIterator i;
  for (i = fileParams.begin(); i !=fileParams.end(); ++i) {

    if (fileParams.isType<int>(fileParams.name(i)))
      fileParams.get<int>(fileParams.name(i));

    if (fileParams.isType<double>(fileParams.name(i)))
      fileParams.get<double>(fileParams.name(i));
  }

  return true;
}
