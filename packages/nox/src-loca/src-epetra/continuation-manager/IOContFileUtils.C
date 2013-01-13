/*
//@HEADER
// ************************************************************************
// 
//            LOCA: Library of Continuation Algorithms Package
//                 Copyright (2005) Sandia Corporation
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
// Questions? Contact Roger Pawlowski (rppawlo@sandia.gov) or 
// Eric Phipps (etphipp@sandia.gov), Sandia National Laboratories.
// ************************************************************************
//  CVS Information
//  $Source$
//  $Author$
//  $Date$
//  $Revision$
// ************************************************************************
//@HEADER
*/

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
  system (command.c_str());

  // String of the line to cut from
  std::ostringstream os;
  os << idStep + 1; 
  std::string lineNumber = os.str();

  // Cutting the file
  command =  "sed '" + lineNumber + ",$ d' " + fileNameBak + " > " + fileName;
  system (command.c_str());

  return true;
}

bool TouchContFileParameters( Teuchos::ParameterList & fileParams )
{
  // Either int or double type
  int dummyInt;
  double dummyDouble;

  // Looping the list
  Teuchos::ParameterList::ConstIterator i;
  for (i = fileParams.begin(); i !=fileParams.end(); ++i) {

    if (fileParams.isType<int>(fileParams.name(i)))
      dummyInt = fileParams.get<int>(fileParams.name(i));

    if (fileParams.isType<double>(fileParams.name(i)))
      dummyDouble = fileParams.get<double>(fileParams.name(i));
  }

  return true;
}
