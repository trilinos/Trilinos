// @HEADER
// ***********************************************************************
//
//                    Teuchos: Common Tools Package
//                 Copyright (2004) Sandia Corporation
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
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
//
// ***********************************************************************
// @HEADER


#include "Teuchos_YamlParameterListHelpers.hpp"
#include "Teuchos_FileInputSource.hpp"
#include "Teuchos_CommHelpers.hpp"

void Teuchos::updateParametersFromYamlFileAndBroadcast(
  const std::string &yamlFileName, 
  const Teuchos::Ptr<Teuchos::ParameterList> &paramList, 
  const Teuchos::Comm<int> &comm, 
  bool overwrite)
{
  struct SafeFile
  {
    SafeFile(const char* fname, const char* options)
    {
      handle = fopen(fname, options);
    }
    ~SafeFile()
    {
      if(handle)
        fclose(handle);
    }
    FILE* handle;
  };
  //BMK note: see teuchos/comm/src/Teuchos_XMLParameterListHelpers.cpp
  if(comm.getSize() == 1)
  {
    updateParametersFromYamlFile(yamlFileName, paramList);
  }
  else
  {
    if(comm.getRank() == 0)
    {
      //BMK: TODO! //reader.setAllowsDuplicateSublists(false);
      //create a string and load file contents into it
      //C way for readability and speed, same thing with C++ streams is slow & ugly
      SafeFile yamlFile(yamlFileName.c_str(), "rb");
      if(!yamlFile.handle)
      {
        throw std::runtime_error(std::string("Failed to open YAML file \"") + yamlFileName + "\"for reading.");
      }
      fseek(yamlFile.handle, 0, SEEK_END);
      int strsize = ftell(yamlFile.handle) + 1;
      rewind(yamlFile.handle);
      //Make the array raii
      Teuchos::ArrayRCP<char> contents(new char[strsize], 0, strsize, true);
      fread((void*) contents.get(), strsize - 1, 1, yamlFile.handle);
      contents.get()[strsize - 1] = 0;
      Teuchos::broadcast<int, int>(comm, 0, &strsize);
      Teuchos::broadcast<int, char>(comm, 0, strsize, contents.get());
      updateParametersFromYamlCString(contents.get(), paramList, overwrite);
    }
    else
    {
      int strsize;
      Teuchos::broadcast<int, int>(comm, 0, &strsize);
      Teuchos::ArrayRCP<char> contents(new char[strsize], 0, strsize, true);
      Teuchos::broadcast<int, char>(comm, 0, strsize, contents.get());
      updateParametersFromYamlCString(contents.get(), paramList, overwrite);
    }
  }
}
