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

#include <string>
#include <fstream>
#include <streambuf>

void Teuchos::updateParametersFromYamlFileAndBroadcast(
  const std::string &yamlFileName,
  const Ptr<ParameterList> &paramList,
  const Comm<int> &comm,
  bool overwrite
  )
{
  if (comm.getSize()==1)
    updateParametersFromYamlFile(yamlFileName, paramList);
  else {
    if (comm.getRank()==0) {
      std::ifstream stream(yamlFileName.c_str());
      TEUCHOS_TEST_FOR_EXCEPTION(!stream.is_open(),
          std::runtime_error,
          "Could not open YAML file " << yamlFileName);
      std::istreambuf_iterator<char> stream_iter(stream);
      std::istreambuf_iterator<char> stream_end;
      std::string yamlString(stream_iter, stream_end);
      int strsize = yamlString.size();
      broadcast<int, int>(comm, 0, &strsize);
      char* ptr = (strsize) ? (&yamlString[0]) : 0;
      broadcast<int, char>(comm, 0, strsize, ptr);
      updateParametersFromYamlString(yamlString, paramList,overwrite, yamlFileName);
    }
    else {
      int strsize;
      broadcast<int, int>(comm, 0, &strsize);
      std::string yamlString;
      yamlString.resize(strsize);
      char* ptr = (strsize) ? (&yamlString[0]) : 0;
      broadcast<int, char>(comm, 0, strsize, ptr);
      updateParametersFromYamlString(yamlString, paramList,overwrite);
    }
  }
}
