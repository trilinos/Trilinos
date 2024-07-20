// @HEADER
// *****************************************************************************
//                    Teuchos: Common Tools Package
//
// Copyright 2004 NTESS and the Teuchos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
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
