// @HEADER
//
// ***********************************************************************
//
//        MueLu: A package for multigrid based preconditioning
//                  Copyright 2012 Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
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
// Questions? Contact
//                    Jonathan Hu       (jhu@sandia.gov)
//                    Andrey Prokopenko (aprokop@sandia.gov)
//                    Ray Tuminaro      (rstumin@sandia.gov)
//
// ***********************************************************************
//
// @HEADER

#ifndef MUELU_PARAMETERLISTUTILS_HPP
#define MUELU_PARAMETERLISTUTILS_HPP

#include <string>
#include <sstream>

#include <Teuchos_ParameterList.hpp>
#include <Teuchos_RCPDecl.hpp>
#include "MueLu_ConfigDefs.hpp"

namespace MueLu {

/* See also: ML_Epetra::UpdateList */
void MergeParameterList(const Teuchos::ParameterList& source, Teuchos::ParameterList& dest, bool overWrite);

void CreateSublists(const Teuchos::ParameterList& List, Teuchos::ParameterList& newList);

// Usage: GetMLSubList(paramList, "smoother", 2);
const Teuchos::ParameterList& GetMLSubList(const Teuchos::ParameterList& paramList, const std::string& type, int levelID);

// Extract all the parameters that begin with "str:" (but skip sublist)
Teuchos::RCP<Teuchos::ParameterList> ExtractSetOfParameters(const Teuchos::ParameterList& paramList, const std::string& str);

//! replace all string occurrences "from" with "to" in "str"
//!
//! @param str: input and output string
//! @param from: search string
//! @param to: replace with "to"
void replaceAll(std::string& str, const std::string& from, const std::string& to);

//! templated version to replace placeholder by data in "str"
template <typename Type>
bool replacePlaceholder(std::string& str, const std::string& placeholder, Type data) {
  std::stringstream s;
  s << data;
  replaceAll(str, placeholder, s.str());
  return true;
}

template <typename Type>
bool actionInterpretParameter(Teuchos::ParameterList& mlParams, const std::string& paramName, std::string& str) {
  // MUELU_READ_PARAM(mlParams, paramName, int, 0, data);

  Type varName;  // = defaultValue; // extract from master list
  if (mlParams.isParameter(paramName)) varName = mlParams.get<Type>(paramName);

  std::stringstream placeholder;
  placeholder << "$" << paramName << "$";

  return MueLu::replacePlaceholder<Type>(str, placeholder.str(), varName);
}

}  // namespace MueLu

#endif  // MUELU_PARAMETERLISTUTILS_HPP
