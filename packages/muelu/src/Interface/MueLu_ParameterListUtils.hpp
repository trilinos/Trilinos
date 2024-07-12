// @HEADER
// *****************************************************************************
//        MueLu: A package for multigrid based preconditioning
//
// Copyright 2012 NTESS and the MueLu contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
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
