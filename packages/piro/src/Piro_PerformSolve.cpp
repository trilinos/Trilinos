// @HEADER
// *****************************************************************************
//        Piro: Strategy package for embedded analysis capabilitites
//
// Copyright 2010 NTESS and the Piro contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Piro_PerformSolve.hpp"

#include "Teuchos_TestForException.hpp"

namespace Piro {

namespace Detail {

Teuchos::Array<bool> createResponseTable(
    int count,
    const std::string selectionType,
    int index,
    const Teuchos::ArrayView<const int> &list)
{
  Teuchos::Array<bool> result;

  if (count > 0) {
    if (selectionType == "All") {
      result.resize(count, true);
    } else if (selectionType == "Last") {
      result = createResponseTableFromIndex(count - 1, count);
    } else if (selectionType == "AllButLast") {
      result.reserve(count);
      result.resize(count - 1, true);
      result.push_back(false);
    } else if (selectionType == "Index") {
      result = createResponseTableFromIndex(index, count);
    } else if (selectionType == "List") {
      result.resize(count, false);
      for (Teuchos::ArrayView<const int>::const_iterator it = list.begin(), it_end = list.end(); it != it_end; ++it) {
        result.at(*it) = true;
      }
    } else {
      TEUCHOS_TEST_FOR_EXCEPT(false);
    }
  }

  return result;
}

Teuchos::Array<bool> parseResponseParameters(Teuchos::ParameterList &params, int responseCount)
{
  const std::string selectionType = params.get("Response Selection", "All");
  const int userIndex = parseResponseIndex(params);
  const Teuchos::Array<int> userList = params.get("Response List", Teuchos::Array<int>());
  return createResponseTable(responseCount, selectionType, userIndex, userList);
}

int parseResponseIndex(Teuchos::ParameterList &params)
{
  return params.get("Response Index", 0);
}

bool parseSensitivityParameters(Teuchos::ParameterList &params)
{
  return params.get("Compute Sensitivities", false);
}

Teuchos::Array<bool> createResponseTableFromIndex(int index, int responseCount)
{
  Teuchos::Array<bool> result(responseCount, false);
  result.at(index) = true;
  return result;
}

} // namespace Detail

} // namespace Piro
