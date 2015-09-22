// @HEADER
// ************************************************************************
//
//        Piro: Strategy package for embedded analysis capabilitites
//                  Copyright (2010) Sandia Corporation
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
// Questions? Contact Andy Salinger (agsalin@sandia.gov), Sandia
// National Laboratories.
//
// ************************************************************************
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
