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
#ifndef MUELU_STRUCTUREDLINEDETECTIONFACTORY_DEF_HPP
#define MUELU_STRUCTUREDLINEDETECTIONFACTORY_DEF_HPP

#include <Xpetra_Matrix.hpp>

#include "MueLu_StructuredLineDetectionFactory_decl.hpp"

#include "MueLu_Level.hpp"

namespace MueLu {

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
RCP<const ParameterList> StructuredLineDetectionFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::GetValidParameterList() const {
  RCP<ParameterList> validParamList = rcp(new ParameterList());

  validParamList->set<RCP<const FactoryBase> >("A", Teuchos::null, "Generating factory of the matrix A");
  validParamList->set<std::string>("orientation", "Z", "Lines orientation");
  validParamList->set<RCP<const FactoryBase> >("lNodesPerDim", Teuchos::null, "Number of nodes per spatial dimension provided by CoordinatesTransferFactory.");

  return validParamList;
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void StructuredLineDetectionFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::DeclareInput(Level& currentLevel) const {
  Input(currentLevel, "A");
  // Request the global number of nodes per dimensions
  if (currentLevel.GetLevelID() == 0) {
    if (currentLevel.IsAvailable("lNodesPerDim", NoFactory::get())) {
      currentLevel.DeclareInput("lNodesPerDim", NoFactory::get(), this);
    } else {
      TEUCHOS_TEST_FOR_EXCEPTION(currentLevel.IsAvailable("gNodesPerDim", NoFactory::get()),
                                 Exceptions::RuntimeError,
                                 "lNodesPerDim was not provided by the user on level0!");
    }
  } else {
    Input(currentLevel, "lNodesPerDim");
  }
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void StructuredLineDetectionFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Build(Level& currentLevel) const {
  // The following three variables are needed by the line smoothers in Ifpack/Ifpack2
  LO NumZDir                       = 0;
  Teuchos::ArrayRCP<LO> VertLineId = Teuchos::arcp<LO>(0);

  // collect information provided by user
  const ParameterList& pL           = GetParameterList();
  const std::string lineOrientation = pL.get<std::string>("orientation");

  // Extract data from currentLevel
  RCP<Matrix> A          = Get<RCP<Matrix> >(currentLevel, "A");
  Array<LO> lNodesPerDir = Get<Array<LO> >(currentLevel, "lNodesPerDim");
  LO numNodes            = lNodesPerDir[0] * lNodesPerDir[1] * lNodesPerDir[2];
  VertLineId.resize(numNodes);
  if (lineOrientation == "X") {
    NumZDir = lNodesPerDir[0];
  } else if (lineOrientation == "Y") {
    NumZDir = lNodesPerDir[1];
  } else if (lineOrientation == "Z") {
    NumZDir = lNodesPerDir[2];
  }

  for (LO k = 0; k < lNodesPerDir[2]; ++k) {
    for (LO j = 0; j < lNodesPerDir[1]; ++j) {
      for (LO i = 0; i < lNodesPerDir[0]; ++i) {
        if (lineOrientation == "X") {
          VertLineId[k * lNodesPerDir[1] * lNodesPerDir[0] + j * lNodesPerDir[0] + i] = k * lNodesPerDir[1] + j;
        } else if (lineOrientation == "Y") {
          VertLineId[k * lNodesPerDir[1] * lNodesPerDir[0] + j * lNodesPerDir[0] + i] = k * lNodesPerDir[0] + i;
        } else if (lineOrientation == "Z") {
          VertLineId[k * lNodesPerDir[1] * lNodesPerDir[0] + j * lNodesPerDir[0] + i] = j * lNodesPerDir[0] + i;
        }
      }
    }
  }

  Set(currentLevel, "CoarseNumZLayers", NumZDir);
  Set(currentLevel, "LineDetection_VertLineIds", VertLineId);
}

}  // namespace MueLu

#endif  // MUELU_STRUCTUREDLINEDETECTIONFACTORY_DEF_HPP
