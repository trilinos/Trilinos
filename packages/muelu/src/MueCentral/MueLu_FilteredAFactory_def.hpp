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
//                    Jeremie Gaidamour (jngaida@sandia.gov)
//                    Jonathan Hu       (jhu@sandia.gov)
//                    Ray Tuminaro      (rstumin@sandia.gov)
//
// ***********************************************************************
//
// @HEADER
#ifndef MUELU_FILTEREDAFACTORY_DEF_HPP
#define MUELU_FILTEREDAFACTORY_DEF_HPP

#include <Xpetra_Matrix.hpp>
#include <Xpetra_MatrixFactory.hpp>

#include "MueLu_FilteredAFactory_decl.hpp"

#include "MueLu_Level.hpp"
#include "MueLu_Monitor.hpp"
#include "MueLu_FactoryManager.hpp"

namespace MueLu {

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  RCP<const ParameterList> FilteredAFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::GetValidParameterList(const ParameterList& paramList) const {
    RCP<ParameterList> validParamList = rcp(new ParameterList());

    validParamList->set< RCP<const FactoryBase> >("A",              Teuchos::null, "Generating factory of the matrix A used for filtering");
    validParamList->set< RCP<const FactoryBase> >("Graph",          Teuchos::null, "Generating fatory for coalesced filtered graph");
    validParamList->set< bool >                  ("lumping",                false, "Use lumping for dropped values");

    return validParamList;
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  void FilteredAFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::DeclareInput(Level& currentLevel) const {
    Input(currentLevel, "A");
    Input(currentLevel, "Graph");
    // NOTE: we do this DeclareInput in such complicated fashion because this is not a part of the parameter list
    currentLevel.DeclareInput("Filtering", currentLevel.GetFactoryManager()->GetFactory("Filtering").get());
  }

  // TODO: rewrite the function using AmalgamationInfo
  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  void FilteredAFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::Build(Level& currentLevel) const {
    using Teuchos::as;

    FactoryMonitor m(*this, "Matrix filtering", currentLevel);

    RCP<Matrix> A = Get< RCP<Matrix> >(currentLevel, "A");
    if (currentLevel.Get<bool>("Filtering", currentLevel.GetFactoryManager()->GetFactory("Filtering").get()) == false) {
      GetOStream(Runtime0,0) << "Filtered matrix is not being constructed as no filtering is being done" << std::endl;
      Set(currentLevel, "A", A);
      return;
    }

    const ParameterList& pL = GetParameterList();
    RCP<GraphBase>  G = Get< RCP<GraphBase> >(currentLevel, "Graph");
    bool      lumping = pL.get<bool>("lumping");
    size_t    blkSize = A->GetFixedBlockSize();

    if (lumping)
      GetOStream(Runtime0,0) << "Lumping dropped entries" << std::endl;

    ArrayView<const GO> GIDs = A->getColMap()->getNodeElementList();

    // NOTE: the good thing is that we mostly deal with local IDs

    // Calculate max entries per row
    RCP<Matrix> filteredA = MatrixFactory::Build(A->getRowMap(), A->getColMap(), A->getNodeMaxNumRowEntries(), Xpetra::StaticProfile);

    Array<GO>   newInds;
    Array<SC>   newVals;
    Array<char> filter(blkSize*G->GetImportMap()->getNodeNumElements(), 0);

    size_t numGRows = G->GetNodeNumVertices(), numInds = 0, diagIndex;
    SC diagExtra;
    for (size_t i = 0; i < numGRows; i++) {
      // Set up filtering array
      Teuchos::ArrayView<const LO> indsG = G->getNeighborVertices(i);
      for (size_t j = 0; j < as<size_t> (indsG.size()); j++)
        for (size_t k = 0; k < blkSize; k++)
          filter[indsG[j]*blkSize+k] = 1;

      for (size_t k = 0; k < blkSize; k++) {
        LocalOrdinal row = i*blkSize+k;
        ArrayView<const LO> oldInds;
        ArrayView<const SC> oldVals;
        A->getLocalRowView(row, oldInds, oldVals);

        diagIndex = (size_t)(-1);
        diagExtra = Teuchos::ScalarTraits<SC>::zero();

        newInds.resize(oldInds.size());
        newVals.resize(oldVals.size());
        numInds = 0;
        for (size_t j = 0; j < as<size_t> (oldInds.size()); j++)
          if (filter[oldInds[j]]) {
            newInds[numInds] = oldInds[j];
            newVals[numInds] = oldVals[j];

            // Remember diagonal position
            if (newInds[numInds] == row)
              diagIndex = numInds;
            numInds++;

          } else {
            diagExtra += oldVals[j];
          }
        // Lump dropped entries
        // NOTE
        //  * Does it make sense to lump for elasticity?
        //  * Is it different for diffusion and elasticity?
        if (lumping)
          newVals[diagIndex] += diagExtra;

        newInds.resize(numInds);
        newVals.resize(numInds);

        // NOTE: this is the only place where we do need GIDs
        for (size_t j = 0; j < numInds; j++)
          newInds[j] = GIDs[newInds[j]];

        filteredA->insertGlobalValues(GIDs[row], newInds, newVals);
      }

      // Clean up filtering array
      for (size_t j = 0; j < as<size_t> (indsG.size()); j++)
        for (size_t k = 0; k < blkSize; k++)
          filter[indsG[j]*blkSize+k] = 0;
    }
    RCP<ParameterList> fillCompleteParams(new ParameterList);;
    fillCompleteParams->set("No Nonlocal Changes", true);
    filteredA->fillComplete(A->getDomainMap(), A->getRangeMap(), fillCompleteParams);

    filteredA->SetFixedBlockSize(blkSize);

    // TODO: Can we reuse max eigenvalue from A?
    // filteredA->SetMaxEigenvalueEstimate(A->GetMaxEigenvalueEstimate());

    Set(currentLevel, "A", filteredA);
  }

} //namespace MueLu

#endif // MUELU_FILTEREDAFACTORY_DEF_HPP
