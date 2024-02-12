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
#ifndef MUELU_REPLICATEPFACTORY_DEF_HPP
#define MUELU_REPLICATEPFACTORY_DEF_HPP

#include <stdlib.h>
#include <iomanip>

// #include <Teuchos_LAPACK.hpp>
#include <Teuchos_SerialDenseMatrix.hpp>
#include <Teuchos_SerialDenseVector.hpp>
#include <Teuchos_SerialDenseSolver.hpp>

#include <Xpetra_CrsMatrixWrap.hpp>
#include <Xpetra_ImportFactory.hpp>
#include <Xpetra_Matrix.hpp>
#include <Xpetra_MapFactory.hpp>
#include <Xpetra_MultiVectorFactory.hpp>
#include <Xpetra_VectorFactory.hpp>

#include <Xpetra_IO.hpp>

#include "MueLu_ReplicatePFactory_decl.hpp"

#include "MueLu_Monitor.hpp"

namespace MueLu {

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
RCP<const ParameterList> ReplicatePFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::GetValidParameterList() const {
  RCP<ParameterList> validParamList = rcp(new ParameterList());
  validParamList->setEntry("replicate: npdes", ParameterEntry(1));

  return validParamList;
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void ReplicatePFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    DeclareInput(Level& fineLevel, Level& /* coarseLevel */) const {
  //    Input(fineLevel, "Psubblock");
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void ReplicatePFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Build(Level& fineLevel,
                                                                         Level& coarseLevel) const {
  return BuildP(fineLevel, coarseLevel);
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void ReplicatePFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::BuildP(Level& fineLevel,
                                                                          Level& coarseLevel) const {
  FactoryMonitor m(*this, "Build", coarseLevel);

  RCP<Matrix> Psubblock   = coarseLevel.Get<RCP<Matrix> >("Psubblock", NoFactory::get());
  const ParameterList& pL = GetParameterList();
  const LO dofPerNode     = as<LO>(pL.get<int>("replicate: npdes"));

  Teuchos::ArrayRCP<const size_t> amalgRowPtr(Psubblock->getLocalNumRows());
  Teuchos::ArrayRCP<const LocalOrdinal> amalgCols(Psubblock->getLocalNumEntries());
  Teuchos::ArrayRCP<const Scalar> amalgVals(Psubblock->getLocalNumEntries());
  Teuchos::RCP<CrsMatrixWrap> Psubblockwrap = Teuchos::rcp_dynamic_cast<CrsMatrixWrap>(Psubblock);
  Teuchos::RCP<CrsMatrix> Psubblockcrs      = Psubblockwrap->getCrsMatrix();
  Psubblockcrs->getAllValues(amalgRowPtr, amalgCols, amalgVals);

  size_t paddedNrows = Psubblock->getRowMap()->getLocalNumElements() * Teuchos::as<size_t>(dofPerNode);
  Teuchos::ArrayRCP<size_t> newPRowPtr(paddedNrows + 1);
  Teuchos::ArrayRCP<LocalOrdinal> newPCols(Psubblock->getLocalNumEntries() * dofPerNode);
  Teuchos::ArrayRCP<Scalar> newPVals(Psubblock->getLocalNumEntries() * dofPerNode);
  size_t cnt = 0;  // local id counter
  for (decltype(amalgRowPtr.size()) i = 0; i < amalgRowPtr.size() - 1; i++) {
    size_t rowLength = amalgRowPtr[i + 1] - amalgRowPtr[i];
    for (int j = 0; j < dofPerNode; j++) {
      newPRowPtr[i * dofPerNode + j] = cnt;
      for (size_t k = 0; k < rowLength; k++) {
        newPCols[cnt]   = amalgCols[k + amalgRowPtr[i]] * dofPerNode + j;
        newPVals[cnt++] = amalgVals[k + amalgRowPtr[i]];
      }
    }
  }

  newPRowPtr[paddedNrows] = cnt;  // close row CSR array
  std::vector<size_t> stridingInfo(1, dofPerNode);

  GlobalOrdinal nCoarseDofs      = Psubblock->getDomainMap()->getLocalNumElements() * dofPerNode;
  GlobalOrdinal indexBase        = Psubblock->getDomainMap()->getIndexBase();
  RCP<const Map> coarseDomainMap = StridedMapFactory::Build(Psubblock->getDomainMap()->lib(),
                                                            Teuchos::OrdinalTraits<Xpetra::global_size_t>::invalid(),
                                                            nCoarseDofs,
                                                            indexBase,
                                                            stridingInfo,
                                                            Psubblock->getDomainMap()->getComm(),
                                                            -1 /* stridedBlockId */,
                                                            0 /*domainGidOffset */);

  size_t nColCoarseDofs = Teuchos::as<size_t>(Psubblock->getColMap()->getLocalNumElements() * dofPerNode);
  Teuchos::Array<GlobalOrdinal> unsmooshColMapGIDs(nColCoarseDofs);
  for (size_t c = 0; c < Psubblock->getColMap()->getLocalNumElements(); ++c) {
    GlobalOrdinal gid = (Psubblock->getColMap()->getGlobalElement(c) - indexBase) * dofPerNode + indexBase;

    for (int i = 0; i < dofPerNode; ++i) {
      unsmooshColMapGIDs[c * dofPerNode + i] = gid + i;
    }
  }
  Teuchos::RCP<Map> coarseColMap = MapFactory::Build(Psubblock->getDomainMap()->lib(),
                                                     Teuchos::OrdinalTraits<Xpetra::global_size_t>::invalid(),
                                                     unsmooshColMapGIDs(),  // View,
                                                     indexBase,
                                                     Psubblock->getDomainMap()->getComm());

  Teuchos::Array<GlobalOrdinal> unsmooshRowMapGIDs(paddedNrows);
  for (size_t c = 0; c < Psubblock->getRowMap()->getLocalNumElements(); ++c) {
    GlobalOrdinal gid = (Psubblock->getRowMap()->getGlobalElement(c) - indexBase) * dofPerNode + indexBase;

    for (int i = 0; i < dofPerNode; ++i) {
      unsmooshRowMapGIDs[c * dofPerNode + i] = gid + i;
    }
  }
  Teuchos::RCP<Map> fineRowMap = MapFactory::Build(Psubblock->getDomainMap()->lib(),
                                                   Teuchos::OrdinalTraits<Xpetra::global_size_t>::invalid(),
                                                   unsmooshRowMapGIDs(),  // View,
                                                   indexBase,
                                                   Psubblock->getDomainMap()->getComm());

  Teuchos::RCP<CrsMatrix> bigPCrs = CrsMatrixFactory::Build(fineRowMap, coarseColMap,
                                                            dofPerNode * Psubblock->getLocalMaxNumRowEntries());
  for (size_t i = 0; i < paddedNrows; i++) {
    bigPCrs->insertLocalValues(i,
                               newPCols.view(newPRowPtr[i], newPRowPtr[i + 1] - newPRowPtr[i]),
                               newPVals.view(newPRowPtr[i], newPRowPtr[i + 1] - newPRowPtr[i]));
  }
  bigPCrs->fillComplete(coarseDomainMap, fineRowMap);

  Teuchos::RCP<Matrix> bigP = Teuchos::rcp(new CrsMatrixWrap(bigPCrs));

  Set(coarseLevel, "P", bigP);
}

}  // namespace MueLu

#define MUELU_REPLICATEPFACTORY_SHORT
#endif  // MUELU_REPLICATEPFACTORY_DEF_HPP
