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
#ifndef PACKAGES_MUELU_SRC_REBALANCING_MUELU_CLONEREPARTITIONINTERFACE_DEF_HPP_
#define PACKAGES_MUELU_SRC_REBALANCING_MUELU_CLONEREPARTITIONINTERFACE_DEF_HPP_

#include <Xpetra_MultiVectorFactory.hpp>
#include <Xpetra_VectorFactory.hpp>

#include "MueLu_CloneRepartitionInterface_decl.hpp"

#include "MueLu_Level.hpp"
#include "MueLu_Exceptions.hpp"
#include "MueLu_Monitor.hpp"

namespace MueLu {

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
Teuchos::RCP<const ParameterList> CloneRepartitionInterface<Scalar, LocalOrdinal, GlobalOrdinal, Node>::GetValidParameterList() const {
  Teuchos::RCP<ParameterList> validParamList = rcp(new ParameterList());
  validParamList->set<Teuchos::RCP<const FactoryBase> >("A", Teuchos::null, "Factory of the matrix A");
  validParamList->set<Teuchos::RCP<const FactoryBase> >("number of partitions", Teuchos::null, "Instance of RepartitionHeuristicFactory.");
  validParamList->set<Teuchos::RCP<const FactoryBase> >("Partition", Teuchos::null, "Factory generating the Partition array (e.g. ZoltanInterface)");

  return validParamList;
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void CloneRepartitionInterface<Scalar, LocalOrdinal, GlobalOrdinal, Node>::DeclareInput(Level &currentLevel) const {
  Input(currentLevel, "A");
  Input(currentLevel, "Partition");
}  // DeclareInput()

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void CloneRepartitionInterface<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Build(Level &currentLevel) const {
  FactoryMonitor m(*this, "Build", currentLevel);
  currentLevel.print(GetOStream(Statistics0, 0));
  // extract blocked operator A from current level
  Teuchos::RCP<Matrix> A                       = Get<Teuchos::RCP<Matrix> >(currentLevel, "A");
  Teuchos::RCP<const Teuchos::Comm<int> > comm = A->getRowMap()->getComm();

  // number of Partitions only used for a shortcut.
  if (currentLevel.IsAvailable("number of partitions")) {
    GetOStream(Warnings0) << "Using user-provided \"number of partitions\", the performance is unknown" << std::endl;
  }

  // ======================================================================================================
  // Construct decomposition vector
  // ======================================================================================================
  RCP<GOVector> decomposition = Teuchos::null;

  // extract decomposition vector
  decomposition                    = Get<RCP<GOVector> >(currentLevel, "Partition");
  ArrayRCP<const GO> decompEntries = decomposition->getData(0);

  if (decomposition.is_null()) {
    GetOStream(Warnings0) << "No repartitioning necessary: partitions were left unchanged by the repartitioner" << std::endl;
    Set<RCP<const Import> >(currentLevel, "Importer", Teuchos::null);
    return;
  }

  // create new decomposition vector
  Teuchos::RCP<GOVector> ret    = Xpetra::VectorFactory<GO, LO, GO, NO>::Build(A->getRowMap(), false);
  ArrayRCP<GO> retDecompEntries = ret->getDataNonConst(0);

  // block size of output vector
  LocalOrdinal blkSize = 1;

  // check for blocking/striding information
  if (A->IsView("stridedMaps") &&
      Teuchos::rcp_dynamic_cast<const StridedMap>(A->getRowMap("stridedMaps")) != Teuchos::null) {
    Xpetra::viewLabel_t oldView  = A->SwitchToView("stridedMaps");  // note: "stridedMaps are always non-overlapping (correspond to range and domain maps!)
    RCP<const StridedMap> strMap = Teuchos::rcp_dynamic_cast<const StridedMap>(A->getRowMap());
    TEUCHOS_TEST_FOR_EXCEPTION(strMap == Teuchos::null, Exceptions::BadCast, "MueLu::CloneRepartitionInterface::Build: cast to strided row map failed.");
    LocalOrdinal stridedBlock = strMap->getStridedBlockId();
    if (stridedBlock == -1)
      blkSize = strMap->getFixedBlockSize();
    else {
      std::vector<size_t> strInfo = strMap->getStridingData();
      blkSize                     = strInfo[stridedBlock];
    }
    oldView = A->SwitchToView(oldView);
    GetOStream(Statistics1) << "CloneRepartitionInterface::Build():"
                            << " found blockdim=" << blkSize << " from strided maps." << std::endl;
  } else {
    GetOStream(Statistics1) << "CloneRepartitionInterface::Build(): no striding information available. Use blockdim=" << blkSize << " (DofsPerNode)." << std::endl;
    blkSize = A->GetFixedBlockSize();
  }

  // plausibility check!
  size_t inLocalLength  = decomposition->getLocalLength();
  size_t outLocalLength = A->getRowMap()->getLocalNumElements();

  // only for non-strided maps
  size_t numLocalNodes = outLocalLength / blkSize;
  TEUCHOS_TEST_FOR_EXCEPTION(Teuchos::as<size_t>(outLocalLength % blkSize) != 0, MueLu::Exceptions::RuntimeError, "CloneRepartitionInterface: inconsistent number of local DOFs (" << outLocalLength << ") and degrees of freedoms (" << blkSize << ")");

  if (numLocalNodes > 0) {
    TEUCHOS_TEST_FOR_EXCEPTION(Teuchos::as<size_t>(inLocalLength % numLocalNodes) != 0, MueLu::Exceptions::RuntimeError, "CloneRepartitionInterface: inconsistent number of local DOFs (" << inLocalLength << ") and number of local nodes (" << numLocalNodes << ")");
    LocalOrdinal inBlkSize = Teuchos::as<LocalOrdinal>(inLocalLength / numLocalNodes);
    // TEUCHOS_TEST_FOR_EXCEPTION(blkSize != inBlkSize, MueLu::Exceptions::RuntimeError,"CloneRepartitionInterface: input block size = " << inBlkSize << " outpub block size = " << blkSize << ". They should be the same.");

    for (LO i = 0; i < Teuchos::as<LO>(numLocalNodes); i++) {
      for (LO j = 0; j < blkSize; j++) {
        retDecompEntries[i * blkSize + j] = Teuchos::as<GO>(decompEntries[i * inBlkSize]);
      }
    }
  }  // end if numLocalNodes > 0
  Set(currentLevel, "Partition", ret);
}  // Build()

}  // namespace MueLu

#endif /* PACKAGES_MUELU_SRC_REBALANCING_MUELU_CLONEREPARTITIONINTERFACE_DEF_HPP_ */
