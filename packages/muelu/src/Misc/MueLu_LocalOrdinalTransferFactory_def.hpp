// @HEADER
// *****************************************************************************
//        MueLu: A package for multigrid based preconditioning
//
// Copyright 2012 NTESS and the MueLu contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef MUELU_LOCALORDINALTRANSFER_FACTORY_DEF_HPP
#define MUELU_LOCALORDINALTRANSFER_FACTORY_DEF_HPP

#include "Xpetra_ImportFactory.hpp"
#include "Xpetra_VectorFactory.hpp"
#include "Xpetra_MapFactory.hpp"
#include "Xpetra_CrsGraph.hpp"

#include "Xpetra_IO.hpp"

#include "MueLu_Aggregates.hpp"
#include "MueLu_LocalOrdinalTransferFactory_decl.hpp"

#include "MueLu_Level.hpp"
#include "MueLu_Monitor.hpp"

namespace MueLu {

template <class LocalOrdinal, class GlobalOrdinal, class Node>
RCP<const ParameterList> LocalOrdinalTransferFactory<LocalOrdinal, GlobalOrdinal, Node>::GetValidParameterList() const {
  RCP<ParameterList> validParamList = rcp(new ParameterList());

  validParamList->set<RCP<const FactoryBase> >(TransferVecName_, Teuchos::null, "Factory for TransferVec generation");
  validParamList->set<RCP<const FactoryBase> >("P Graph", Teuchos::null, "Factory for P generation");
  validParamList->set<RCP<const FactoryBase> >("Aggregates", Teuchos::null, "Factory for aggregates generation");
  validParamList->set<RCP<const FactoryBase> >("CoarseMap", Teuchos::null, "Generating factory of the coarse map");

  return validParamList;
}

template <class LocalOrdinal, class GlobalOrdinal, class Node>
void LocalOrdinalTransferFactory<LocalOrdinal, GlobalOrdinal, Node>::DeclareInput(Level &fineLevel, Level &coarseLevel) const {
  static bool isAvailableXfer = false;
  if (coarseLevel.GetRequestMode() == Level::REQUEST) {
    isAvailableXfer = coarseLevel.IsAvailable(TransferVecName_, this);
    if (isAvailableXfer == false) {
      Input(fineLevel, TransferVecName_);
      Input(fineLevel, "CoarseMap");

      if (useAggregatesMode_)
        Input(fineLevel, "Aggregates");
      else {
        Input(coarseLevel, "P Graph");
      }
    }
  }
}

template <class LocalOrdinal, class GlobalOrdinal, class Node>
void LocalOrdinalTransferFactory<LocalOrdinal, GlobalOrdinal, Node>::Build(Level &fineLevel, Level &coarseLevel) const {
  if (useAggregatesMode_)
    BuildAggregates(fineLevel, coarseLevel);
  else
    BuildFC(fineLevel, coarseLevel);
}

template <class LocalOrdinal, class GlobalOrdinal, class Node>
void LocalOrdinalTransferFactory<LocalOrdinal, GlobalOrdinal, Node>::BuildFC(Level &fineLevel, Level &coarseLevel) const {
  FactoryMonitor m(*this, "Build", coarseLevel);

  GetOStream(Runtime0) << "Transferring " << TransferVecName_ << std::endl;
  LO LO_INVALID = Teuchos::OrdinalTraits<LO>::invalid();

  if (coarseLevel.IsAvailable(TransferVecName_, this)) {
    GetOStream(Runtime0) << "Reusing " << TransferVecName_ << std::endl;
    return;
  }

  // Get everything we need
  RCP<const CrsGraph> P          = Get<RCP<const CrsGraph> >(coarseLevel, "P Graph");
  RCP<LocalOrdinalVector> fineTV = Get<RCP<LocalOrdinalVector> >(fineLevel, TransferVecName_);
  RCP<const Map> coarseMap       = Get<RCP<const Map> >(fineLevel, "CoarseMap");
  RCP<const Map> uniqueMap       = fineTV->getMap();
  ArrayRCP<const LO> fineData    = fineTV->getData(0);

  // Allocate new LO Vector
  RCP<LocalOrdinalVector> coarseTV = LocalOrdinalVectorFactory::Build(coarseMap, 1);
  ArrayRCP<LO> coarseData          = coarseTV->getDataNonConst(0);

  // Invalidate everything first, to check for errors
  for (LO i = 0; i < coarseData.size(); i++)
    coarseData[i] = LO_INVALID;

  // Fill in coarse TV
  LO domMapNumElements = P->getDomainMap()->getLocalNumElements();
  for (LO row = 0; row < (LO)P->getLocalNumRows(); row++) {
    LO fineNumber = fineData[row];
    ArrayView<const LO> indices;
    P->getLocalRowView(row, indices);

    for (LO j = 0; j < (LO)indices.size(); j++) {
      LO col = indices[j];
      if (col >= domMapNumElements) {
        // skip off rank entries of P
      } else {
        coarseData[col] = fineNumber;
      }
    }
  }

#ifdef HAVE_MUELU_DEBUG
  size_t error_count = 0;
  {
    RCP<LocalOrdinalVector> coarseTVghosted;
    RCP<const Import> importer = P->getImporter();
    if (!importer.is_null()) {
      coarseTVghosted = LocalOrdinalVectorFactory::Build(P->getColMap(), 1);
      coarseTVghosted->doImport(*coarseTV, *importer, Xpetra::INSERT);
    } else {
      coarseTVghosted = coarseTV;
    }
    ArrayRCP<LO> coarseDataGhosted = coarseTVghosted->getDataNonConst(0);
    for (LO col = 0; col < (LO)P->getColMap()->getLocalNumElements(); col++) {
      if (coarseDataGhosted[col] == LO_INVALID)
        error_count++;
    }
    for (LO row = 0; row < (LO)P->getLocalNumRows(); row++) {
      LO fineNumber = fineData[row];
      ArrayView<const LO> indices;
      P->getLocalRowView(row, indices);
      for (LO j = 0; j < (LO)indices.size(); j++) {
        if (coarseDataGhosted[indices[j]] != fineNumber)
          error_count++;
      }
    }
  }

  // Error checking:  All nodes in an aggregate must share a local ordinal
  if (error_count > 0) {
    std::ostringstream ofs;
    ofs << "LocalOrdinalTransferFactory(" << TransferVecName_ << "): ERROR:  Each coarse dof must have a unique LO value.  We had " << std::to_string(error_count) << " unknowns that did not match.";
    throw std::runtime_error(ofs.str());
  }
#endif

  Set<RCP<LocalOrdinalVector> >(coarseLevel, TransferVecName_, coarseTV);
}

template <class LocalOrdinal, class GlobalOrdinal, class Node>
void LocalOrdinalTransferFactory<LocalOrdinal, GlobalOrdinal, Node>::BuildAggregates(Level &fineLevel, Level &coarseLevel) const {
  FactoryMonitor m(*this, "Build", coarseLevel);

  GetOStream(Runtime0) << "Transferring " << TransferVecName_ << std::endl;
  RCP<LocalOrdinalVector> coarseTV;
  RCP<LocalOrdinalVector> fineTV;
  LO LO_INVALID = Teuchos::OrdinalTraits<LO>::invalid();

  if (coarseLevel.IsAvailable(TransferVecName_, this)) {
    GetOStream(Runtime0) << "Reusing " << TransferVecName_ << std::endl;
    return;
  }

  RCP<Aggregates> aggregates = Get<RCP<Aggregates> >(fineLevel, "Aggregates");
  fineTV                     = Get<RCP<LocalOrdinalVector> >(fineLevel, TransferVecName_);
  RCP<const Map> coarseMap   = Get<RCP<const Map> >(fineLevel, "CoarseMap");
  RCP<const Map> uniqueMap   = fineTV->getMap();

  ArrayView<const GO> elementAList = coarseMap->getLocalElementList();

  coarseTV = LocalOrdinalVectorFactory::Build(coarseMap, 1);

  // Create overlapped fine TV to reduce global communication
  RCP<LocalOrdinalVector> ghostedTV = fineTV;
  if (aggregates->AggregatesCrossProcessors()) {
    RCP<const Map> nonUniqueMap = aggregates->GetMap();
    RCP<const Import> importer  = ImportFactory::Build(uniqueMap, nonUniqueMap);

    ghostedTV = LocalOrdinalVectorFactory::Build(nonUniqueMap, 1);
    ghostedTV->doImport(*fineTV, *importer, Xpetra::INSERT);
  }

  // Get some info about aggregates
  int myPID                             = uniqueMap->getComm()->getRank();
  ArrayRCP<LO> aggSizes                 = aggregates->ComputeAggregateSizesArrayRCP();
  const ArrayRCP<const LO> vertex2AggID = aggregates->GetVertex2AggId()->getData(0);
  const ArrayRCP<const LO> procWinner   = aggregates->GetProcWinner()->getData(0);

  ArrayRCP<const LO> fineData = ghostedTV->getData(0);
  ArrayRCP<LO> coarseData     = coarseTV->getDataNonConst(0);

  // Invalidate everything first, to check for errors
  for (LO i = 0; i < coarseData.size(); i++)
    coarseData[i] = LO_INVALID;

  // Fill in coarse TV
  size_t error_count = 0;
  for (LO lnode = 0; lnode < vertex2AggID.size(); lnode++) {
    if (procWinner[lnode] == myPID &&
        // lnode < vertex2AggID.size() &&
        lnode < fineData.size() &&  // TAW do not access off-processor data
        vertex2AggID[lnode] < coarseData.size()) {
      if (coarseData[vertex2AggID[lnode]] == LO_INVALID)
        coarseData[vertex2AggID[lnode]] = fineData[lnode];
      if (coarseData[vertex2AggID[lnode]] != fineData[lnode])
        error_count++;
    }
  }

  // Error checking:  All nodes in an aggregate must share a local ordinal
  if (error_count > 0) {
    std::ostringstream ofs;
    ofs << "LocalOrdinalTransferFactory: ERROR:  Each aggregate must have a unique LO value.  We had " << std::to_string(error_count) << " unknowns that did not match.";
    throw std::runtime_error(ofs.str());
  }

  Set<RCP<LocalOrdinalVector> >(coarseLevel, TransferVecName_, coarseTV);
}

}  // namespace MueLu

#endif  // MUELU_LOCALORDINALTRANSFER_FACTORY_DEF_HPP
