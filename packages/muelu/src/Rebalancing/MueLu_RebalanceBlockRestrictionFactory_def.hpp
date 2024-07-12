// @HEADER
// *****************************************************************************
//        MueLu: A package for multigrid based preconditioning
//
// Copyright 2012 NTESS and the MueLu contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef MUELU_REBALANCEBLOCKRESTRICTIONFACTORY_DEF_HPP_
#define MUELU_REBALANCEBLOCKRESTRICTIONFACTORY_DEF_HPP_

#include <Teuchos_Tuple.hpp>

#include "Xpetra_MultiVector.hpp"
#include "Xpetra_MultiVectorFactory.hpp"
#include "Xpetra_Vector.hpp"
#include "Xpetra_VectorFactory.hpp"
#include <Xpetra_Matrix.hpp>
#include <Xpetra_BlockedCrsMatrix.hpp>
#include <Xpetra_MapFactory.hpp>
#include <Xpetra_MapExtractor.hpp>
#include <Xpetra_MapExtractorFactory.hpp>
#include <Xpetra_MatrixFactory.hpp>
#include <Xpetra_Import.hpp>
#include <Xpetra_ImportFactory.hpp>

#include "MueLu_RebalanceBlockRestrictionFactory_decl.hpp"

#include "MueLu_HierarchyUtils.hpp"
#include "MueLu_FactoryManagerBase.hpp"
#include "MueLu_Level.hpp"
#include "MueLu_MasterList.hpp"
#include "MueLu_Monitor.hpp"
#include "MueLu_PerfUtils.hpp"

namespace MueLu {

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
RCP<const ParameterList> RebalanceBlockRestrictionFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::GetValidParameterList() const {
  RCP<ParameterList> validParamList = rcp(new ParameterList());

#define SET_VALID_ENTRY(name) validParamList->setEntry(name, MasterList::getEntry(name))
  SET_VALID_ENTRY("repartition: use subcommunicators");
#undef SET_VALID_ENTRY

  // validParamList->set< RCP<const FactoryBase> >("repartition: use subcommunicators", Teuchos::null, "test");

  validParamList->set<RCP<const FactoryBase> >("R", Teuchos::null, "Factory of the restriction operator that need to be rebalanced (only used if type=Restriction)");
  validParamList->set<RCP<const FactoryBase> >("Importer", Teuchos::null, "Generating factory of the matrix Importer for rebalancing");
  validParamList->set<RCP<const FactoryBase> >("SubImporters", Teuchos::null, "Generating factory of the matrix sub-Importers for rebalancing");
  validParamList->set<RCP<const FactoryBase> >("Nullspace", Teuchos::null, "Generating factory of the Nullspace operator");

  return validParamList;
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void RebalanceBlockRestrictionFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::AddFactoryManager(RCP<const FactoryManagerBase> FactManager) {
  FactManager_.push_back(FactManager);
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void RebalanceBlockRestrictionFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::DeclareInput(Level &fineLevel, Level &coarseLevel) const {
  Input(coarseLevel, "R");

  std::vector<Teuchos::RCP<const FactoryManagerBase> >::const_iterator it;
  for (it = FactManager_.begin(); it != FactManager_.end(); ++it) {
    SetFactoryManager fineSFM(rcpFromRef(fineLevel), *it);
    SetFactoryManager coarseSFM(rcpFromRef(coarseLevel), *it);

    if (!UseSingleSourceImporters_) coarseLevel.DeclareInput("Importer", (*it)->GetFactory("Importer").get(), this);
    coarseLevel.DeclareInput("Nullspace", (*it)->GetFactory("Nullspace").get(), this);
  }

  // Use the non-manager path if the maps / importers are generated in one place
  if (UseSingleSourceImporters_) {
    Input(coarseLevel, "SubImporters");
  }
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void RebalanceBlockRestrictionFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Build(Level &fineLevel, Level &coarseLevel) const {
  FactoryMonitor m(*this, "Build", coarseLevel);

  RCP<Teuchos::FancyOStream> out = Teuchos::fancyOStream(Teuchos::rcpFromRef(std::cout));

  Teuchos::RCP<Matrix> originalTransferOp = Teuchos::null;
  originalTransferOp                      = Get<RCP<Matrix> >(coarseLevel, "R");

  RCP<Xpetra::BlockedCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> > bOriginalTransferOp =
      Teuchos::rcp_dynamic_cast<Xpetra::BlockedCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> >(originalTransferOp);
  TEUCHOS_TEST_FOR_EXCEPTION(bOriginalTransferOp == Teuchos::null, Exceptions::BadCast, "MueLu::RebalanceBlockTransferFactory::Build: input matrix P or R is not of type BlockedCrsMatrix! error.");

  RCP<const MapExtractor> rangeMapExtractor  = bOriginalTransferOp->getRangeMapExtractor();
  RCP<const MapExtractor> domainMapExtractor = bOriginalTransferOp->getDomainMapExtractor();

  // restrict communicator?
  bool bRestrictComm      = false;
  const ParameterList &pL = GetParameterList();
  if (pL.get<bool>("repartition: use subcommunicators") == true)
    bRestrictComm = true;

  // check if GIDs for full maps have to be sorted:
  // For the Thyra mode ordering they do not have to be sorted since the GIDs are
  // numbered as 0...n1,0...,n2 (starting with zero for each subblock). The MapExtractor
  // generates unique GIDs during the construction.
  // For Xpetra style, the GIDs have to be reordered. Such that one obtains a ordered
  // list of GIDs in an increasing ordering. In Xpetra, the GIDs are all unique through
  // out all submaps.
  bool bThyraRangeGIDs  = rangeMapExtractor->getThyraMode();
  bool bThyraDomainGIDs = domainMapExtractor->getThyraMode();

  // rebuild rebalanced blocked P operator
  std::vector<GO> fullRangeMapVector;
  std::vector<GO> fullDomainMapVector;
  std::vector<RCP<const Map> > subBlockRRangeMaps;
  std::vector<RCP<const Map> > subBlockRDomainMaps;
  subBlockRRangeMaps.reserve(bOriginalTransferOp->Rows());   // reserve size for block P operators
  subBlockRDomainMaps.reserve(bOriginalTransferOp->Cols());  // reserve size for block P operators

  std::vector<Teuchos::RCP<Matrix> > subBlockRebR;
  subBlockRebR.reserve(bOriginalTransferOp->Cols());

  std::vector<RCP<const Import> > importers = std::vector<RCP<const Import> >(bOriginalTransferOp->Rows(), Teuchos::null);
  if (UseSingleSourceImporters_) {
    importers = Get<std::vector<RCP<const Import> > >(coarseLevel, "SubImporters");
  }

  int curBlockId = 0;
  Teuchos::RCP<const Import> rebalanceImporter;
  std::vector<Teuchos::RCP<const FactoryManagerBase> >::const_iterator it;
  for (it = FactManager_.begin(); it != FactManager_.end(); ++it) {
    // begin SubFactoryManager environment
    SetFactoryManager fineSFM(rcpFromRef(fineLevel), *it);
    SetFactoryManager coarseSFM(rcpFromRef(coarseLevel), *it);

    if (UseSingleSourceImporters_)
      rebalanceImporter = importers[curBlockId];
    else
      rebalanceImporter = coarseLevel.Get<Teuchos::RCP<const Import> >("Importer", (*it)->GetFactory("Importer").get());

    // extract matrix block
    Teuchos::RCP<Matrix> Rii = bOriginalTransferOp->getMatrix(curBlockId, curBlockId);

    // TODO run this only in the debug version
    TEUCHOS_TEST_FOR_EXCEPTION(bThyraRangeGIDs == true && Rii->getRowMap()->getMinAllGlobalIndex() != 0,
                               Exceptions::RuntimeError,
                               "MueLu::RebalanceBlockRestrictionFactory::Build: inconsistent Thyra GIDs. Thyra global ids for block range " << curBlockId << " start with " << Rii->getRowMap()->getMinAllGlobalIndex() << " but should start with 0");
    TEUCHOS_TEST_FOR_EXCEPTION(bThyraDomainGIDs == true && Rii->getColMap()->getMinAllGlobalIndex() != 0,
                               Exceptions::RuntimeError,
                               "MueLu::RebalanceBlockRestrictionFactory::Build: inconsistent Thyra GIDs. Thyra global ids for block domain " << curBlockId << " start with " << Rii->getColMap()->getMinAllGlobalIndex() << " but should start with 0");

    Teuchos::RCP<Matrix> rebRii;
    if (rebalanceImporter != Teuchos::null) {
      std::stringstream ss;
      ss << "Rebalancing restriction block R(" << curBlockId << "," << curBlockId << ")";
      SubFactoryMonitor m1(*this, ss.str(), coarseLevel);
      {
        SubFactoryMonitor subM(*this, "Rebalancing restriction -- fusedImport", coarseLevel);
        // Note: The 3rd argument says to use originalR's domain map.

        RCP<Map> dummy;
        rebRii = MatrixFactory::Build(Rii, *rebalanceImporter, dummy, rebalanceImporter->getTargetMap());
      }

      RCP<ParameterList> params = rcp(new ParameterList());
      params->set("printLoadBalancingInfo", true);
      std::stringstream ss2;
      ss2 << "R(" << curBlockId << "," << curBlockId << ") rebalanced:";
      GetOStream(Statistics0) << PerfUtils::PrintMatrixInfo(*rebRii, ss2.str(), params);
    } else {
      rebRii                    = Rii;
      RCP<ParameterList> params = rcp(new ParameterList());
      params->set("printLoadBalancingInfo", true);
      std::stringstream ss2;
      ss2 << "R(" << curBlockId << "," << curBlockId << ") not rebalanced:";
      GetOStream(Statistics0) << PerfUtils::PrintMatrixInfo(*rebRii, ss2.str(), params);
    }

    // fix striding information for rebalanced diagonal block rebRii
    Teuchos::RCP<const StridedMap> orig_stridedRgMap = Teuchos::rcp_dynamic_cast<const StridedMap>(rangeMapExtractor->getMap(Teuchos::as<size_t>(curBlockId), rangeMapExtractor->getThyraMode()));
    Teuchos::RCP<const Map> stridedRgMap             = Teuchos::null;
    if (orig_stridedRgMap != Teuchos::null) {
      std::vector<size_t> stridingData                       = orig_stridedRgMap->getStridingData();
      Teuchos::ArrayView<const GlobalOrdinal> nodeRangeMapii = rebRii->getRangeMap()->getLocalElementList();
      stridedRgMap                                           = StridedMapFactory::Build(
                                                    originalTransferOp->getRangeMap()->lib(),
                                                    Teuchos::OrdinalTraits<Xpetra::global_size_t>::invalid(),
                                                    nodeRangeMapii,
                                                    rebRii->getRangeMap()->getIndexBase(),
                                                    stridingData,
                                                    originalTransferOp->getRangeMap()->getComm(),
                                                    orig_stridedRgMap->getStridedBlockId(),
                                                    orig_stridedRgMap->getOffset());
    } else
      stridedRgMap = Rii->getRangeMap();

    Teuchos::RCP<const StridedMap> orig_stridedDoMap = Teuchos::rcp_dynamic_cast<const StridedMap>(domainMapExtractor->getMap(Teuchos::as<size_t>(curBlockId), domainMapExtractor->getThyraMode()));
    Teuchos::RCP<const Map> stridedDoMap             = Teuchos::null;
    if (orig_stridedDoMap != Teuchos::null) {
      std::vector<size_t> stridingData                        = orig_stridedDoMap->getStridingData();
      Teuchos::ArrayView<const GlobalOrdinal> nodeDomainMapii = rebRii->getDomainMap()->getLocalElementList();
      stridedDoMap                                            = StridedMapFactory::Build(
                                                     originalTransferOp->getDomainMap()->lib(),
                                                     Teuchos::OrdinalTraits<Xpetra::global_size_t>::invalid(),
                                                     nodeDomainMapii,
                                                     rebRii->getDomainMap()->getIndexBase(),
                                                     stridingData,
                                                     originalTransferOp->getDomainMap()->getComm(),
                                                     orig_stridedDoMap->getStridedBlockId(),
                                                     orig_stridedDoMap->getOffset());
    } else
      stridedDoMap = Rii->getDomainMap();

    if (bRestrictComm) {
      stridedRgMap->removeEmptyProcesses();
      stridedDoMap->removeEmptyProcesses();
    }

    TEUCHOS_TEST_FOR_EXCEPTION(stridedRgMap == Teuchos::null, Exceptions::RuntimeError, "MueLu::RebalanceBlockRestrictionFactory::Build: failed to generate striding information. error.");
    TEUCHOS_TEST_FOR_EXCEPTION(stridedDoMap == Teuchos::null, Exceptions::RuntimeError, "MueLu::RebalanceBlockRestrictionFactory::Build: failed to generate striding information. error.");

    // replace stridedMaps view in diagonal sub block
    if (rebRii->IsView("stridedMaps")) rebRii->RemoveView("stridedMaps");
    rebRii->CreateView("stridedMaps", stridedRgMap, stridedDoMap);

    // store rebalanced subblock
    subBlockRebR.push_back(rebRii);

    // append strided row map (= range map) to list of range maps.
    Teuchos::RCP<const Map> rangeMapii = rebRii->getRowMap("stridedMaps");
    subBlockRRangeMaps.push_back(rangeMapii);
    Teuchos::ArrayView<const GlobalOrdinal> nodeRangeMapii = rebRii->getRangeMap()->getLocalElementList();
    // append the GIDs in the end. Do not sort if we have Thyra style GIDs
    fullRangeMapVector.insert(fullRangeMapVector.end(), nodeRangeMapii.begin(), nodeRangeMapii.end());
    if (bThyraRangeGIDs == false)
      sort(fullRangeMapVector.begin(), fullRangeMapVector.end());

    // append strided col map (= domain map) to list of range maps.
    Teuchos::RCP<const Map> domainMapii = rebRii->getColMap("stridedMaps");
    subBlockRDomainMaps.push_back(domainMapii);
    Teuchos::ArrayView<const GlobalOrdinal> nodeDomainMapii = rebRii->getDomainMap()->getLocalElementList();
    // append the GIDs in the end. Do not sort if we have Thyra style GIDs
    fullDomainMapVector.insert(fullDomainMapVector.end(), nodeDomainMapii.begin(), nodeDomainMapii.end());
    if (bThyraDomainGIDs == false)
      sort(fullDomainMapVector.begin(), fullDomainMapVector.end());

    ////////////////////////////////////////////////////////////

    // rebalance null space
    // This rebalances the null space partial vector associated with the current block (generated by the NullspaceFactory
    // associated with the block)
    if (rebalanceImporter != Teuchos::null) {  // rebalance null space
      std::stringstream ss2;
      ss2 << "Rebalancing nullspace block(" << curBlockId << "," << curBlockId << ")";
      SubFactoryMonitor subM(*this, ss2.str(), coarseLevel);

      RCP<MultiVector> nullspace         = coarseLevel.Get<RCP<MultiVector> >("Nullspace", (*it)->GetFactory("Nullspace").get());
      RCP<MultiVector> permutedNullspace = MultiVectorFactory::Build(rebalanceImporter->getTargetMap(), nullspace->getNumVectors());
      permutedNullspace->doImport(*nullspace, *rebalanceImporter, Xpetra::INSERT);

      // TODO subcomm enabled everywhere or nowhere
      if (bRestrictComm)
        permutedNullspace->replaceMap(permutedNullspace->getMap()->removeEmptyProcesses());

      coarseLevel.Set<RCP<MultiVector> >("Nullspace", permutedNullspace, (*it)->GetFactory("Nullspace").get());

    }       // end rebalance null space
    else {  // do nothing
      RCP<MultiVector> nullspace = coarseLevel.Get<RCP<MultiVector> >("Nullspace", (*it)->GetFactory("Nullspace").get());
      coarseLevel.Set<RCP<MultiVector> >("Nullspace", nullspace, (*it)->GetFactory("Nullspace").get());
    }

    ////////////////////////////////////////////////////////////

    curBlockId++;
  }  // end for loop

  // extract map index base from maps of blocked P
  GO rangeIndexBase  = originalTransferOp->getRangeMap()->getIndexBase();
  GO domainIndexBase = originalTransferOp->getDomainMap()->getIndexBase();

  // check this
  Teuchos::ArrayView<GO> fullRangeMapGIDs(fullRangeMapVector.size() ? &fullRangeMapVector[0] : 0, fullRangeMapVector.size());
  Teuchos::RCP<const StridedMap> stridedRgFullMap = Teuchos::rcp_dynamic_cast<const StridedMap>(rangeMapExtractor->getFullMap());
  Teuchos::RCP<const Map> fullRangeMap            = Teuchos::null;
  if (stridedRgFullMap != Teuchos::null) {
    std::vector<size_t> stridedData = stridedRgFullMap->getStridingData();
    fullRangeMap =
        StridedMapFactory::Build(
            originalTransferOp->getRangeMap()->lib(),
            Teuchos::OrdinalTraits<Xpetra::global_size_t>::invalid(),
            fullRangeMapGIDs,
            rangeIndexBase,
            stridedData,
            originalTransferOp->getRangeMap()->getComm(),
            stridedRgFullMap->getStridedBlockId(),
            stridedRgFullMap->getOffset());
  } else {
    fullRangeMap =
        MapFactory::Build(
            originalTransferOp->getRangeMap()->lib(),
            Teuchos::OrdinalTraits<Xpetra::global_size_t>::invalid(),
            fullRangeMapGIDs,
            rangeIndexBase,
            originalTransferOp->getRangeMap()->getComm());
  }

  Teuchos::ArrayView<GO> fullDomainMapGIDs(fullDomainMapVector.size() ? &fullDomainMapVector[0] : 0, fullDomainMapVector.size());
  Teuchos::RCP<const StridedMap> stridedDoFullMap = Teuchos::rcp_dynamic_cast<const StridedMap>(domainMapExtractor->getFullMap());
  Teuchos::RCP<const Map> fullDomainMap           = Teuchos::null;
  if (stridedDoFullMap != Teuchos::null) {
    TEUCHOS_TEST_FOR_EXCEPTION(stridedDoFullMap == Teuchos::null, Exceptions::BadCast, "MueLu::BlockedPFactory::Build: full map in domain map extractor has no striding information! error.");
    std::vector<size_t> stridedData2 = stridedDoFullMap->getStridingData();
    fullDomainMap =
        StridedMapFactory::Build(
            originalTransferOp->getDomainMap()->lib(),
            Teuchos::OrdinalTraits<Xpetra::global_size_t>::invalid(),
            fullDomainMapGIDs,
            domainIndexBase,
            stridedData2,
            originalTransferOp->getDomainMap()->getComm(),
            stridedDoFullMap->getStridedBlockId(),
            stridedDoFullMap->getOffset());
  } else {
    fullDomainMap =
        MapFactory::Build(
            originalTransferOp->getDomainMap()->lib(),
            Teuchos::OrdinalTraits<Xpetra::global_size_t>::invalid(),
            fullDomainMapGIDs,
            domainIndexBase,
            originalTransferOp->getDomainMap()->getComm());
  }

  if (bRestrictComm) {
    fullRangeMap->removeEmptyProcesses();
    fullDomainMap->removeEmptyProcesses();
  }

  // build map extractors
  Teuchos::RCP<const MapExtractor> rebrangeMapExtractor =
      MapExtractorFactory::Build(fullRangeMap, subBlockRRangeMaps, bThyraRangeGIDs);
  Teuchos::RCP<const MapExtractor> rebdomainMapExtractor =
      MapExtractorFactory::Build(fullDomainMap, subBlockRDomainMaps, bThyraDomainGIDs);

  Teuchos::RCP<BlockedCrsMatrix> bRebR = Teuchos::rcp(new BlockedCrsMatrix(rebrangeMapExtractor, rebdomainMapExtractor, 10));
  for (size_t i = 0; i < subBlockRRangeMaps.size(); i++) {
    Teuchos::RCP<CrsMatrixWrap> crsOpii = Teuchos::rcp_dynamic_cast<CrsMatrixWrap>(subBlockRebR[i]);
    bRebR->setMatrix(i, i, crsOpii);
  }

  bRebR->fillComplete();

  Set(coarseLevel, "R", Teuchos::rcp_dynamic_cast<Matrix>(bRebR));  // do nothing  // TODO remove this!

}  // Build

}  // namespace MueLu

#endif /* MUELU_REBALANCEBLOCKRESTRICTIONFACTORY_DEF_HPP_ */
