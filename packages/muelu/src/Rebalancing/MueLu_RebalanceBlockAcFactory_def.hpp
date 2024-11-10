// @HEADER
// *****************************************************************************
//        MueLu: A package for multigrid based preconditioning
//
// Copyright 2012 NTESS and the MueLu contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef MUELU_REBALANCEBLOCKACFACTORY_DEF_HPP_
#define MUELU_REBALANCEBLOCKACFACTORY_DEF_HPP_

#include <Xpetra_Matrix.hpp>
#include <Xpetra_CrsMatrix.hpp>
#include <Xpetra_CrsMatrixWrap.hpp>
#include <Xpetra_MatrixFactory.hpp>
#include <Xpetra_MapExtractor.hpp>
#include <Xpetra_MapExtractorFactory.hpp>
#include <Xpetra_StridedMap.hpp>
#include <Xpetra_StridedMapFactory.hpp>
#include <Xpetra_BlockedCrsMatrix.hpp>

#include <Xpetra_VectorFactory.hpp>

#include "MueLu_RebalanceBlockAcFactory_decl.hpp"

#include "MueLu_FactoryManagerBase.hpp"
#include "MueLu_HierarchyUtils.hpp"
#include "MueLu_MasterList.hpp"
#include "MueLu_Monitor.hpp"
#include "MueLu_PerfUtils.hpp"

namespace MueLu {

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
RCP<const ParameterList> RebalanceBlockAcFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::GetValidParameterList() const {
  RCP<ParameterList> validParamList = rcp(new ParameterList());

#define SET_VALID_ENTRY(name) validParamList->setEntry(name, MasterList::getEntry(name))
  SET_VALID_ENTRY("repartition: use subcommunicators");
#undef SET_VALID_ENTRY

  validParamList->set<RCP<const FactoryBase> >("A", Teuchos::null, "Generating factory of the matrix A for rebalancing");
  validParamList->set<RCP<const FactoryBase> >("Importer", Teuchos::null, "Generating factory of the matrix Importer for rebalancing");
  validParamList->set<RCP<const FactoryBase> >("SubImporters", Teuchos::null, "Generating factory of the matrix sub-Importers for rebalancing");

  return validParamList;
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void RebalanceBlockAcFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::AddFactoryManager(RCP<const FactoryManagerBase> FactManager) {
  FactManager_.push_back(FactManager);
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void RebalanceBlockAcFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::DeclareInput(Level &fineLevel, Level &coarseLevel) const {
  Input(coarseLevel, "A");

  std::vector<Teuchos::RCP<const FactoryManagerBase> >::const_iterator it;
  for (it = FactManager_.begin(); it != FactManager_.end(); ++it) {
    SetFactoryManager fineSFM(rcpFromRef(fineLevel), *it);
    SetFactoryManager coarseSFM(rcpFromRef(coarseLevel), *it);

    coarseLevel.DeclareInput("Importer", (*it)->GetFactory("Importer").get(), this);
  }

  // Use the non-manager path if the maps / importers are generated in one place
  if (FactManager_.size() == 0) {
    Input(coarseLevel, "SubImporters");
  }
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void RebalanceBlockAcFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Build(Level &fineLevel, Level &coarseLevel) const {
  FactoryMonitor m(*this, "Computing blocked Ac", coarseLevel);

  RCP<Teuchos::FancyOStream> out = Teuchos::fancyOStream(Teuchos::rcpFromRef(std::cout));

  RCP<Matrix> originalAc = Get<RCP<Matrix> >(coarseLevel, "A");

  RCP<Xpetra::BlockedCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> > bA = Teuchos::rcp_dynamic_cast<Xpetra::BlockedCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> >(originalAc);
  TEUCHOS_TEST_FOR_EXCEPTION(bA == Teuchos::null, Exceptions::BadCast, "MueLu::RebalanceBlockAcFactory::Build: input matrix A is not of type BlockedCrsMatrix! error.");
  TEUCHOS_TEST_FOR_EXCEPTION(bA->Rows() != bA->Cols(), Exceptions::RuntimeError, "MueLu::RebalanceBlockAcFactory::Build: Blocked operator has " << bA->Rows() << " and " << bA->Cols() << ". We only support square matrices (with same number of blocks and columns).");

  // Variables to set up map extractors for blocked operators
  std::vector<GO> fullRangeMapVector;
  std::vector<GO> fullDomainMapVector;
  std::vector<RCP<const Map> > subBlockARangeMaps;
  std::vector<RCP<const Map> > subBlockADomainMaps;
  subBlockARangeMaps.reserve(bA->Rows());
  subBlockADomainMaps.reserve(bA->Cols());

  // store map extractors
  RCP<const MapExtractor> rangeMapExtractor  = bA->getRangeMapExtractor();
  RCP<const MapExtractor> domainMapExtractor = bA->getDomainMapExtractor();

  // check if GIDs for full maps have to be sorted:
  // For the Thyra mode ordering they do not have to be sorted since the GIDs are
  // numbered as 0...n1,0...,n2 (starting with zero for each subblock). The MapExtractor
  // generates unique GIDs during the construction.
  // For Xpetra style, the GIDs have to be reordered. Such that one obtains a ordered
  // list of GIDs in an increasing ordering. In Xpetra, the GIDs are all unique through
  // out all submaps.
  bool bThyraRangeGIDs  = rangeMapExtractor->getThyraMode();
  bool bThyraDomainGIDs = domainMapExtractor->getThyraMode();

  // vector containing rebalanced blocks (for final output)
  std::vector<RCP<Matrix> > subBlockRebA =
      std::vector<RCP<Matrix> >(bA->Cols() * bA->Rows(), Teuchos::null);

  // vector with Import objects from the different
  // RepartitionFactory instances
  std::vector<RCP<const Import> > importers = std::vector<RCP<const Import> >(bA->Rows(), Teuchos::null);
  std::vector<Teuchos::RCP<const FactoryManagerBase> >::const_iterator it;
  size_t idx = 0;
  for (it = FactManager_.begin(); it != FactManager_.end(); ++it) {
    SetFactoryManager fineSFM(rcpFromRef(fineLevel), *it);
    SetFactoryManager coarseSFM(rcpFromRef(coarseLevel), *it);

    RCP<const Import> rebalanceImporter = coarseLevel.Get<RCP<const Import> >("Importer", (*it)->GetFactory("Importer").get());
    importers[idx]                      = rebalanceImporter;
    idx++;
  }
  if (FactManager_.size() == 0) {
    importers = Get<std::vector<RCP<const Import> > >(coarseLevel, "SubImporters");
  }

  // restrict communicator?
  bool bRestrictComm      = false;
  const ParameterList &pL = GetParameterList();
  if (pL.get<bool>("repartition: use subcommunicators") == true)
    bRestrictComm = true;

  RCP<ParameterList> XpetraList = Teuchos::rcp(new ParameterList());
  if (bRestrictComm)
    XpetraList->set("Restrict Communicator", true);
  else
    XpetraList->set("Restrict Communicator", false);

  // communicator for final (rebalanced) operator.
  // If communicator is not restricted it should be identical to the communicator in bA
  RCP<const Teuchos::Comm<int> > rebalancedComm = Teuchos::null;

  // loop through all blocks and rebalance blocks
  // Note: so far we do not support rebalancing of nested operators
  //       TODO add a check for this
  for (size_t i = 0; i < bA->Rows(); i++) {
    for (size_t j = 0; j < bA->Cols(); j++) {
      // extract matrix block
      RCP<Matrix> Aij = bA->getMatrix(i, j);

      std::stringstream ss;
      ss << "Rebalancing matrix block A(" << i << "," << j << ")";
      SubFactoryMonitor subM(*this, ss.str(), coarseLevel);

      RCP<Matrix> rebAij = Teuchos::null;
      // General rebalancing
      if (importers[i] != Teuchos::null &&
          importers[j] != Teuchos::null &&
          Aij != Teuchos::null) {
        RCP<const Map> targetRangeMap  = importers[i]->getTargetMap();
        RCP<const Map> targetDomainMap = importers[j]->getTargetMap();

        // Copy the block Aij
        // TAW: Do we need a copy or can we do in-place rebalancing?
        //      If we do in-place rebalancing the original distribution is lost
        //      We don't really need it any more, though.
        // RCP<Matrix> cAij = MatrixFactory::BuildCopy(Aij);
        RCP<Matrix> cAij = Aij;  // do not copy the matrix data (just an rcp pointer)

        // create a new importer for column map needed for rebalanced Aij
        Teuchos::RCP<const Import> rebAijImport = ImportFactory::Build(importers[j]->getTargetMap(), cAij->getColMap());
        TEUCHOS_TEST_FOR_EXCEPTION(rebAijImport.is_null() == true, Exceptions::RuntimeError, "MueLu::RebalanceBlockAcFactory::Build: Importer associated with block " << j << " is null.");

        Teuchos::RCP<const CrsMatrixWrap> cAwij = Teuchos::rcp_dynamic_cast<const CrsMatrixWrap>(cAij);
        TEUCHOS_TEST_FOR_EXCEPTION(cAwij.is_null() == true, Exceptions::RuntimeError, "MueLu::RebalanceBlockAcFactory::Build: Block (" << i << "," << j << ") is not of type CrsMatrix. We cannot rebalanced (nested) operators.");
        Teuchos::RCP<CrsMatrix> cAmij = cAwij->getCrsMatrix();

        // change domain map to rebalanced domain map (in-place). Update the importer to represent the column map
        // cAmij->replaceDomainMapAndImporter(importers[j]->getTargetMap(),rebAijImport);

        // rebalance rows of matrix block. Don't change the domain map (-> Teuchos::null)
        // NOTE: If the communicator is restricted away, Build returns Teuchos::null.
        rebAij = MatrixFactory::Build(cAij, *(importers[i]), *(importers[j]), targetDomainMap, targetRangeMap, XpetraList);
      }  // rebalance matrix block A(i,i)
      else {
        rebAij = Aij;  // no rebalancing or empty block!
      }

      // store new block in output
      subBlockRebA[i * bA->Cols() + j] = rebAij;

      if (!rebAij.is_null()) {
        // store communicator
        if (rebalancedComm.is_null()) rebalancedComm = rebAij->getRowMap()->getComm();

        // printout rebalancing information
        RCP<ParameterList> params = rcp(new ParameterList());
        params->set("printLoadBalancingInfo", true);
        std::stringstream ss2;
        ss2 << "A(" << i << "," << j << ") rebalanced:";
        GetOStream(Statistics0) << PerfUtils::PrintMatrixInfo(*rebAij, ss2.str(), params);
      }
    }  // loop over columns j

    // fix striding information of diagonal blocks
    // Note: we do not care about the off-diagonal blocks. We just make sure, that the
    // diagonal blocks have the corresponding striding information from the map extractors
    // Note: the diagonal block never should be zero.
    // TODO what if a diagonal block is Teuchos::null?
    if (subBlockRebA[i * bA->Cols() + i].is_null() == false) {
      RCP<Matrix> rebAii                               = subBlockRebA[i * bA->Cols() + i];
      Teuchos::RCP<const StridedMap> orig_stridedRgMap = Teuchos::rcp_dynamic_cast<const StridedMap>(rangeMapExtractor->getMap(i, rangeMapExtractor->getThyraMode()));
      Teuchos::RCP<const Map> stridedRgMap             = Teuchos::null;
      if (orig_stridedRgMap != Teuchos::null) {
        std::vector<size_t> stridingData                       = orig_stridedRgMap->getStridingData();
        Teuchos::ArrayView<const GlobalOrdinal> nodeRangeMapii = rebAii->getRangeMap()->getLocalElementList();
        stridedRgMap                                           = StridedMapFactory::Build(
                                                      bA->getRangeMap()->lib(),
                                                      Teuchos::OrdinalTraits<Xpetra::global_size_t>::invalid(),
                                                      nodeRangeMapii,
                                                      rebAii->getRangeMap()->getIndexBase(),
                                                      stridingData,
                                                      rebalancedComm,
                                                      orig_stridedRgMap->getStridedBlockId(),
                                                      orig_stridedRgMap->getOffset());
      }
      Teuchos::RCP<const StridedMap> orig_stridedDoMap = Teuchos::rcp_dynamic_cast<const StridedMap>(domainMapExtractor->getMap(i, domainMapExtractor->getThyraMode()));
      Teuchos::RCP<const Map> stridedDoMap             = Teuchos::null;
      if (orig_stridedDoMap != Teuchos::null) {
        std::vector<size_t> stridingData                        = orig_stridedDoMap->getStridingData();
        Teuchos::ArrayView<const GlobalOrdinal> nodeDomainMapii = rebAii->getDomainMap()->getLocalElementList();
        stridedDoMap                                            = StridedMapFactory::Build(
                                                       bA->getDomainMap()->lib(),
                                                       Teuchos::OrdinalTraits<Xpetra::global_size_t>::invalid(),
                                                       nodeDomainMapii,
                                                       rebAii->getDomainMap()->getIndexBase(),
                                                       stridingData,
                                                       rebalancedComm,
                                                       orig_stridedDoMap->getStridedBlockId(),
                                                       orig_stridedDoMap->getOffset());
      }

      if (bRestrictComm) {
        stridedRgMap->removeEmptyProcesses();
        stridedDoMap->removeEmptyProcesses();
      }

      TEUCHOS_TEST_FOR_EXCEPTION(stridedRgMap == Teuchos::null, Exceptions::RuntimeError, "MueLu::RebalanceBlockAcFactory::Build: failed to generate striding information. error.");
      TEUCHOS_TEST_FOR_EXCEPTION(stridedDoMap == Teuchos::null, Exceptions::RuntimeError, "MueLu::RebalanceBlockAcFactory::Build: failed to generate striding information. error.");

      // replace stridedMaps view in diagonal sub block
      if (rebAii->IsView("stridedMaps")) rebAii->RemoveView("stridedMaps");
      rebAii->CreateView("stridedMaps", stridedRgMap, stridedDoMap);
      // collect Xpetra-based global row ids for map extractors
      subBlockARangeMaps.push_back(rebAii->getRowMap("stridedMaps"));
      Teuchos::ArrayView<const GlobalOrdinal> nodeRangeMap = rebAii->getRangeMap()->getLocalElementList();
      // append the GIDs in the end. Do not sort if we have Thyra style GIDs
      fullRangeMapVector.insert(fullRangeMapVector.end(), nodeRangeMap.begin(), nodeRangeMap.end());
      if (bThyraRangeGIDs == false)
        sort(fullRangeMapVector.begin(), fullRangeMapVector.end());

      subBlockADomainMaps.push_back(rebAii->getColMap("stridedMaps"));
      Teuchos::ArrayView<const GlobalOrdinal> nodeDomainMap = rebAii->getDomainMap()->getLocalElementList();
      // append the GIDs in the end. Do not sort if we have Thyra style GIDs
      fullDomainMapVector.insert(fullDomainMapVector.end(), nodeDomainMap.begin(), nodeDomainMap.end());
      if (bThyraDomainGIDs == false)
        sort(fullDomainMapVector.begin(), fullDomainMapVector.end());
    }  // end if rebAii != Teuchos::null
  }    // loop over rows i

  // all sub blocks are rebalanced (if available)

  // Short cut if this processor is not in the list of active processors
  if (rebalancedComm == Teuchos::null) {
    GetOStream(Debug, -1) << "RebalanceBlockedAc: deactivate proc " << originalAc->getRowMap()->getComm()->getRank() << std::endl;
    // TAW: it is important that we create a dummy object of type BlockedCrsMatrix (even if we set it to Teuchos::null)
    Teuchos::RCP<BlockedCrsMatrix> reb_bA = Teuchos::null;
    coarseLevel.Set("A", Teuchos::rcp_dynamic_cast<Matrix>(reb_bA), this);
    return;
  }

  // now, subBlockRebA contains all rebalanced matrix blocks
  // extract map index base from maps of blocked A
  GO rangeIndexBase  = bA->getRangeMap()->getIndexBase();
  GO domainIndexBase = bA->getDomainMap()->getIndexBase();

  Teuchos::ArrayView<GO> fullRangeMapGIDs(fullRangeMapVector.size() ? &fullRangeMapVector[0] : 0, fullRangeMapVector.size());
  Teuchos::RCP<const StridedMap> stridedRgFullMap = Teuchos::rcp_dynamic_cast<const StridedMap>(rangeMapExtractor->getFullMap());
  Teuchos::RCP<const Map> fullRangeMap            = Teuchos::null;
  if (stridedRgFullMap != Teuchos::null) {
    std::vector<size_t> stridedData = stridedRgFullMap->getStridingData();
    fullRangeMap =
        StridedMapFactory::Build(
            bA->getRangeMap()->lib(),
            Teuchos::OrdinalTraits<Xpetra::global_size_t>::invalid(),
            fullRangeMapGIDs,
            rangeIndexBase,
            stridedData,
            rebalancedComm,
            stridedRgFullMap->getStridedBlockId(),
            stridedRgFullMap->getOffset());
  } else {
    fullRangeMap =
        MapFactory::Build(
            bA->getRangeMap()->lib(),
            Teuchos::OrdinalTraits<Xpetra::global_size_t>::invalid(),
            fullRangeMapGIDs,
            rangeIndexBase,
            rebalancedComm);
  }
  Teuchos::ArrayView<GO> fullDomainMapGIDs(fullDomainMapVector.size() ? &fullDomainMapVector[0] : 0, fullDomainMapVector.size());

  Teuchos::RCP<const StridedMap> stridedDoFullMap = Teuchos::rcp_dynamic_cast<const StridedMap>(domainMapExtractor->getFullMap());
  Teuchos::RCP<const Map> fullDomainMap           = Teuchos::null;
  if (stridedDoFullMap != Teuchos::null) {
    TEUCHOS_TEST_FOR_EXCEPTION(stridedDoFullMap == Teuchos::null, Exceptions::BadCast, "MueLu::RebalanceBlockedAc::Build: full map in domain map extractor has no striding information! error.");
    std::vector<size_t> stridedData2 = stridedDoFullMap->getStridingData();
    fullDomainMap =
        StridedMapFactory::Build(
            bA->getDomainMap()->lib(),
            Teuchos::OrdinalTraits<Xpetra::global_size_t>::invalid(),
            fullDomainMapGIDs,
            domainIndexBase,
            stridedData2,
            rebalancedComm,
            stridedDoFullMap->getStridedBlockId(),
            stridedDoFullMap->getOffset());
  } else {
    fullDomainMap =
        MapFactory::Build(
            bA->getDomainMap()->lib(),
            Teuchos::OrdinalTraits<Xpetra::global_size_t>::invalid(),
            fullDomainMapGIDs,
            domainIndexBase,
            rebalancedComm);
  }

  if (bRestrictComm) {
    fullRangeMap->removeEmptyProcesses();
    fullDomainMap->removeEmptyProcesses();
  }

  // build map extractors
  RCP<const MapExtractor> rebRangeMapExtractor  = MapExtractorFactory::Build(fullRangeMap, subBlockARangeMaps, bThyraRangeGIDs);
  RCP<const MapExtractor> rebDomainMapExtractor = MapExtractorFactory::Build(fullDomainMap, subBlockADomainMaps, bThyraDomainGIDs);

  TEUCHOS_TEST_FOR_EXCEPTION(rangeMapExtractor->NumMaps() != rebRangeMapExtractor->NumMaps(), Exceptions::RuntimeError,
                             "MueLu::RebalanceBlockedAc::Build: Rebalanced RangeMapExtractor has " << rebRangeMapExtractor->NumMaps()
                                                                                                   << " sub maps. Original RangeMapExtractor has " << rangeMapExtractor->NumMaps() << ". They must match!");
  TEUCHOS_TEST_FOR_EXCEPTION(domainMapExtractor->NumMaps() != rebDomainMapExtractor->NumMaps(), Exceptions::RuntimeError,
                             "MueLu::RebalanceBlockedAc::Build: Rebalanced DomainMapExtractor has " << rebDomainMapExtractor->NumMaps()
                                                                                                    << " sub maps. Original DomainMapExtractor has " << domainMapExtractor->NumMaps() << ". They must match!");

  Teuchos::RCP<BlockedCrsMatrix> reb_bA = Teuchos::rcp(new BlockedCrsMatrix(rebRangeMapExtractor, rebDomainMapExtractor, 10));
  for (size_t i = 0; i < bA->Rows(); i++) {
    for (size_t j = 0; j < bA->Cols(); j++) {
      reb_bA->setMatrix(i, j, subBlockRebA[i * bA->Cols() + j]);
    }
  }

  reb_bA->fillComplete();

  coarseLevel.Set("A", Teuchos::rcp_dynamic_cast<Matrix>(reb_bA), this);
  // rebalance additional data:
  // be aware, that we just call the rebalance factories without switching to local
  // factory managers, i.e. the rebalance factories have to be defined with the appropriate
  // factories by the user!
  if (rebalanceFacts_.begin() != rebalanceFacts_.end()) {
    SubFactoryMonitor m2(*this, "Rebalance additional data", coarseLevel);

    // call Build of all user-given transfer factories
    for (std::vector<RCP<const FactoryBase> >::const_iterator it2 = rebalanceFacts_.begin(); it2 != rebalanceFacts_.end(); ++it2) {
      GetOStream(Runtime0) << "RebalanceBlockedAc: call rebalance factory " << (*it2).get() << ": " << (*it2)->description() << std::endl;
      (*it2)->CallBuild(coarseLevel);
    }
  }
}  // Build()

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void RebalanceBlockAcFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::AddRebalanceFactory(const RCP<const FactoryBase> &factory) {
  rebalanceFacts_.push_back(factory);
}  // AddRebalanceFactory()

}  // namespace MueLu

#endif /* MUELU_REBALANCEBLOCKACFACTORY_DEF_HPP_ */
