// @HEADER
// *****************************************************************************
//        MueLu: A package for multigrid based preconditioning
//
// Copyright 2012 NTESS and the MueLu contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef MUELU_REBALANCEBLOCKINTERPOLATIONFACTORY_DEF_HPP_
#define MUELU_REBALANCEBLOCKINTERPOLATIONFACTORY_DEF_HPP_

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

#include "MueLu_RebalanceBlockInterpolationFactory_decl.hpp"

#include "MueLu_FactoryManagerBase.hpp"
#include "MueLu_HierarchyUtils.hpp"
#include "MueLu_Level.hpp"
#include "MueLu_Monitor.hpp"
#include "MueLu_PerfUtils.hpp"

namespace MueLu {

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
RCP<const ParameterList> RebalanceBlockInterpolationFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::GetValidParameterList() const {
  RCP<ParameterList> validParamList = rcp(new ParameterList());

  validParamList->set<RCP<const FactoryBase> >("P", Teuchos::null, "Factory of the prolongation operator that need to be rebalanced (only used if type=Interpolation)");
  validParamList->set<RCP<const FactoryBase> >("A", Teuchos::null, "Factory for generating the non-rebalanced coarse level A. We need this to make sure the non-rebalanced coarse A is calculated first before rebalancing takes place.");

  validParamList->set<RCP<const FactoryBase> >("Coordinates", Teuchos::null, "Factory for generating the non-rebalanced Coordinates.");
  validParamList->set<RCP<const FactoryBase> >("Importer", Teuchos::null, "Generating factory of the matrix Importer for rebalancing");
  validParamList->set<RCP<const FactoryBase> >("SubImporters", Teuchos::null, "Generating factory of the matrix sub-Importers for rebalancing");

#define SET_VALID_ENTRY(name) validParamList->setEntry(name, MasterList::getEntry(name))
  // SET_VALID_ENTRY("repartition: use subcommunicators");
#undef SET_VALID_ENTRY

  // TODO validation: "P" parameter valid only for type="Interpolation" and "R" valid only for type="Restriction". Like so:
  // if (paramList.isEntry("type") && paramList.get("type) == "Interpolation) {
  //     validParamList->set< RCP<const FactoryBase> >("P",              Teuchos::null, "Factory of the prolongation operator that need to be rebalanced (only used if type=Interpolation)");

  return validParamList;
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void RebalanceBlockInterpolationFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::AddFactoryManager(RCP<const FactoryManagerBase> FactManager) {
  FactManager_.push_back(FactManager);
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void RebalanceBlockInterpolationFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::DeclareInput(Level &fineLevel, Level &coarseLevel) const {
  Input(coarseLevel, "P");
  Input(coarseLevel, "A");  // we request the non-rebalanced coarse level A since we have to make sure it is calculated before rebalancing starts!

  std::vector<Teuchos::RCP<const FactoryManagerBase> >::const_iterator it;
  for (it = FactManager_.begin(); it != FactManager_.end(); ++it) {
    SetFactoryManager fineSFM(rcpFromRef(fineLevel), *it);
    SetFactoryManager coarseSFM(rcpFromRef(coarseLevel), *it);

    // Request Importer and Coordinates (if defined in xml file)
    // Note, that we have to use the Level::DeclareInput routine in order to use the FactoryManager *it (rather than the main factory manager)
    coarseLevel.DeclareInput("Importer", (*it)->GetFactory("Importer").get(), this);
    if ((*it)->hasFactory("Coordinates") == true)
      coarseLevel.DeclareInput("Coordinates", (*it)->GetFactory("Coordinates").get(), this);
  }

  // Use the non-manager path if the maps / importers are generated in one place
  if (FactManager_.size() == 0) {
    Input(coarseLevel, "Importer");
    Input(coarseLevel, "SubImporters");
    Input(coarseLevel, "Coordinates");
  }
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void RebalanceBlockInterpolationFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Build(Level &fineLevel, Level &coarseLevel) const {
  FactoryMonitor m(*this, "Build", coarseLevel);
  typedef Xpetra::MultiVector<typename Teuchos::ScalarTraits<Scalar>::magnitudeType, LO, GO, NO> xdMV;
  typedef Xpetra::BlockedMultiVector<typename Teuchos::ScalarTraits<Scalar>::magnitudeType, LO, GO, NO> xdBV;

  bool UseSingleSource = FactManager_.size() == 0;
  // RCP<Teuchos::FancyOStream> out = Teuchos::fancyOStream(Teuchos::rcpFromRef(std::cout));

  Teuchos::RCP<Matrix> nonrebCoarseA = Get<RCP<Matrix> >(coarseLevel, "A");

  Teuchos::RCP<Matrix> originalTransferOp = Teuchos::null;
  originalTransferOp                      = Get<RCP<Matrix> >(coarseLevel, "P");

  RCP<Xpetra::BlockedCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> > bOriginalTransferOp =
      Teuchos::rcp_dynamic_cast<Xpetra::BlockedCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> >(originalTransferOp);
  TEUCHOS_TEST_FOR_EXCEPTION(bOriginalTransferOp == Teuchos::null, Exceptions::BadCast, "MueLu::RebalanceBlockTransferFactory::Build: input matrix P or R is not of type BlockedCrsMatrix! error.");

  RCP<const MapExtractor> rangeMapExtractor  = bOriginalTransferOp->getRangeMapExtractor();
  RCP<const MapExtractor> domainMapExtractor = bOriginalTransferOp->getDomainMapExtractor();

  // check if GIDs for full maps have to be sorted:
  // For the Thyra mode ordering they do not have to be sorted since the GIDs are
  // numbered as 0...n1,0...,n2 (starting with zero for each subblock). The MapExtractor
  // generates unique GIDs during the construction.
  // For Xpetra style, the GIDs have to be reordered. Such that one obtains a ordered
  // list of GIDs in an increasing ordering. In Xpetra, the GIDs are all unique through
  // out all submaps.
  bool bThyraRangeGIDs  = rangeMapExtractor->getThyraMode();
  bool bThyraDomainGIDs = domainMapExtractor->getThyraMode();

  // declare variables for maps of blocked rebalanced prolongation operator
  std::vector<GO> fullRangeMapVector;   // contains all range GIDs on current processor
  std::vector<GO> fullDomainMapVector;  // contains all domain GIDs on current processor
  std::vector<RCP<const Map> > subBlockPRangeMaps;
  std::vector<RCP<const Map> > subBlockPDomainMaps;
  subBlockPRangeMaps.reserve(bOriginalTransferOp->Rows());   // reserve size for block P operators
  subBlockPDomainMaps.reserve(bOriginalTransferOp->Cols());  // reserve size for block P operators

  std::vector<Teuchos::RCP<Matrix> > subBlockRebP;
  subBlockRebP.reserve(bOriginalTransferOp->Rows());

  // For use in single-source mode only
  std::vector<RCP<const Import> > importers = std::vector<RCP<const Import> >(bOriginalTransferOp->Rows(), Teuchos::null);
  std::vector<RCP<xdMV> > newCoordinates(bOriginalTransferOp->Rows());
  RCP<xdBV> oldCoordinates;
  if (UseSingleSource) {
    importers      = Get<std::vector<RCP<const Import> > >(coarseLevel, "SubImporters");
    oldCoordinates = Get<RCP<xdBV> >(coarseLevel, "Coordinates");
  }

  int curBlockId                               = 0;
  Teuchos::RCP<const Import> rebalanceImporter = Teuchos::null;
  std::vector<Teuchos::RCP<const FactoryManagerBase> >::const_iterator it;
  for (it = FactManager_.begin(); it != FactManager_.end(); ++it) {
    // begin SubFactoryManager environment
    SetFactoryManager fineSFM(rcpFromRef(fineLevel), *it);
    SetFactoryManager coarseSFM(rcpFromRef(coarseLevel), *it);

    // TAW: use the Level::Get routine in order to access the data declared in (*it) factory manager (rather than the main factory manager)
    if (UseSingleSource)
      rebalanceImporter = importers[curBlockId];
    else
      rebalanceImporter = coarseLevel.Get<Teuchos::RCP<const Import> >("Importer", (*it)->GetFactory("Importer").get());

    // extract diagonal matrix block
    Teuchos::RCP<Matrix> Pii         = bOriginalTransferOp->getMatrix(curBlockId, curBlockId);
    Teuchos::RCP<CrsMatrixWrap> Pwii = Teuchos::rcp_dynamic_cast<CrsMatrixWrap>(Pii);
    TEUCHOS_TEST_FOR_EXCEPTION(Pwii == Teuchos::null, Xpetra::Exceptions::BadCast, "MueLu::RebalanceBlockTransferFactory::Build: block " << curBlockId << " is not of type CrsMatrixWrap. We need an underlying CsrMatrix to replace domain map and importer!");

    MUELU_TEST_FOR_EXCEPTION(bThyraRangeGIDs == true && Pii->getRowMap()->getMinAllGlobalIndex() != 0,
                             Exceptions::RuntimeError,
                             "MueLu::RebalanceBlockInterpolationFactory::Build: inconsistent Thyra GIDs. Thyra global ids for block range " << curBlockId << " start with " << Pii->getRowMap()->getMinAllGlobalIndex() << " but should start with 0");
    MUELU_TEST_FOR_EXCEPTION(bThyraDomainGIDs == true && Pii->getColMap()->getMinAllGlobalIndex() != 0,
                             Exceptions::RuntimeError,
                             "MueLu::RebalanceBlockInterpolationFactory::Build: inconsistent Thyra GIDs. Thyra global ids for block domain " << curBlockId << " start with " << Pii->getColMap()->getMinAllGlobalIndex() << " but should start with 0");

    // rebalance P11
    if (rebalanceImporter != Teuchos::null) {
      std::stringstream ss;
      ss << "Rebalancing prolongator block P(" << curBlockId << "," << curBlockId << ")";
      SubFactoryMonitor m1(*this, ss.str(), coarseLevel);

      // P is the transfer operator from the coarse grid to the fine grid.
      // P must transfer the data from the newly reordered coarse A to the (unchanged) fine A.
      // This means that the domain map (coarse) of P must be changed according to the new partition. The range map (fine) is kept unchanged.
      //
      // The domain map of P must match the range map of R.
      // See also note below about domain/range map of R and its implications for P.
      //
      // To change the domain map of P, P needs to be fillCompleted again with the new domain map.
      // To achieve this, P is copied into a new matrix that is not fill-completed.
      // The doImport() operation is just used here to make a copy of P: the importer is trivial and there is no data movement involved.
      // The reordering actually happens during the fillComplete() with domainMap == rebalanceImporter->getTargetMap().
      RCP<const Import> newImporter;
      {
        SubFactoryMonitor subM(*this, "Rebalancing prolongator  -- fast map replacement", coarseLevel);
        newImporter = ImportFactory::Build(rebalanceImporter->getTargetMap(), Pii->getColMap());
        Pwii->getCrsMatrix()->replaceDomainMapAndImporter(rebalanceImporter->getTargetMap(), newImporter);
      }

      RCP<ParameterList> params = rcp(new ParameterList());
      params->set("printLoadBalancingInfo", true);
      std::stringstream ss2;
      ss2 << "P(" << curBlockId << "," << curBlockId << ") rebalanced:";
      GetOStream(Statistics0) << PerfUtils::PrintMatrixInfo(*Pii, ss2.str(), params);

      // store rebalanced P block
      subBlockRebP.push_back(Pii);

      // rebalance coordinates
      // TAW: Note, that each sub-block manager overwrites the Coordinates. So far we only support one set of Coordinates
      //      for a multiphysics problem (i.e., we only support volume coupled problems with the same mesh)
      if (UseSingleSource) {
        RCP<xdMV> localCoords = oldCoordinates->getMultiVector(curBlockId);

        // FIXME: This should be extended to work with blocking
        RCP<const Import> coordImporter = rebalanceImporter;

        RCP<xdMV> permutedLocalCoords = Xpetra::MultiVectorFactory<typename Teuchos::ScalarTraits<Scalar>::magnitudeType, LO, GO, NO>::Build(coordImporter->getTargetMap(), localCoords->getNumVectors());
        permutedLocalCoords->doImport(*localCoords, *coordImporter, Xpetra::INSERT);

        newCoordinates[curBlockId] = permutedLocalCoords;
      } else if ((*it)->hasFactory("Coordinates") == true && coarseLevel.IsAvailable("Coordinates", (*it)->GetFactory("Coordinates").get()) == true) {
        RCP<xdMV> coords = coarseLevel.Get<RCP<xdMV> >("Coordinates", (*it)->GetFactory("Coordinates").get());

        // This line must be after the Get call
        SubFactoryMonitor subM(*this, "Rebalancing coordinates", coarseLevel);

        LO nodeNumElts = coords->getMap()->getLocalNumElements();

        // If a process has no matrix rows, then we can't calculate blocksize using the formula below.
        LO myBlkSize = 0, blkSize = 0;

        if (nodeNumElts > 0) {
          MUELU_TEST_FOR_EXCEPTION(rebalanceImporter->getSourceMap()->getLocalNumElements() % nodeNumElts != 0,
                                   Exceptions::RuntimeError,
                                   "MueLu::RebalanceBlockInterpolationFactory::Build: block size. " << rebalanceImporter->getSourceMap()->getLocalNumElements() << " not divisable by " << nodeNumElts);
          myBlkSize = rebalanceImporter->getSourceMap()->getLocalNumElements() / nodeNumElts;
        }

        MueLu_maxAll(coords->getMap()->getComm(), myBlkSize, blkSize);

        RCP<const Import> coordImporter = Teuchos::null;
        if (blkSize == 1) {
          coordImporter = rebalanceImporter;
        } else {
          // NOTE: there is an implicit assumption here: we assume that dof any node are enumerated consequently
          // Proper fix would require using decomposition similar to how we construct importer in the
          // RepartitionFactory
          RCP<const Map> origMap = coords->getMap();
          GO indexBase           = origMap->getIndexBase();

          ArrayView<const GO> OEntries = rebalanceImporter->getTargetMap()->getLocalElementList();
          LO numEntries                = OEntries.size() / blkSize;
          ArrayRCP<GO> Entries(numEntries);
          for (LO i = 0; i < numEntries; i++)
            Entries[i] = (OEntries[i * blkSize] - indexBase) / blkSize + indexBase;

          RCP<const Map> targetMap = MapFactory::Build(origMap->lib(), origMap->getGlobalNumElements(), Entries(), indexBase, origMap->getComm());
          coordImporter            = ImportFactory::Build(origMap, targetMap);
        }

        RCP<xdMV> permutedCoords = Xpetra::MultiVectorFactory<typename Teuchos::ScalarTraits<Scalar>::magnitudeType, LO, GO, NO>::Build(coordImporter->getTargetMap(), coords->getNumVectors());
        permutedCoords->doImport(*coords, *coordImporter, Xpetra::INSERT);

        const ParameterList &pL = GetParameterList();
        if (pL.isParameter("repartition: use subcommunicators") == true && pL.get<bool>("repartition: use subcommunicators") == true)
          permutedCoords->replaceMap(permutedCoords->getMap()->removeEmptyProcesses());

        Set(coarseLevel, "Coordinates", permutedCoords);
      }
    }  // end rebalance P(1,1)
    else {
      RCP<ParameterList> params = rcp(new ParameterList());
      params->set("printLoadBalancingInfo", true);
      std::stringstream ss;
      ss << "P(" << curBlockId << "," << curBlockId << ") not rebalanced:";
      GetOStream(Statistics0) << PerfUtils::PrintMatrixInfo(*Pii, ss.str(), params);
      // store rebalanced P block
      subBlockRebP.push_back(Pii);

      // Store Coordinates on coarse level (generated by this)
      // TAW: Note, that each sub-block manager overwrites the Coordinates. So far we only support one set of Coordinates
      //      for a multiphysics problem (i.e., we only support volume coupled problems with the same mesh)
      if ((*it)->hasFactory("Coordinates") == true && coarseLevel.IsAvailable("Coordinates", (*it)->GetFactory("Coordinates").get()) == true) {
        coarseLevel.Set("Coordinates", coarseLevel.Get<RCP<xdMV> >("Coordinates", (*it)->GetFactory("Coordinates").get()), this);
      }
    }

    // fix striding information for rebalanced diagonal block Pii
    // RCP<const Xpetra::MapExtractor<Scalar, LocalOrdinal, GlobalOrdinal, Node> > rgPMapExtractor = bOriginalTransferOp->getRangeMapExtractor(); // original map extractor
    Teuchos::RCP<const StridedMap> orig_stridedRgMap = Teuchos::rcp_dynamic_cast<const StridedMap>(rangeMapExtractor->getMap(Teuchos::as<size_t>(curBlockId), rangeMapExtractor->getThyraMode()));
    Teuchos::RCP<const Map> stridedRgMap             = Teuchos::null;
    if (orig_stridedRgMap != Teuchos::null) {
      std::vector<size_t> stridingData                       = orig_stridedRgMap->getStridingData();
      Teuchos::ArrayView<const GlobalOrdinal> nodeRangeMapii = Pii->getRangeMap()->getLocalElementList();
      stridedRgMap                                           = StridedMapFactory::Build(
                                                    originalTransferOp->getRangeMap()->lib(),
                                                    Teuchos::OrdinalTraits<Xpetra::global_size_t>::invalid(),
                                                    nodeRangeMapii,
                                                    Pii->getRangeMap()->getIndexBase(),
                                                    stridingData,
                                                    originalTransferOp->getRangeMap()->getComm(),
                                                    orig_stridedRgMap->getStridedBlockId(),
                                                    orig_stridedRgMap->getOffset());
    } else
      stridedRgMap = Pii->getRangeMap();

    // RCP<const Xpetra::MapExtractor<Scalar, LocalOrdinal, GlobalOrdinal, Node> > doPMapExtractor = bOriginalTransferOp->getDomainMapExtractor(); // original map extractor
    Teuchos::RCP<const StridedMap> orig_stridedDoMap = Teuchos::rcp_dynamic_cast<const StridedMap>(domainMapExtractor->getMap(Teuchos::as<size_t>(curBlockId), domainMapExtractor->getThyraMode()));

    Teuchos::RCP<const Map> stridedDoMap = Teuchos::null;
    if (orig_stridedDoMap != Teuchos::null) {
      std::vector<size_t> stridingData                        = orig_stridedDoMap->getStridingData();
      Teuchos::ArrayView<const GlobalOrdinal> nodeDomainMapii = Pii->getDomainMap()->getLocalElementList();
      stridedDoMap                                            = StridedMapFactory::Build(originalTransferOp->getDomainMap()->lib(),
                                                                                         Teuchos::OrdinalTraits<Xpetra::global_size_t>::invalid(),
                                                                                         nodeDomainMapii,
                                                                                         Pii->getDomainMap()->getIndexBase(),
                                                                                         stridingData,
                                                                                         originalTransferOp->getDomainMap()->getComm(),
                                                                                         orig_stridedDoMap->getStridedBlockId(),
                                                                                         orig_stridedDoMap->getOffset());
    } else
      stridedDoMap = Pii->getDomainMap();

    TEUCHOS_TEST_FOR_EXCEPTION(stridedRgMap == Teuchos::null, Exceptions::RuntimeError, "MueLu::RebalanceBlockInterpolationFactory::Build: failed to generate striding information. error.");
    TEUCHOS_TEST_FOR_EXCEPTION(stridedDoMap == Teuchos::null, Exceptions::RuntimeError, "MueLu::RebalanceBlockInterpolationFactory::Build: failed to generate striding information. error.");

    // replace stridedMaps view in diagonal sub block
    if (Pii->IsView("stridedMaps")) Pii->RemoveView("stridedMaps");
    Pii->CreateView("stridedMaps", stridedRgMap, stridedDoMap);

    // append strided row map (= range map) to list of range maps.
    Teuchos::RCP<const Map> rangeMapii = Pii->getRowMap("stridedMaps");  // strided range map
    subBlockPRangeMaps.push_back(rangeMapii);
    Teuchos::ArrayView<const GlobalOrdinal> nodeRangeMapii = Pii->getRangeMap()->getLocalElementList();
    // append the GIDs in the end. Do not sort if we have Thyra style GIDs
    fullRangeMapVector.insert(fullRangeMapVector.end(), nodeRangeMapii.begin(), nodeRangeMapii.end());
    if (bThyraRangeGIDs == false)
      sort(fullRangeMapVector.begin(), fullRangeMapVector.end());

    // append strided col map (= domain map) to list of range maps.
    Teuchos::RCP<const Map> domainMapii = Pii->getColMap("stridedMaps");  // strided domain map
    subBlockPDomainMaps.push_back(domainMapii);
    Teuchos::ArrayView<const GlobalOrdinal> nodeDomainMapii = Pii->getDomainMap()->getLocalElementList();
    // append the GIDs in the end. Do not sort if we have Thyra style GIDs
    fullDomainMapVector.insert(fullDomainMapVector.end(), nodeDomainMapii.begin(), nodeDomainMapii.end());
    if (bThyraDomainGIDs == false)
      sort(fullDomainMapVector.begin(), fullDomainMapVector.end());

    curBlockId++;  // increase block id index

  }  // end SubFactoryManager environment

  // extract map index base from maps of blocked P
  GO rangeIndexBase  = originalTransferOp->getRangeMap()->getIndexBase();
  GO domainIndexBase = originalTransferOp->getDomainMap()->getIndexBase();

  RCP<const Xpetra::MapExtractor<Scalar, LocalOrdinal, GlobalOrdinal, Node> > rangePMapExtractor = bOriginalTransferOp->getRangeMapExtractor();  // original map extractor
  Teuchos::ArrayView<GO> fullRangeMapGIDs(fullRangeMapVector.size() ? &fullRangeMapVector[0] : 0, fullRangeMapVector.size());
  Teuchos::RCP<const StridedMap> stridedRgFullMap = Teuchos::rcp_dynamic_cast<const StridedMap>(rangePMapExtractor->getFullMap());
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

  RCP<const Xpetra::MapExtractor<Scalar, LocalOrdinal, GlobalOrdinal, Node> > domainAMapExtractor = bOriginalTransferOp->getDomainMapExtractor();
  Teuchos::ArrayView<GO> fullDomainMapGIDs(fullDomainMapVector.size() ? &fullDomainMapVector[0] : 0, fullDomainMapVector.size());
  Teuchos::RCP<const StridedMap> stridedDoFullMap = Teuchos::rcp_dynamic_cast<const StridedMap>(domainAMapExtractor->getFullMap());
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

  // build map extractors
  Teuchos::RCP<const MapExtractor> rebrangeMapExtractor =
      MapExtractorFactory::Build(fullRangeMap, subBlockPRangeMaps, bThyraRangeGIDs);
  Teuchos::RCP<const MapExtractor> rebdomainMapExtractor =
      MapExtractorFactory::Build(fullDomainMap, subBlockPDomainMaps, bThyraDomainGIDs);

  Teuchos::RCP<BlockedCrsMatrix> bRebP = Teuchos::rcp(new BlockedCrsMatrix(rebrangeMapExtractor, rebdomainMapExtractor, 10));
  for (size_t i = 0; i < subBlockPRangeMaps.size(); i++) {
    Teuchos::RCP<CrsMatrixWrap> crsOpii = Teuchos::rcp_dynamic_cast<CrsMatrixWrap>(subBlockRebP[i]);
    TEUCHOS_TEST_FOR_EXCEPTION(crsOpii == Teuchos::null, Xpetra::Exceptions::BadCast, "MueLu::RebalanceBlockTransferFactory::Build: block P" << i << " is not of type CrsMatrixWrap.");
    bRebP->setMatrix(i, i, crsOpii);
  }
  bRebP->fillComplete();

  Set(coarseLevel, "P", Teuchos::rcp_dynamic_cast<Matrix>(bRebP));

  // Finish up the coordinates (single source only)
  if (UseSingleSource) {
    RCP<xdBV> bcoarseCoords = rcp(new xdBV(rebrangeMapExtractor->getBlockedMap(), newCoordinates));
    Set(coarseLevel, "Coordinates", bcoarseCoords);
  }

}  // Build

}  // namespace MueLu

#endif /* MUELU_REBALANCEBLOCKINTERPOLATIONFACTORY_DEF_HPP_ */
