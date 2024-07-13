// @HEADER
// *****************************************************************************
//        MueLu: A package for multigrid based preconditioning
//
// Copyright 2012 NTESS and the MueLu contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef MUELU_SETUPREGIONHIERARCHY_DEF_HPP
#define MUELU_SETUPREGIONHIERARCHY_DEF_HPP

#include <vector>
#include <iostream>

#include <Tpetra_KokkosCompat_DefaultNode.hpp>

#include <Teuchos_RCP.hpp>

#include <Xpetra_ConfigDefs.hpp>
#include <Xpetra_Export.hpp>
#include <Xpetra_Import.hpp>
#include <Xpetra_Map.hpp>
#include <Xpetra_MultiVector.hpp>
#include <Xpetra_Vector.hpp>
#include <Xpetra_CrsMatrixWrap.hpp>

#include <MueLu_CreateXpetraPreconditioner.hpp>
#include <MueLu_Utilities.hpp>

#if defined(HAVE_MUELU_ZOLTAN2) && defined(HAVE_MPI)
#include <MueLu_RepartitionFactory.hpp>
#include <MueLu_RepartitionHeuristicFactory.hpp>
#include <MueLu_Zoltan2Interface.hpp>
#endif

#include "SetupRegionVector_def.hpp"
#include "SetupRegionMatrix_def.hpp"
#include "SetupRegionSmoothers_def.hpp"

#if defined(HAVE_MUELU_AMESOS2)
#include <Amesos2_config.h>
#include <Amesos2.hpp>
#endif

using Teuchos::Array;
using Teuchos::ArrayRCP;
using Teuchos::ArrayView;
using Teuchos::ParameterList;
using Teuchos::RCP;

/*! \brief Create coarse level maps with continuous GIDs
 *
 *  The direct solver requires maps with continuous GIDs. Starting from the
 *  coarse level composite maps with discontinuous GIDs, we create a new row map
 *  and a matching column map.
 *
 *  Range and Domain map happen to correspond to the Row map, so we don't have
 *  to deal with them in particular.
 */
template <class LocalOrdinal, class GlobalOrdinal, class Node>
void createContinuousCoarseLevelMaps(const RCP<const Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node>> rowMap,  ///< row map
                                     const RCP<const Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node>> colMap,  ///< column map
                                     RCP<Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node>>& contRowMap,         ///< row map with continuous GIDs
                                     RCP<Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node>>& contColMap          ///< column map with continuous GIDs
) {
#include "Xpetra_UseShortNamesOrdinal.hpp"
  //   /!\ This function is pure ordinal, no scalar type is passed as input
  //       This means that use only three template paramters and that we
  //       do not use the Scalar dependent short names!

  // Create row map with continuous GIDs
  contRowMap = MapFactory::Build(rowMap->lib(),
                                 rowMap->getGlobalNumElements(),
                                 rowMap->getLocalNumElements(),
                                 rowMap->getIndexBase(),
                                 rowMap->getComm());

  return;
}  // createContinuousCoarseLevelMaps

/* Reconstruct coarse-level maps (assuming fully structured grids)
 *
 * We know the regional map on the coarse levels since they are just the
 * row maps of the coarse level operators. Though, we need to re-construct
 * the quasiRegional and composite maps ourselves.
 *
 * We ultimately are only interested in the composite map on the coarsest level.
 * Intermediate levels are dealt with along the way, because we go through all
 * levels recursively.
 *
 * Assumptions:
 * - fully structured grid
 * - only on region per proc and group
 */
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void MakeCoarseLevelMaps(const int maxRegPerGID,
                         Teuchos::RCP<MueLu::Hierarchy<Scalar, LocalOrdinal, GlobalOrdinal, Node>> regHierarchy) {
#include "Xpetra_UseShortNames.hpp"
#include "MueLu_UseShortNames.hpp"

  using MT = typename Teuchos::ScalarTraits<SC>::magnitudeType;

  RCP<Level> level0 = regHierarchy->GetLevel(0);

  const GO GO_INV     = Teuchos::OrdinalTraits<GO>::invalid();
  const int numLevels = regHierarchy->GetNumLevels();

  // RCP<Teuchos::FancyOStream> fancy = Teuchos::fancyOStream(Teuchos::rcpFromRef(std::cout));
  // Teuchos::FancyOStream& out = *fancy;
  const int myRank = level0->GetComm()->getRank();

  Teuchos::Array<LO> coarseCompositeToRegionLIDs;
  Teuchos::ArrayView<LO> compositeToRegionLIDs = level0->Get<ArrayView<LO>>("compositeToRegionLIDs");

  for (int currentLevel = 1; currentLevel < numLevels; ++currentLevel) {
    RCP<Level> level         = regHierarchy->GetLevel(currentLevel);
    RCP<Matrix> regProlong   = level->Get<RCP<Matrix>>("P");
    RCP<const Map> regRowMap = regProlong->getColMap();

    RCP<Level> levelFine         = regHierarchy->GetLevel(currentLevel - 1);
    RCP<Import> regRowImportFine = levelFine->Get<RCP<Import>>("rowImport");
    RCP<Matrix> regMatFine       = levelFine->Get<RCP<Matrix>>("A");
    RCP<const Map> regRowMapFine = regMatFine->getRowMap();

    // Extracting some basic information about local mesh in composite/region format
    const size_t numFineRegionNodes    = regProlong->getLocalNumRows();
    const size_t numFineCompositeNodes = compositeToRegionLIDs.size();
    const size_t numFineDuplicateNodes = numFineRegionNodes - numFineCompositeNodes;

    const size_t numCoarseRegionNodes = regProlong->getColMap()->getLocalNumElements();

    // Find the regionLIDs associated with local duplicated nodes
    // This will allow us to later loop only on duplicated nodes
    size_t countComposites = 0, countDuplicates = 0;
    Array<LO> fineDuplicateLIDs(numFineDuplicateNodes);
    for (size_t regionIdx = 0; regionIdx < numFineRegionNodes; ++regionIdx) {
      if (compositeToRegionLIDs[countComposites] == static_cast<LO>(regionIdx)) {
        ++countComposites;
      } else {
        fineDuplicateLIDs[countDuplicates] = regionIdx;
        ++countDuplicates;
      }
    }

    // We gather the coarse GIDs associated with each fine point in the local composite mesh part.
    RCP<Xpetra::Vector<MT, LO, GO, NO>> coarseCompositeGIDs = Xpetra::VectorFactory<MT, LO, GO, NO>::Build(regRowImportFine->getSourceMap(), false);
    Teuchos::ArrayRCP<MT> coarseCompositeGIDsData           = coarseCompositeGIDs->getDataNonConst(0);

    for (size_t compositeNodeIdx = 0; compositeNodeIdx < numFineCompositeNodes; ++compositeNodeIdx) {
      ArrayView<const LO> coarseRegionLID;  // Should contain a single value
      ArrayView<const SC> dummyData;        // Should contain a single value
      regProlong->getLocalRowView(compositeToRegionLIDs[compositeNodeIdx],
                                  coarseRegionLID,
                                  dummyData);
      if (coarseRegionLID.size() == 1) {
        coarseCompositeGIDsData[compositeNodeIdx] = regProlong->getColMap()->getGlobalElement(coarseRegionLID[0]);
      } else {
        coarseCompositeGIDsData[compositeNodeIdx] = -1;
      }
    }

    // We communicate the above GIDs to their duplicate so that we can replace GIDs of the region
    // column map and form the quasiRegion column map.
    RCP<Xpetra::Vector<MT, LO, GO, NO>> coarseQuasiregionGIDs;
    RCP<Xpetra::Vector<MT, LO, GO, NO>> coarseRegionGIDs;
    compositeToRegional(coarseCompositeGIDs,
                        coarseQuasiregionGIDs,
                        coarseRegionGIDs,
                        regRowMapFine,
                        regRowImportFine);

    RCP<Xpetra::MultiVector<LocalOrdinal, LocalOrdinal, GlobalOrdinal, Node>> regionsPerGIDWithGhosts =
        Xpetra::MultiVectorFactory<LO, LO, GO, NO>::Build(regRowMapFine,
                                                          maxRegPerGID,
                                                          false);
    RCP<Xpetra::MultiVector<GO, LO, GO, NO>> interfaceGIDs = Xpetra::MultiVectorFactory<GO, LO, GO, NO>::Build(regRowMapFine,
                                                                                                               maxRegPerGID,
                                                                                                               false);

    Array<ArrayRCP<const LO>> regionsPerGIDWithGhostsFine(maxRegPerGID);
    Array<ArrayRCP<LO>> regionsPerGIDWithGhostsCoarse(maxRegPerGID);
    Array<ArrayRCP<GO>> interfaceGIDsCoarse(maxRegPerGID);
    for (size_t idx = 0; idx < static_cast<size_t>(maxRegPerGID); ++idx) {
      regionsPerGIDWithGhostsFine[idx]   = levelFine->Get<RCP<Xpetra::MultiVector<LO, LO, GO, NO>>>("regionsPerGIDWithGhosts")->getData(idx);
      regionsPerGIDWithGhostsCoarse[idx] = regionsPerGIDWithGhosts->getDataNonConst(idx);
      interfaceGIDsCoarse[idx]           = interfaceGIDs->getDataNonConst(idx);
      for (size_t coarseIdx = 0;
           coarseIdx < regionsPerGIDWithGhosts->getLocalLength(); ++coarseIdx) {
        regionsPerGIDWithGhostsCoarse[idx][coarseIdx] = -1;
        interfaceGIDsCoarse[idx][coarseIdx]           = 0;
      }
    }

    for (size_t fineIdx = 0; fineIdx < numFineRegionNodes; ++fineIdx) {
      ArrayView<const LO> coarseRegionLID;  // Should contain a single value
      ArrayView<const SC> dummyData;        // Should contain a single value
      regProlong->getLocalRowView(fineIdx,
                                  coarseRegionLID,
                                  dummyData);
      const LO coarseIdx = coarseRegionLID[0];

      // Now fill regionPerGIDWithGhostsCoarse[:][coarseRegionLID]
      // with data from regionPerGIDWithGhostsFine[:][fineIdx].
      // The problem is we might have more then maxRegPerGID on currentLevel
      // then on currentLevel-1... we might need to do a union or something?
      // I guess technically using the restriction operator here would be more
      // helpful than the prolongator, this way we could find all the fine interface
      // points easily and compute the union of the PIDs they belong too.
      // This can actually be done extracting the local matrix and compute the transpose.
      // For now let us assume that maxRegPerGID is constant and hope for the best.
      LO countFinePIDs   = 0;
      LO countCoarsePIDs = 0;
      for (LO idx = 0; idx < maxRegPerGID; ++idx) {
        if (-1 < regionsPerGIDWithGhostsFine[idx][fineIdx]) {
          ++countFinePIDs;
        }
        if (-1 < regionsPerGIDWithGhostsCoarse[idx][coarseIdx]) {
          ++countCoarsePIDs;
        }
      }

      if (countCoarsePIDs < countFinePIDs) {
        for (LO idx = 0; idx < countFinePIDs; ++idx) {
          regionsPerGIDWithGhostsCoarse[idx][coarseIdx] = regionsPerGIDWithGhostsFine[idx][fineIdx];
          if (regionsPerGIDWithGhostsCoarse[idx][coarseIdx] == myRank) {
            interfaceGIDsCoarse[idx][coarseIdx] = regRowMap->getGlobalElement(coarseIdx);
          }
        }
      }
    }

    Array<GO> fineRegionDuplicateCoarseLIDs(numFineDuplicateNodes);
    Array<GO> fineRegionDuplicateCoarseGIDs(numFineDuplicateNodes);
    for (size_t duplicateIdx = 0; duplicateIdx < numFineDuplicateNodes; ++duplicateIdx) {
      ArrayView<const LO> coarseRegionLID;  // Should contain a single value
      ArrayView<const SC> dummyData;        // Should contain a single value
      regProlong->getLocalRowView(fineDuplicateLIDs[duplicateIdx],
                                  coarseRegionLID,
                                  dummyData);
      fineRegionDuplicateCoarseLIDs[duplicateIdx] = regProlong->getColMap()->getGlobalElement(coarseRegionLID[0]);
      fineRegionDuplicateCoarseGIDs[duplicateIdx] = (coarseQuasiregionGIDs->getDataNonConst(0))[fineDuplicateLIDs[duplicateIdx]];
    }

    // Create the coarseQuasiregRowMap, it will be based on the coarseRegRowMap
    LO countCoarseComposites = 0;
    coarseCompositeToRegionLIDs.resize(numCoarseRegionNodes);
    Array<GO> coarseQuasiregRowMapData = regProlong->getColMap()->getLocalElementList();
    Array<GO> coarseCompRowMapData(numCoarseRegionNodes, -1);
    for (size_t regionIdx = 0; regionIdx < numCoarseRegionNodes; ++regionIdx) {
      const GO initialValue = coarseQuasiregRowMapData[regionIdx];
      for (size_t duplicateIdx = 0; duplicateIdx < numFineDuplicateNodes; ++duplicateIdx) {
        if ((initialValue == fineRegionDuplicateCoarseLIDs[duplicateIdx]) &&
            (fineRegionDuplicateCoarseGIDs[duplicateIdx] < coarseQuasiregRowMapData[regionIdx]) &&
            (-1 < fineRegionDuplicateCoarseGIDs[duplicateIdx])) {
          coarseQuasiregRowMapData[regionIdx] = fineRegionDuplicateCoarseGIDs[duplicateIdx];
        }
      }
      if (initialValue == coarseQuasiregRowMapData[regionIdx]) {
        coarseCompRowMapData[countCoarseComposites]        = coarseQuasiregRowMapData[regionIdx];
        coarseCompositeToRegionLIDs[countCoarseComposites] = regionIdx;
        ++countCoarseComposites;
      }
    }
    coarseCompRowMapData.resize(countCoarseComposites);
    coarseCompositeToRegionLIDs.resize(countCoarseComposites);

    // We are now ready to fill up the outputs
    RCP<const Map> regRowMapCurrent = regProlong->getColMap();

    RCP<Map> quasiRegRowMap = MapFactory::Build(regProlong->getColMap()->lib(),
                                                GO_INV,
                                                coarseQuasiregRowMapData(),
                                                regProlong->getColMap()->getIndexBase(),
                                                regProlong->getColMap()->getComm());

    RCP<Map> compRowMap = MapFactory::Build(regProlong->getColMap()->lib(),
                                            GO_INV,
                                            coarseCompRowMapData(),
                                            regProlong->getColMap()->getIndexBase(),
                                            regProlong->getColMap()->getComm());

    RCP<Import> regRowImportCurrent = ImportFactory::Build(compRowMap, quasiRegRowMap);

    // Now generate matvec data
    Teuchos::ArrayRCP<LocalOrdinal> regionMatVecLIDs;
    Teuchos::RCP<Xpetra::Import<LocalOrdinal, GlobalOrdinal, Node>> regionInterfaceImporter;
    SetupMatVec(interfaceGIDs, regionsPerGIDWithGhosts,
                regRowMapCurrent, regRowImportCurrent,
                regionMatVecLIDs, regionInterfaceImporter);

    // Fill level with the outputs
    level->Set<Teuchos::ArrayRCP<LO>>("regionMatVecLIDs", regionMatVecLIDs);
    level->Set<RCP<Xpetra::Import<LocalOrdinal, GlobalOrdinal, Node>>>("regionInterfaceImporter", regionInterfaceImporter);
    level->Set<RCP<Xpetra::MultiVector<GlobalOrdinal, LocalOrdinal, GlobalOrdinal, Node>>>("interfaceGIDs", interfaceGIDs);
    level->Set<RCP<Xpetra::MultiVector<LO, LO, GO, NO>>>("regionsPerGIDWithGhosts", regionsPerGIDWithGhosts);
    level->Set<RCP<Xpetra::Import<LocalOrdinal, GlobalOrdinal, Node>>>("rowImport", regRowImportCurrent);

    // Finally reset compositeToRegionLIDs
    compositeToRegionLIDs = coarseCompositeToRegionLIDs();
  }  // Loop over numLevels
}  // MakeCoarseLevelMaps

// Form the composite coarse level operator
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void MakeCoarseCompositeOperator(RCP<const Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node>>& compRowMap,
                                 RCP<const Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node>>& quasiRegRowMap,
                                 RCP<const Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node>>& quasiRegColMap,
                                 RCP<Xpetra::Import<LocalOrdinal, GlobalOrdinal, Node>>& regRowImporter,
                                 RCP<Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>>& regMatrix,
                                 RCP<Xpetra::MultiVector<typename Teuchos::ScalarTraits<Scalar>::coordinateType, LocalOrdinal, GlobalOrdinal, Node>>& regCoarseCoordinates,
                                 RCP<Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>>& coarseCompOp,
                                 RCP<Xpetra::MultiVector<typename Teuchos::ScalarTraits<Scalar>::coordinateType, LocalOrdinal, GlobalOrdinal, Node>>& compCoarseCoordinates,
                                 const bool makeCompCoords) {
#include "Xpetra_UseShortNames.hpp"
  using CoordType = typename Teuchos::ScalarTraits<Scalar>::coordinateType;
  coarseCompOp    = MatrixFactory::Build(compRowMap,
                                         // This estimate is very conservative and probably costs us lots of memory...
                                         8 * regMatrix->getCrsGraph()->getLocalMaxNumRowEntries());
  regionalToComposite(regMatrix,
                      quasiRegRowMap,
                      quasiRegColMap,
                      regRowImporter, Xpetra::ADD,
                      coarseCompOp);

  const int dofsPerNode = regMatrix->GetFixedBlockSize();
  coarseCompOp->SetFixedBlockSize(dofsPerNode);
  coarseCompOp->setObjectLabel("coarse composite operator");

  // Create coarse composite coordinates for repartitioning
  if (makeCompCoords) {
    const int check = regMatrix->getRowMap()->getLocalNumElements() % regCoarseCoordinates->getMap()->getLocalNumElements();
    TEUCHOS_ASSERT(check == 0);

    RCP<const Map> compCoordMap;
    RCP<Xpetra::Import<LocalOrdinal, GlobalOrdinal, Node>> regCoordImporter;
    if (dofsPerNode == 1) {
      compCoordMap     = compRowMap;
      regCoordImporter = regRowImporter;
    } else {
      using size_type = typename Teuchos::Array<GlobalOrdinal>::size_type;
      Array<GlobalOrdinal> compCoordMapData(compRowMap->getLocalNumElements() / dofsPerNode);
      ArrayView<const GlobalOrdinal> compRowMapData = compRowMap->getLocalElementList();
      for (size_type nodeIdx = 0; nodeIdx < compCoordMapData.size(); ++nodeIdx) {
        compCoordMapData[nodeIdx] = compRowMapData[nodeIdx * dofsPerNode] / dofsPerNode;
      }
      compCoordMap = MapFactory::Build(compRowMap->lib(),
                                       compRowMap->getGlobalNumElements() / dofsPerNode,
                                       compCoordMapData(),
                                       compRowMap->getIndexBase(),
                                       compRowMap->getComm());

      RCP<Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node>> quasiRegCoordMap;
      Array<GlobalOrdinal> quasiRegCoordMapData(quasiRegRowMap->getLocalNumElements() / dofsPerNode);
      ArrayView<const GlobalOrdinal> quasiRegRowMapData = quasiRegRowMap->getLocalElementList();
      for (size_type nodeIdx = 0; nodeIdx < quasiRegCoordMapData.size(); ++nodeIdx) {
        quasiRegCoordMapData[nodeIdx] = quasiRegRowMapData[nodeIdx * dofsPerNode] / dofsPerNode;
      }
      quasiRegCoordMap = MapFactory::Build(quasiRegRowMap->lib(),
                                           quasiRegRowMap->getGlobalNumElements() / dofsPerNode,
                                           quasiRegCoordMapData(),
                                           quasiRegRowMap->getIndexBase(),
                                           quasiRegRowMap->getComm());
      regCoordImporter = ImportFactory::Build(compCoordMap, quasiRegCoordMap);
    }
    compCoarseCoordinates = Xpetra::MultiVectorFactory<CoordType, LocalOrdinal, GlobalOrdinal, Node>::Build(compCoordMap, regCoarseCoordinates->getNumVectors());
    TEUCHOS_ASSERT(Teuchos::nonnull(compCoarseCoordinates));

    // The following looks like regionalToComposite for Vector
    // but it is a bit different since we do not want to add
    // entries in the coordinate MultiVector as we would for
    // a solution or residual vector.
    RCP<Xpetra::MultiVector<CoordType, LocalOrdinal, GlobalOrdinal, Node>> quasiRegCoarseCoordinates;
    quasiRegCoarseCoordinates = regCoarseCoordinates;
    TEUCHOS_ASSERT(Teuchos::nonnull(quasiRegCoarseCoordinates));
    quasiRegCoarseCoordinates->replaceMap(regCoordImporter->getTargetMap());
    compCoarseCoordinates->doExport(*quasiRegCoarseCoordinates, *(regCoordImporter), Xpetra::INSERT);
  }
}  // MakeCoarseCompositeOperator

/* Create a direct solver for a composite operator
 *
 * Create the solver object and compute symbolic and numeric factorization.
 * Finally, the solver object will be ready to be applied during the V-cycle call.
 *
 * \note For now, we're limited to Tpetra/Amesos2. From Amesos2, we use KLU as direct solver.
 */
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
RCP<Amesos2::Solver<Tpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>, Tpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>>>
MakeCompositeDirectSolver(RCP<Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>>& compOp) {
  using Tpetra_CrsMatrix   = Tpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>;
  using Tpetra_MultiVector = Tpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>;
  using Utilities          = MueLu::Utilities<Scalar, LocalOrdinal, GlobalOrdinal, Node>;
  using Teuchos::TimeMonitor;

  RCP<Amesos2::Solver<Tpetra_CrsMatrix, Tpetra_MultiVector>> coarseSolver;
  {
    RCP<TimeMonitor> tm = rcp(new TimeMonitor(*TimeMonitor::getNewTimer("MakeCompositeDirectSolver: 1 - Setup")));

    // convert matrix to Tpetra
    RCP<Tpetra_CrsMatrix> tMat = Utilities::Op2NonConstTpetraCrs(compOp);

    // Amesos2-specific key phrase that denote smoother type
    std::string amesos2SolverName = "KLU2";
    TEUCHOS_ASSERT(Amesos2::query(amesos2SolverName));
    coarseSolver = Amesos2::create<Tpetra_CrsMatrix, Tpetra_MultiVector>(amesos2SolverName, tMat);

    Teuchos::ParameterList amesos2_params("Amesos2");
    amesos2_params.sublist(amesos2SolverName).set("IsContiguous", false, "Are GIDs Contiguous");
    coarseSolver->setParameters(Teuchos::rcpFromRef(amesos2_params));
  }

  {
    RCP<TimeMonitor> tm = rcp(new TimeMonitor(*TimeMonitor::getNewTimer("MakeCompositeDirectSolver: 2 - Factorization")));

    coarseSolver->symbolicFactorization();
    coarseSolver->numericFactorization();
  }

  return coarseSolver;
}  // MakeCorseCompositeDirectSolver

/* Rebalance coarse operator
 *
 */
#if defined(HAVE_MUELU_ZOLTAN2) && defined(HAVE_MPI)
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void RebalanceCoarseCompositeOperator(const int rebalanceNumPartitions,
                                      RCP<Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>>& coarseCompOp,
                                      RCP<Xpetra::MultiVector<typename Teuchos::ScalarTraits<Scalar>::coordinateType, LocalOrdinal, GlobalOrdinal, Node>>& compCoarseCoordinates,
                                      RCP<Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>>& rebalancedCompOp,
                                      RCP<Xpetra::MultiVector<typename Teuchos::ScalarTraits<Scalar>::coordinateType, LocalOrdinal, GlobalOrdinal, Node>>& rebalancedCoordinates,
                                      RCP<const Xpetra::Import<LocalOrdinal, GlobalOrdinal, Node>>& rebalanceImporter) {
#include "MueLu_UseShortNames.hpp"
  using CoordType = typename Teuchos::ScalarTraits<Scalar>::coordinateType;
  using Teuchos::TimeMonitor;

  RCP<TimeMonitor> tm = rcp(new TimeMonitor(*TimeMonitor::getNewTimer("RebalanceCoarseCompositeOperator: ")));

  RCP<Teuchos::FancyOStream> fos = Teuchos::fancyOStream(Teuchos::rcpFromRef(std::cout));
  fos->setOutputToRootOnly(0);
  *fos << "Rebalancing coarse composite operator." << std::endl;

  const int numPartitions = rebalanceNumPartitions;

  // We build a fake level
  Level level;
  level.SetLevelID(1);

  RCP<FactoryManagerBase> factoryHandler = rcp(new FactoryManager());
  level.SetFactoryManager(factoryHandler);

  RCP<TimeMonitor> tmLocal = rcp(new TimeMonitor(*TimeMonitor::getNewTimer("rebalanceCoarse: Zoltan construction:")));

  RCP<Zoltan2Interface> zoltan = rcp(new Zoltan2Interface());

  level.Set<RCP<Matrix>>("A", coarseCompOp);
  level.Set<RCP<MultiVector>>("Coordinates", compCoarseCoordinates);
  //  int numPartitions = Get<int>(level, "number of partitions");

  RCP<RepartitionFactory> repart = rcp(new RepartitionFactory());
  Teuchos::ParameterList paramList;
  paramList.set("repartition: remap parts", false);
  if (numPartitions > 0) {  // If number of coarse rebalance partitions was provided by the user.
    level.Set<int>("number of partitions", numPartitions);
  } else {
    Teuchos::ParameterList paramListHeuristic;
    paramListHeuristic.set("repartition: start level", 1);
    RCP<RepartitionHeuristicFactory> repartHeuristic = rcp(new RepartitionHeuristicFactory());
    repartHeuristic->SetParameterList(paramListHeuristic);
    repart->SetFactory("number of partitions", repartHeuristic);
  }
  repart->SetParameterList(paramList);
  repart->SetFactory("Partition", zoltan);

  // Build
  level.Request("Importer", repart.get());
  repart->Build(level);

  tmLocal = Teuchos::null;

  // Build importer for rebalancing
  level.Get("Importer", rebalanceImporter, repart.get());

  ParameterList XpetraList;
  XpetraList.set("Restrict Communicator", true);
  XpetraList.set("Timer Label", "MueLu::RebalanceAc-for-coarseAMG");

  // Build rebalanced coarse composite operator
  rebalancedCompOp = MatrixFactory::Build(coarseCompOp, *rebalanceImporter, *rebalanceImporter, rebalanceImporter->getTargetMap(), rebalanceImporter->getTargetMap(), rcp(&XpetraList, false));
  if (!rebalancedCompOp.is_null()) {
    rebalancedCompOp->SetFixedBlockSize(coarseCompOp->GetFixedBlockSize());
  }

  // Build rebalanced coarse coordinates (The following code is borrowed from MueLu_RebalanceTransferFactory_def.hpp)
  LO blkSize = coarseCompOp->GetFixedBlockSize();
  RCP<const Import> coordImporter;
  if (blkSize == 1) {
    coordImporter = rebalanceImporter;

  } else {
    // NOTE: there is an implicit assumption here: we assume that dof any node are enumerated consequently
    // Proper fix would require using decomposition similar to how we construct importer in the
    // RepartitionFactory
    RCP<const Map> origMap = compCoarseCoordinates->getMap();
    GO indexBase           = origMap->getIndexBase();

    ArrayView<const GO> OEntries = rebalanceImporter->getTargetMap()->getLocalElementList();
    LO numEntries                = OEntries.size() / blkSize;
    ArrayRCP<GO> Entries(numEntries);
    for (LO i = 0; i < numEntries; i++)
      Entries[i] = (OEntries[i * blkSize] - indexBase) / blkSize + indexBase;

    RCP<const Map> targetMap = MapFactory::Build(origMap->lib(), origMap->getGlobalNumElements(), Entries(), indexBase, origMap->getComm());
    coordImporter            = ImportFactory::Build(origMap, targetMap);
  }
  rebalancedCoordinates = Xpetra::MultiVectorFactory<CoordType, LocalOrdinal, GlobalOrdinal, Node>::Build(coordImporter->getTargetMap(), compCoarseCoordinates->getNumVectors());
  rebalancedCoordinates->doImport(*compCoarseCoordinates, *coordImporter, Xpetra::INSERT);
  rebalancedCoordinates->replaceMap(rebalancedCoordinates->getMap()->removeEmptyProcesses());

  return;
}  // RebalanceCoarseCompositeOperator
#endif

/* Create an AMG hierarchy for a composite operator
 *
 * Create the hierarchy object and perform the multigrid setup.
 * Finally, the hierarhcy object will be ready to be applied during the region MG V-cycle call.
 */
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
RCP<MueLu::Hierarchy<Scalar, LocalOrdinal, GlobalOrdinal, Node>>
MakeCompositeAMGHierarchy(RCP<Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>>& compOp,
                          const std::string& xmlFileName,
                          RCP<Xpetra::MultiVector<typename Teuchos::ScalarTraits<Scalar>::coordinateType, LocalOrdinal, GlobalOrdinal, Node>> coordinates) {
#include "MueLu_UseShortNames.hpp"
  using coordinates_type = typename Teuchos::ScalarTraits<Scalar>::coordinateType;

  const Scalar one = Teuchos::ScalarTraits<Scalar>::one();

  RCP<Teuchos::FancyOStream> fos = Teuchos::fancyOStream(Teuchos::rcpFromRef(std::cout));
  fos->setOutputToRootOnly(0);
  *fos << "Attempting to setup AMG hierarchy for the composite coarse grid problem" << std::endl;

  // Get parameter list for AMG hierarchy
  RCP<ParameterList> mueluParams = Teuchos::rcp(new ParameterList());
  Teuchos::updateParametersFromXmlFileAndBroadcast(xmlFileName, mueluParams.ptr(),
                                                   *compOp->getRowMap()->getComm());

  // Get the user data sublist
  const std::string userName            = "user data";
  Teuchos::ParameterList& userParamList = mueluParams->sublist(userName);

  // Add nullspace information
  {
    // Compute nullspace
    RCP<MultiVector> nullspace;
    if ((compOp->GetFixedBlockSize() == 1) || Teuchos::is_null(coordinates)) {  // Scalar problem, constant nullspace
      nullspace = MultiVectorFactory::Build(compOp->getRowMap(), 1);
      nullspace->putScalar(one);
    } else if (compOp->GetFixedBlockSize() == 2) {  // 2D Elasticity
      nullspace = MultiVectorFactory::Build(compOp->getRowMap(), 3);
      Array<ArrayRCP<SC>> nullspaceData(3);
      Array<ArrayRCP<const coordinates_type>> coordinateData(2);

      // Calculate center
      const coordinates_type cx = coordinates->getVector(0)->meanValue();
      const coordinates_type cy = coordinates->getVector(1)->meanValue();

      coordinateData[0] = coordinates->getData(0);
      coordinateData[1] = coordinates->getData(1);

      for (int vecIdx = 0; vecIdx < 3; ++vecIdx) {
        nullspaceData[vecIdx] = nullspace->getDataNonConst(vecIdx);
      }

      for (size_t nodeIdx = 0; nodeIdx < coordinates->getLocalLength(); ++nodeIdx) {
        // translations
        nullspaceData[0][2 * nodeIdx + 0] = one;
        nullspaceData[1][2 * nodeIdx + 1] = one;

        // rotation about z axis
        nullspaceData[2][2 * nodeIdx + 0] = -(coordinateData[1][nodeIdx] - cy);
        nullspaceData[2][2 * nodeIdx + 1] = (coordinateData[0][nodeIdx] - cx);
      }

    } else if (compOp->GetFixedBlockSize() == 3) {  // 3D Elasticity
      nullspace = MultiVectorFactory::Build(compOp->getRowMap(), 6);
      Array<ArrayRCP<SC>> nullspaceData(6);
      Array<ArrayRCP<const coordinates_type>> coordinateData(3);

      // Calculate center
      const coordinates_type cx = coordinates->getVector(0)->meanValue();
      const coordinates_type cy = coordinates->getVector(1)->meanValue();
      const coordinates_type cz = coordinates->getVector(2)->meanValue();

      coordinateData[0] = coordinates->getData(0);
      coordinateData[1] = coordinates->getData(1);
      coordinateData[2] = coordinates->getData(2);

      for (int vecIdx = 0; vecIdx < 6; ++vecIdx) {
        nullspaceData[vecIdx] = nullspace->getDataNonConst(vecIdx);
      }

      for (size_t nodeIdx = 0; nodeIdx < coordinates->getLocalLength(); ++nodeIdx) {
        // translations
        nullspaceData[0][3 * nodeIdx + 0] = one;
        nullspaceData[1][3 * nodeIdx + 1] = one;
        nullspaceData[2][3 * nodeIdx + 2] = one;

        // rotation about z axis
        nullspaceData[3][3 * nodeIdx + 0] = -(coordinateData[1][nodeIdx] - cy);
        nullspaceData[3][3 * nodeIdx + 1] = (coordinateData[0][nodeIdx] - cx);

        // rotation about x axis
        nullspaceData[4][3 * nodeIdx + 1] = -(coordinateData[2][nodeIdx] - cz);
        nullspaceData[4][3 * nodeIdx + 2] = (coordinateData[1][nodeIdx] - cy);

        // rotation about y axis
        nullspaceData[5][3 * nodeIdx + 0] = (coordinateData[2][nodeIdx] - cz);
        nullspaceData[5][3 * nodeIdx + 2] = -(coordinateData[0][nodeIdx] - cx);
      }
    }

    // Equalize norms of all vectors to that of the first one
    // We do not normalize them as a vector of ones seems nice
    Teuchos::Array<typename Teuchos::ScalarTraits<Scalar>::magnitudeType> norms2(nullspace->getNumVectors());
    nullspace->norm2(norms2);
    Teuchos::Array<Scalar> norms2scalar(nullspace->getNumVectors());
    for (size_t vectorIdx = 0; vectorIdx < nullspace->getNumVectors(); ++vectorIdx) {
      norms2scalar[vectorIdx] = norms2[0] / norms2[vectorIdx];
    }
    nullspace->scale(norms2scalar);

    // Insert into parameter list
    userParamList.set("Nullspace", nullspace);
  }

  // Add coordinate information for rebalancing
  {
    if (Teuchos::nonnull(coordinates)) {
      userParamList.set("Coordinates", coordinates);
    }
  }

  // Create an AMG hierarchy based on the composite coarse level operator from the region MG scheme
  RCP<Hierarchy> compOpHiearchy = MueLu::CreateXpetraPreconditioner(compOp, *mueluParams);

  // We will use the hiearchy as a solver
  compOpHiearchy->IsPreconditioner(false);
  compOpHiearchy->SetVerbLevel(MueLu::VERB_NONE);

  return compOpHiearchy;
}  // MakeCompositeAMGHierarchy

// Make interface scaling factors recursively
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void MakeInterfaceScalingFactors(const int numLevels,
                                 Teuchos::RCP<MueLu::Hierarchy<Scalar, LocalOrdinal, GlobalOrdinal, Node>> regHierarchy) {
#include "Xpetra_UseShortNames.hpp"
#include "MueLu_UseShortNames.hpp"

  const SC SC_ONE = Teuchos::ScalarTraits<SC>::one();

  TEUCHOS_TEST_FOR_EXCEPT_MSG(!(numLevels > 0), "We require numLevel > 0. Probably, numLevel has not been set, yet.");

  for (int l = 0; l < numLevels; l++) {
    RCP<Level> level                                                       = regHierarchy->GetLevel(l);
    RCP<Matrix> regMat                                                     = level->Get<RCP<Matrix>>("A");
    RCP<const Map> regRowMap                                               = regMat->getRowMap();
    RCP<Xpetra::Import<LocalOrdinal, GlobalOrdinal, Node>> regRowImporters = level->Get<RCP<Xpetra::Import<LocalOrdinal, GlobalOrdinal, Node>>>("rowImport");
    // initialize region vector with all ones.
    RCP<Xpetra::Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node>> regInterfaceScalings = VectorFactory::Build(regRowMap);
    regInterfaceScalings->putScalar(SC_ONE);

    // transform to composite layout while adding interface values via the Export() combine mode
    RCP<Vector> compInterfaceScalingSum = VectorFactory::Build(regRowImporters->getSourceMap(), true);
    regionalToComposite(regInterfaceScalings, compInterfaceScalingSum, regRowImporters);

    /* transform composite layout back to regional layout. Now, GIDs associated
     * with region interface should carry a scaling factor (!= 1).
     */
    RCP<Vector> quasiRegInterfaceScaling;  // Is that vector really needed?
    compositeToRegional(compInterfaceScalingSum, quasiRegInterfaceScaling,
                        regInterfaceScalings,
                        regRowMap, regRowImporters);

    level->Set<RCP<Xpetra::Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node>>>("regInterfaceScalings", regInterfaceScalings);
  }
}  // MakeInterfaceScalingFactors

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void createRegionHierarchy(const int numDimensions,
                           const Array<int> lNodesPerDim,
                           const std::string aggregationRegionType,
                           RCP<Teuchos::ParameterList>& interfaceParams,
                           const int maxRegPerGID,
                           RCP<Teuchos::ParameterList>& coarseSolverData,
                           Array<RCP<Teuchos::ParameterList>>& smootherParams,
                           RCP<Teuchos::ParameterList> hierarchyData,
                           RCP<MueLu::Hierarchy<Scalar, LocalOrdinal, GlobalOrdinal, Node>>& regHierarchy,
                           const bool keepCoarseCoords) {
#include "Xpetra_UseShortNames.hpp"
#include "MueLu_UseShortNames.hpp"
  using Teuchos::TimeMonitor;
  // This monitor times everything and gets the overall setting cost
  RCP<TimeMonitor> tm = rcp(new TimeMonitor(*TimeMonitor::getNewTimer("createRegionHierarchy")));

  using Hierarchy          = MueLu::Hierarchy<SC, LO, GO, NO>;
  using Utilities          = MueLu::Utilities<SC, LO, GO, NO>;
  using DirectCoarseSolver = Amesos2::Solver<Tpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>, Tpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>>;

  // std::cout << mapComp->getComm()->getRank() << " | Setting up MueLu hierarchies ..." << std::endl;
  int numLevels = 0;

  RCP<TimeMonitor> tmLocal = rcp(new TimeMonitor(*TimeMonitor::getNewTimer("createRegionHierarchy: ExtractData")));
  // resize Arrays and vectors
  {
    // resize level containers
    numLevels = regHierarchy->GetNumLevels();
    smootherParams.resize(numLevels);

    // resize group containers on each level
    for (int l = 0; l < numLevels; ++l) {
      // Also doing some initialization in the smootherParams
      if (l > 0) {
        smootherParams[l] = rcp(new Teuchos::ParameterList(*smootherParams[0]));
      }
    }
  }

  RCP<Level> level0            = regHierarchy->GetLevel(0);
  RCP<Matrix> regMat           = level0->Get<RCP<Matrix>>("A");
  RCP<const Map> revisedRowMap = regMat->getRowMap();
  /* Get coarse level matrices and prolongators from MueLu hierarchy
   * Note: fine level has been dealt with previously, so we start at level 1 here.
   */
  using real_type                  = typename Teuchos::ScalarTraits<Scalar>::coordinateType;
  using realvaluedmultivector_type = Xpetra::MultiVector<real_type, LocalOrdinal, GlobalOrdinal, Node>;
  RCP<realvaluedmultivector_type> regCoarseCoordinates;
  for (int l = 1; l < numLevels; ++l) {  // Note: we start at level 1 (which is the first coarse level)
    RCP<Level> level = regHierarchy->GetLevel(l);

    if (keepCoarseCoords && (l == numLevels - 1)) {
      regCoarseCoordinates = level->Get<RCP<realvaluedmultivector_type>>("Coordinates2", MueLu::NoFactory::get());
    }

    RCP<Matrix> regMatrices = level->Get<RCP<Matrix>>("A", MueLu::NoFactory::get());

    // Create residual and solution vectors and cache them for vCycle apply
    std::string levelName("level");
    levelName += std::to_string(l);
    ParameterList& levelList = hierarchyData->sublist(levelName, false, "list of data on current level");
    RCP<Vector> regRes       = VectorFactory::Build(revisedRowMap, true);
    RCP<Vector> regSol       = VectorFactory::Build(revisedRowMap, true);

    levelList.set<RCP<Vector>>("residual", regRes, "Cached residual vector");
    levelList.set<RCP<Vector>>("solution", regSol, "Cached solution vector");
  }

  // std::cout << mapComp->getComm()->getRank() << " | MakeCoarseLevelMaps ..." << std::endl;

  tmLocal = Teuchos::null;
  tmLocal = rcp(new TimeMonitor(*TimeMonitor::getNewTimer("createRegionHierarchy: MakeCoarseLevel")));

  MakeCoarseLevelMaps(maxRegPerGID,
                      regHierarchy);

  tmLocal = Teuchos::null;
  tmLocal = rcp(new TimeMonitor(*TimeMonitor::getNewTimer("createRegionHierarchy: MakeInterfaceScaling")));

  MakeInterfaceScalingFactors(numLevels,
                              regHierarchy);

  // std::cout << mapComp->getComm()->getRank() << " | Setup smoothers ..." << std::endl;

  tmLocal = Teuchos::null;
  tmLocal = rcp(new TimeMonitor(*TimeMonitor::getNewTimer("createRegionHierarchy: SmootherSetup")));
  // Only set the smoother up to the last but one level
  // if we want to use a smoother on the coarse level
  // we will handle that separately with "coarse solver type"
  for (int levelIdx = 0; levelIdx < numLevels - 1; ++levelIdx) {
    RCP<Level> level                                                                    = regHierarchy->GetLevel(levelIdx);
    RCP<Matrix> regMatrix                                                               = level->Get<RCP<Matrix>>("A", MueLu::NoFactory::get());
    RCP<const Map> regRowMap                                                            = regMatrix->getRowMap();
    RCP<Xpetra::Import<LocalOrdinal, GlobalOrdinal, Node>> regRowImporter               = level->Get<RCP<Xpetra::Import<LocalOrdinal, GlobalOrdinal, Node>>>("rowImport");
    RCP<Xpetra::Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node>> regInterfaceScalings = level->Get<RCP<Xpetra::Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node>>>("regInterfaceScalings");

    smootherParams[levelIdx]->set("smoother: level", levelIdx);
    smootherSetup(smootherParams[levelIdx], regRowMap,
                  regMatrix, regInterfaceScalings,
                  regRowImporter);
  }

  // std::cout << mapComp->getComm()->getRank() << " | CreateCoarseSolver ..." << std::endl;

  tmLocal                            = Teuchos::null;
  tmLocal                            = rcp(new TimeMonitor(*TimeMonitor::getNewTimer("createRegionHierarchy: CreateCoarseSolver")));
  const std::string coarseSolverType = coarseSolverData->get<std::string>("coarse solver type");
  if (coarseSolverType == "smoother") {
    RCP<Level> level                 = regHierarchy->GetLevel(numLevels - 1);
    RCP<Matrix> regMatrix            = level->Get<RCP<Matrix>>("A", MueLu::NoFactory::get());
    RCP<const Map> regRowMap         = regMatrix->getRowMap();
    RCP<Import> regRowImporter       = level->Get<RCP<Import>>("rowImport");
    RCP<Vector> regInterfaceScalings = level->Get<RCP<Vector>>("regInterfaceScalings");

    // Set the smoother on the coarsest level.
    const std::string smootherXMLFileName   = coarseSolverData->get<std::string>("smoother xml file");
    RCP<ParameterList> coarseSmootherParams = smootherParams[numLevels - 1];
    Teuchos::updateParametersFromXmlFileAndBroadcast(smootherXMLFileName, coarseSmootherParams.ptr(), *level->GetComm());
    coarseSmootherParams->set("smoother: level", numLevels - 1);
    coarseSmootherParams->print();

    smootherSetup(smootherParams[numLevels - 1], regRowMap,
                  regMatrix, regInterfaceScalings,
                  regRowImporter);
  } else if ((coarseSolverType == "direct") || (coarseSolverType == "amg")) {
    // A composite coarse matrix is needed

    // std::cout << mapComp->getComm()->getRank() << " | MakeCoarseCompositeOperator ..." << std::endl;

    RCP<Level> level              = regHierarchy->GetLevel(numLevels - 1);
    RCP<Matrix> regMatrix         = level->Get<RCP<Matrix>>("A", MueLu::NoFactory::get());
    RCP<Import> regRowImporter    = level->Get<RCP<Import>>("rowImport");
    RCP<const Map> compRowMap     = regRowImporter->getSourceMap();
    RCP<const Map> quasiRegRowMap = regRowImporter->getTargetMap();
    RCP<const Map> quasiRegColMap = regRowImporter->getTargetMap();  // col map same as row map.

    RCP<Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>> coarseCompOp;
    RCP<realvaluedmultivector_type> compCoarseCoordinates;
    MakeCoarseCompositeOperator(compRowMap,
                                quasiRegRowMap,
                                quasiRegColMap,
                                regRowImporter,
                                regMatrix,
                                regCoarseCoordinates,
                                coarseCompOp,
                                compCoarseCoordinates,
                                keepCoarseCoords);

    coarseSolverData->set<RCP<const Map>>("compCoarseRowMap", coarseCompOp->getRowMap());

    // std::cout << mapComp->getComm()->getRank() << " | MakeCoarseCompositeSolver ..." << std::endl;
    if (coarseSolverType == "direct") {
      RCP<DirectCoarseSolver> coarseDirectSolver = MakeCompositeDirectSolver(coarseCompOp);
      coarseSolverData->set<RCP<DirectCoarseSolver>>("direct solver object", coarseDirectSolver);
    } else if (coarseSolverType == "amg") {
      if (keepCoarseCoords == false) {
        std::cout << "WARNING: you requested a coarse AMG solver but you did not request coarse coordinates to be kept, repartitioning is not possible!" << std::endl;
      }

      RCP<Hierarchy> coarseAMGHierarchy;
      std::string amgXmlFileName = coarseSolverData->get<std::string>("amg xml file");
#if defined(HAVE_MUELU_ZOLTAN2) && defined(HAVE_MPI)
      const bool coarseSolverRebalance = coarseSolverData->get<bool>("coarse solver rebalance");
      if (keepCoarseCoords == true && coarseSolverRebalance == true) {
        RCP<Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>> rebalancedCompOp;
        RCP<Xpetra::MultiVector<typename Teuchos::ScalarTraits<Scalar>::coordinateType, LocalOrdinal, GlobalOrdinal, Node>> rebalancedCoordinates;
        RCP<const Import> rebalanceImporter;

        const int rebalanceNumPartitions = coarseSolverData->get<int>("coarse rebalance num partitions");
        RebalanceCoarseCompositeOperator(rebalanceNumPartitions,
                                         coarseCompOp,
                                         compCoarseCoordinates,
                                         rebalancedCompOp,
                                         rebalancedCoordinates,
                                         rebalanceImporter);
        coarseSolverData->set<RCP<const Import>>("rebalanceImporter", rebalanceImporter);

        if (!rebalancedCompOp.is_null())
          coarseAMGHierarchy = MakeCompositeAMGHierarchy(rebalancedCompOp, amgXmlFileName, rebalancedCoordinates);

      } else {
        coarseAMGHierarchy = MakeCompositeAMGHierarchy(coarseCompOp, amgXmlFileName, compCoarseCoordinates);
      }
#else
      coarseAMGHierarchy = MakeCompositeAMGHierarchy(coarseCompOp, amgXmlFileName, compCoarseCoordinates);
#endif
      coarseSolverData->set<RCP<Hierarchy>>("amg hierarchy object", coarseAMGHierarchy);
    }
  } else {
    TEUCHOS_TEST_FOR_EXCEPT_MSG(false, "Unknown coarse solver type.");
  }

}  // createRegionHierarchy

#endif  // MUELU_SETUPREGIONHIERARCHY_DEF_HPP
