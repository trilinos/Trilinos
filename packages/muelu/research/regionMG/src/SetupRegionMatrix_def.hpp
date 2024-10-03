// @HEADER
// *****************************************************************************
//        MueLu: A package for multigrid based preconditioning
//
// Copyright 2012 NTESS and the MueLu contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef MUELU_SETUPREGIONMATRIX_DEF_HPP
#define MUELU_SETUPREGIONMATRIX_DEF_HPP

#include <vector>
#include <iostream>

#define RegionsSpanProcs 1
#define MultipleRegionsPerProc 2

#include <Tpetra_KokkosCompat_DefaultNode.hpp>
#include <KokkosSparse_spmv.hpp>

#include <Xpetra_ConfigDefs.hpp>
#include <Xpetra_Export.hpp>
#include <Xpetra_Import.hpp>
#include <Xpetra_Map.hpp>
#include <Xpetra_MultiVector.hpp>
#include <Xpetra_Vector.hpp>
#include <Xpetra_CrsMatrixWrap.hpp>

using Teuchos::Array;
using Teuchos::ArrayRCP;
using Teuchos::RCP;

/*! \brief Find common regions of two nodes
 *
 */
template <class LocalOrdinal, class GlobalOrdinal, class Node>
Teuchos::Array<int> findCommonRegions(const GlobalOrdinal nodeA,                                 ///< GID of first node
                                      const GlobalOrdinal nodeB,                                 ///< GID of second node
                                      const Array<ArrayRCP<const LocalOrdinal>> nodesToRegions,  ///< mapping of nodes to regions
                                      RCP<const Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node>> nodesToRegionsMap) {
#include "Xpetra_UseShortNamesOrdinal.hpp"
  using Teuchos::TimeMonitor;
  RCP<TimeMonitor> tm = rcp(new TimeMonitor(*TimeMonitor::getNewTimer("findCommonRegions: 1 - Extract regions")));

  // extract node-to-regions mapping for both nodes A and B
  Array<int> regionsA, regionsB;
  {
    LO nodeALID = nodesToRegionsMap->getLocalElement(nodeA);
    LO nodeBLID = nodesToRegionsMap->getLocalElement(nodeB);
    for (int i = 0; i < nodesToRegions.size(); ++i) {
      regionsA.push_back(nodesToRegions[i][nodeALID]);
      regionsB.push_back(nodesToRegions[i][nodeBLID]);
    }
  }

  tm = Teuchos::null;
  tm = rcp(new TimeMonitor(*TimeMonitor::getNewTimer("findCommonRegions: 2 - Sort regions")));

  // identify common regions
  std::vector<int> commonRegions(nodesToRegions.size());
  std::sort(regionsA.begin(), regionsA.end());
  std::sort(regionsB.begin(), regionsB.end());

  tm = Teuchos::null;
  tm = rcp(new TimeMonitor(*TimeMonitor::getNewTimer("findCommonRegions: 3 - Find commons")));

  std::vector<int>::iterator it = std::set_intersection(regionsA.begin(),
                                                        regionsA.end(), regionsB.begin(), regionsB.end(), commonRegions.begin());
  commonRegions.resize(it - commonRegions.begin());

  tm = Teuchos::null;
  tm = rcp(new TimeMonitor(*TimeMonitor::getNewTimer("findCommonRegions: 4 - Clean-up output")));

  // remove '-1' entries
  Teuchos::Array<int> finalCommonRegions;
  for (std::size_t i = 0; i < commonRegions.size(); ++i) {
    if (commonRegions[i] != -1)
      finalCommonRegions.push_back(commonRegions[i]);
  }

  tm = Teuchos::null;

  return finalCommonRegions;
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void MakeQuasiregionMatrices(const RCP<Xpetra::CrsMatrixWrap<Scalar, LocalOrdinal, GlobalOrdinal, Node>> AComp,
                             RCP<Xpetra::MultiVector<LocalOrdinal, LocalOrdinal, GlobalOrdinal, Node>> regionsPerGIDWithGhosts,
                             RCP<const Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node>> rowMap,
                             RCP<const Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node>> colMap,
                             RCP<Xpetra::Import<LocalOrdinal, GlobalOrdinal, Node>>& rowImport,
                             RCP<Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>>& quasiRegionMats,
                             const Teuchos::ArrayRCP<LocalOrdinal>& regionMatVecLIDs) {
#include "Xpetra_UseShortNames.hpp"
  using Teuchos::RCP;
  using Teuchos::TimeMonitor;

  Array<ArrayRCP<const LO>> regionPerGIDWithGhostsData(regionsPerGIDWithGhosts->getNumVectors());
  for (size_t vecIdx = 0; vecIdx < regionsPerGIDWithGhosts->getNumVectors(); ++vecIdx) {
    regionPerGIDWithGhostsData[vecIdx] = regionsPerGIDWithGhosts->getData(vecIdx);
  }

  /* We use the edge-based splitting, i.e. we first modify off-diagonal
   * entries in the composite matrix, then decompose it into region matrices
   * and finally take care of diagonal entries by enforcing the nullspace
   * preservation constraint.
   */

  // Import data from AComp into the quasiRegion matrices
  // Since the stencil size does not grow between composite and region format
  // use AComp->getCrsGraph()->getLocalMaxNumRowEntries() to get an upper
  // bound of the number of nonzeros per row in quasiRegionMat
  RCP<TimeMonitor> tm = rcp(new TimeMonitor(*TimeMonitor::getNewTimer("MakeQuasiregionMatrices: 1 - Create Matrix")));

  quasiRegionMats = MatrixFactory::Build(rowMap,
                                         colMap,
                                         AComp->getCrsGraph()->getLocalMaxNumRowEntries());

  tm = Teuchos::null;
  tm = rcp(new TimeMonitor(*TimeMonitor::getNewTimer("MakeQuasiregionMatrices: 2 - Import data")));

  quasiRegionMats->doImport(*AComp,
                            *(rowImport),
                            Xpetra::INSERT);

  tm = Teuchos::null;
  tm = rcp(new TimeMonitor(*TimeMonitor::getNewTimer("MakeQuasiregionMatrices: 3 - Scale interface entries")));

  RCP<CrsMatrixWrap> quasiRegionCrsWrap = Teuchos::rcp_dynamic_cast<CrsMatrixWrap>(quasiRegionMats);
  RCP<CrsMatrix> quasiRegionCrs         = quasiRegionCrsWrap->getCrsMatrix();

  // Grab first and last element of sorted interface LIDs
  Array<LO> interfaceLIDs(regionMatVecLIDs());
  std::sort(interfaceLIDs.begin(), interfaceLIDs.end());
  auto vecEnd   = std::unique(interfaceLIDs.begin(), interfaceLIDs.end());
  auto vecStart = interfaceLIDs.begin();

  GO rowGID;
  LocalOrdinal col;
  GlobalOrdinal colGID;
  std::size_t sizeOfCommonRegions;
  std::size_t numEntries = 0;
  for (auto row = vecStart; row < vecEnd; ++row) {
    rowGID     = rowMap->getGlobalElement(*row);
    numEntries = quasiRegionMats->getNumEntriesInLocalRow(*row);  // number of entries in this row
    Array<SC> values(numEntries);                                 // non-zeros in this row
    Array<LO> colInds(numEntries);                                // local column indices
    quasiRegionMats->getLocalRowCopy(*row, colInds, values, numEntries);

    for (std::size_t entryIdx = 0; entryIdx < numEntries; ++entryIdx) {  // loop over all entries in this row
      col    = colInds[entryIdx];
      colGID = colMap->getGlobalElement(col);
      Array<int> commonRegions;
      if (rowGID != colGID) {  // Skip the diagonal entry. It will be processed later.
        commonRegions = findCommonRegions(rowGID, colGID, regionPerGIDWithGhostsData, regionsPerGIDWithGhosts->getMap());
      }

      sizeOfCommonRegions = commonRegions.size();
      if (1 < sizeOfCommonRegions) {
        values[entryIdx] /= Teuchos::as<double>(sizeOfCommonRegions);
      }
    }
    quasiRegionMats->replaceLocalValues(*row, colInds, values);
  }

  tm = Teuchos::null;
  tm = rcp(new TimeMonitor(*TimeMonitor::getNewTimer("MakeQuasiregionMatrices: 4 - fillComplete")));

  quasiRegionMats->fillComplete();

  tm = Teuchos::null;
}  // MakeQuasiregionMatrices

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void MakeRegionMatrices(const RCP<Xpetra::CrsMatrixWrap<Scalar, LocalOrdinal, GlobalOrdinal, Node>> AComp,
                        const RCP<const Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node>> mapComp,
                        RCP<const Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node>> rowMap,
                        RCP<const Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node>> revisedRowMap,
                        RCP<const Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node>> revisedColMap,
                        RCP<Xpetra::Import<LocalOrdinal, GlobalOrdinal, Node>>& rowImport,
                        RCP<Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>>& quasiRegionMats,
                        RCP<Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>>& regionMats) {
#include "Xpetra_UseShortNames.hpp"
  using Teuchos::RCP;
  using Teuchos::TimeMonitor;
  const SC SC_ONE  = Teuchos::ScalarTraits<SC>::one();
  const SC SC_ZERO = Teuchos::ScalarTraits<SC>::zero();

  RCP<TimeMonitor> tm = rcp(new TimeMonitor(*TimeMonitor::getNewTimer("MakeRegionMatrices: 1 - Create Matrix")));

  // Copy data from quasiRegionMats, but into new map layout
  {
    regionMats = rcp(new CrsMatrixWrap(revisedRowMap, revisedColMap, 9));

    // Extract current region CrsMatrix
    RCP<CrsMatrix> regionCrsMat = Teuchos::rcp_dynamic_cast<CrsMatrixWrap>(regionMats)->getCrsMatrix();

    // Extract current quasi-region CrsMatrix
    RCP<CrsMatrix> quasiRegionCrsMat = Teuchos::rcp_dynamic_cast<CrsMatrixWrap>(quasiRegionMats)->getCrsMatrix();

    // Pull out the data from the quasi-region CrsMatrix
    ArrayRCP<const size_t> rowptrQuasiRegion;
    ArrayRCP<const LocalOrdinal> colindQuasiRegion;
    ArrayRCP<const Scalar> valuesQuasiRegion;
    quasiRegionCrsMat->getAllValues(rowptrQuasiRegion, colindQuasiRegion, valuesQuasiRegion);

    // Do a deep copy of values
    // (at least we've been doing deep copies so far, maybe we could do shallow copies to save time?)
    ArrayRCP<size_t> rowptrRegion(rowptrQuasiRegion.size());
    ArrayRCP<LocalOrdinal> colindRegion(colindQuasiRegion.size());
    ArrayRCP<Scalar> valuesRegion(valuesQuasiRegion.size());

    regionCrsMat->allocateAllValues(valuesRegion.size(), rowptrRegion, colindRegion, valuesRegion);

    for (LocalOrdinal idx = 0; idx < static_cast<LocalOrdinal>(rowptrRegion.size()); ++idx) {
      rowptrRegion[idx] = rowptrQuasiRegion[idx];
    }

    for (LocalOrdinal idx = 0; idx < static_cast<LocalOrdinal>(colindRegion.size()); ++idx) {
      colindRegion[idx] = colindQuasiRegion[idx];
      valuesRegion[idx] = valuesQuasiRegion[idx];
    }

    // Set and fillComplete the region CrsMatrix
    regionCrsMat->setAllValues(rowptrRegion, colindRegion, valuesRegion);
    regionCrsMat->expertStaticFillComplete(revisedRowMap, revisedRowMap);
  }

  tm = Teuchos::null;
  tm = rcp(new TimeMonitor(*TimeMonitor::getNewTimer("MakeRegionMatrices: 2 - Enforce nullspace constraint")));

  // enforce nullspace constraint
  RCP<Vector> regNspViolation;
  {
    // compute violation of nullspace property close to DBCs
    RCP<Vector> nspVec = VectorFactory::Build(mapComp);
    nspVec->putScalar(SC_ONE);
    RCP<Vector> nspViolation = VectorFactory::Build(mapComp, true);
    AComp->apply(*nspVec, *nspViolation);

    // move to regional layout
    RCP<Vector> quasiRegNspViolation = VectorFactory::Build(rowMap, true);
    regNspViolation                  = VectorFactory::Build(revisedRowMap, true);
    compositeToRegional(nspViolation, quasiRegNspViolation, regNspViolation,
                        revisedRowMap, rowImport);

    /* The nullspace violation computed in the composite layout needs to be
     * transfered to the regional layout. Since we use this to compute
     * the splitting of the diagonal values, we need to split the nullspace
     * violation. We'd like to use the 'regInterfaceScaling', though this is
     * not setup at this point. So, let's compute it right now.
     *
     * ToDo: Move setup of 'regInterfaceScaling' up front to use it here.
     */
    {
      // initialize region vector with all ones.
      RCP<Vector> interfaceScaling = VectorFactory::Build(revisedRowMap);
      interfaceScaling->putScalar(SC_ONE);

      // transform to composite layout while adding interface values via the Export() combine mode
      RCP<Vector> compInterfaceScalingSum = VectorFactory::Build(mapComp, true);
      regionalToComposite(interfaceScaling, compInterfaceScalingSum, rowImport);

      /* transform composite layout back to regional layout. Now, GIDs associated
       * with region interface should carry a scaling factor (!= 1).
       */
      RCP<Vector> quasiRegInterfaceScaling;
      compositeToRegional(compInterfaceScalingSum, quasiRegInterfaceScaling,
                          interfaceScaling,
                          revisedRowMap, rowImport);

      // modify its interface entries
      RCP<Vector> inverseInterfaceScaling = VectorFactory::Build(interfaceScaling->getMap(), true);
      inverseInterfaceScaling->reciprocal(*interfaceScaling);
      regNspViolation->elementWiseMultiply(SC_ONE, *regNspViolation, *inverseInterfaceScaling, SC_ZERO);
    }
  }

  tm = Teuchos::null;
  tm = rcp(new TimeMonitor(*TimeMonitor::getNewTimer("MakeRegionMatrices: 3 - Replace diagonal")));

  RCP<Vector> regNsp;
  regNsp = VectorFactory::Build(revisedRowMap);
  regNsp->putScalar(SC_ONE);

  RCP<Vector> regCorrection;
  regCorrection = VectorFactory::Build(revisedRowMap, true);
  regionMats->apply(*regNsp, *regCorrection);
  regionMats->SetFixedBlockSize(AComp->GetFixedBlockSize());

  RCP<Vector> regDiag = Teuchos::null;
  regDiag             = VectorFactory::Build(revisedRowMap, true);
  regionMats->getLocalDiagCopy(*regDiag);
  regDiag->update(-SC_ONE, *regCorrection, SC_ONE, *regNspViolation, SC_ONE);

  // Extract current region matrix in as CrsMatrix
  RCP<CrsMatrix> regionCrsMat = Teuchos::rcp_dynamic_cast<CrsMatrixWrap>(regionMats)->getCrsMatrix();
  regionCrsMat->replaceDiag(*regDiag);

  tm = Teuchos::null;
}  // MakeRegionMatrices

/*! \brief Transform regional matrix to composite layout
 *
 *  Starting from a \c Matrix in regional layout, we
 *  1. copy data from regional matrix into quasiRegional matrix and set all maps
 *     to be quasiRegional maps
 *  2. export it into a \c Matrix with composite layout using the given \c combineMode.
 *     Note: on-process values also have to be added to account for region interfaces inside a process.
 *
 *  \note We also need the capability to add processor-local values. This is not supported by
 *  available CombineMode options in Xpetra/Tpetra, so we use a manual implementation here.
 *
 *  \return Composite matrix that is fill-completed
 */
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void regionalToComposite(const RCP<Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>>& regMat,  ///< Matrix in region layout [in]
                         const RCP<const Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node>> rowMap,        ///< row maps in quasiRegion layout [in]
                         const RCP<const Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node>> colMap,        ///< col maps in quasiRegion layout [in]
                         const RCP<Xpetra::Import<LocalOrdinal, GlobalOrdinal, Node>> rowImport,        ///< row importer in region layout [in]
                         const Xpetra::CombineMode combineMode,                                         ///< Combine mode for import/export [in]
                         RCP<Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>>& compMat        ///< Matrix in composite layout [in/out]
) {
#include "Xpetra_UseShortNames.hpp"
  using std::size_t;
  using Teuchos::rcp;
  using Teuchos::TimeMonitor;

  RCP<TimeMonitor> tm = rcp(new TimeMonitor(*TimeMonitor::getNewTimer("regionalToComposite: Matrix")));

  const SC SC_ZERO = Teuchos::ScalarTraits<Scalar>::zero();
  const SC SC_ONE  = Teuchos::ScalarTraits<Scalar>::one();

  // Make sure we add into an zero composite matrix
  compMat->setAllToScalar(SC_ZERO);

  /* Let's fake an ADD combine mode that also adds local values by
   * 1. exporting quasiRegional matrices to auxiliary composite matrices (1 per group)
   * 2. add all auxiliary matrices together
   */

  // Copy data from regMat into quasiRegMat
  RCP<Matrix> quasiRegMat;
  {
    quasiRegMat = rcp(new CrsMatrixWrap(rowMap,
                                        colMap,
                                        regMat->getCrsGraph()->getLocalMaxNumRowEntries()));

    // Extract current quasi-region CrsMatrix
    RCP<CrsMatrix> quasiRegionCrsMat = Teuchos::rcp_dynamic_cast<CrsMatrixWrap>(quasiRegMat)->getCrsMatrix();

    // Extract current region CrsMatrix
    RCP<CrsMatrix> regionCrsMat = Teuchos::rcp_dynamic_cast<CrsMatrixWrap>(regMat)->getCrsMatrix();

    // Pull out the data from the region CrsMatrix
    Teuchos::ArrayRCP<const size_t> rowptrRegion;
    Teuchos::ArrayRCP<const LocalOrdinal> colindRegion;
    Teuchos::ArrayRCP<const Scalar> valuesRegion;
    regionCrsMat->getAllValues(rowptrRegion, colindRegion, valuesRegion);

    // Do a deep copy of values
    // (at least we've been doing deep copies so far, maybe we could do shallow copies to save time?)
    Teuchos::ArrayRCP<size_t> rowptrQuasiRegion(rowptrRegion.size());
    Teuchos::ArrayRCP<LocalOrdinal> colindQuasiRegion(colindRegion.size());
    Teuchos::ArrayRCP<Scalar> valuesQuasiRegion(valuesRegion.size());

    quasiRegionCrsMat->allocateAllValues(valuesQuasiRegion.size(), rowptrQuasiRegion, colindQuasiRegion, valuesQuasiRegion);

    for (LocalOrdinal idx = 0; idx < static_cast<LocalOrdinal>(rowptrQuasiRegion.size()); ++idx) {
      rowptrQuasiRegion[idx] = rowptrRegion[idx];
    }

    for (LocalOrdinal idx = 0; idx < static_cast<LocalOrdinal>(colindQuasiRegion.size()); ++idx) {
      colindQuasiRegion[idx] = colindRegion[idx];
      valuesQuasiRegion[idx] = valuesRegion[idx];
    }

    // Set and fillComplete the quasiRegion CrsMatrix
    quasiRegionCrsMat->setAllValues(rowptrQuasiRegion, colindQuasiRegion, valuesQuasiRegion);
    quasiRegionCrsMat->expertStaticFillComplete(rowMap, rowMap);
  }

  // Export from quasiRegional format to composite layout
  RCP<Matrix> partialCompMat;
  partialCompMat = MatrixFactory::Build(compMat->getRowMap(),
                                        8 * regMat->getCrsGraph()->getLocalMaxNumRowEntries());
  partialCompMat->doExport(*(quasiRegMat), *(rowImport), Xpetra::INSERT);
  partialCompMat->fillComplete();

  // Add all partialCompMat together
  MatrixMatrix::TwoMatrixAdd(*partialCompMat, false, SC_ONE, *compMat, SC_ONE);

  compMat->fillComplete();

  return;
}  // regionalToComposite

/*! \brief Compute local data needed to perform a MatVec in region format

 */
template <class LocalOrdinal, class GlobalOrdinal, class Node>
void SetupMatVec(const Teuchos::RCP<Xpetra::MultiVector<GlobalOrdinal, LocalOrdinal, GlobalOrdinal, Node>>& interfaceGIDsMV,
                 const Teuchos::RCP<Xpetra::MultiVector<LocalOrdinal, LocalOrdinal, GlobalOrdinal, Node>>& regionsPerGIDWithGhosts,
                 const Teuchos::RCP<const Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node>>& regionRowMap,
                 const Teuchos::RCP<Xpetra::Import<LocalOrdinal, GlobalOrdinal, Node>>& rowImport,
                 Teuchos::ArrayRCP<LocalOrdinal>& regionMatVecLIDs,
                 Teuchos::RCP<Xpetra::Import<LocalOrdinal, GlobalOrdinal, Node>>& regionInterfaceImporter) {
#include "Xpetra_UseShortNamesOrdinal.hpp"
  using Teuchos::TimeMonitor;

  RCP<TimeMonitor> tm;
  tm = rcp(new TimeMonitor(*TimeMonitor::getNewTimer("SetupMatVec: 1 - sumInterfaceValues")));

  const LO maxRegPerGID = static_cast<LO>(regionsPerGIDWithGhosts->getNumVectors());
  const int myRank      = regionRowMap->getComm()->getRank();
  interfaceGIDsMV->replaceMap(regionRowMap);
  RCP<Xpetra::MultiVector<GO, LO, GO, NO>> interfaceGIDs;
  interfaceGIDs = interfaceGIDsMV;
  sumInterfaceValues(interfaceGIDs, regionRowMap, rowImport);

  tm = Teuchos::null;
  tm = rcp(new TimeMonitor(*TimeMonitor::getNewTimer("SetupMatVec: 2 - build regionInterfaceMap")));

  Teuchos::Array<LO> regionMatVecLIDstmp;
  Teuchos::Array<GO> regionMatVecGIDs;
  Array<ArrayRCP<const LO>> regionsPerGIDWithGhostsData(maxRegPerGID);
  Array<ArrayRCP<const GO>> interfaceGIDsData(maxRegPerGID);
  for (LO regionIdx = 0; regionIdx < maxRegPerGID; ++regionIdx) {
    regionsPerGIDWithGhostsData[regionIdx] = regionsPerGIDWithGhosts->getData(regionIdx);
    interfaceGIDsData[regionIdx]           = interfaceGIDs->getData(regionIdx);
    for (LO idx = 0; idx < static_cast<LO>(regionsPerGIDWithGhostsData[regionIdx].size()); ++idx) {
      if ((regionsPerGIDWithGhostsData[regionIdx][idx] != -1) && (regionsPerGIDWithGhostsData[regionIdx][idx] != myRank)) {
        regionMatVecLIDstmp.push_back(idx);
        regionMatVecGIDs.push_back(interfaceGIDsData[regionIdx][idx]);
      }
    }
  }

  // Copy the temporary regionMatVecLIDstmp into an ArrayRCP
  // so we can store it and retrieve it easily later on.
  regionMatVecLIDs.deepCopy(regionMatVecLIDstmp());

  RCP<Map> regionInterfaceMap = Xpetra::MapFactory<LO, GO, Node>::Build(regionRowMap->lib(),
                                                                        Teuchos::OrdinalTraits<GO>::invalid(),
                                                                        regionMatVecGIDs(),
                                                                        regionRowMap->getIndexBase(),
                                                                        regionRowMap->getComm());

  tm = Teuchos::null;
  tm = rcp(new TimeMonitor(*TimeMonitor::getNewTimer("SetupMatVec: 3 - Build importer")));

  regionInterfaceImporter = ImportFactory::Build(regionRowMap, regionInterfaceMap);

  tm = Teuchos::null;
}  // SetupMatVec

/*! \brief Compute the residual \f$r = b - Ax\f$ with pre-computed communication patterns
 *
 *  The residual is computed based on matrices and vectors in a regional layout.
 *  1. Compute and sum interface values of y = A*x in regional layout using the fast MatVev
 *  2. Compute residual via r = b - y
 */
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void computeResidual(RCP<Xpetra::Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node>>& regRes,           ///< residual (to be evaluated)
                     const RCP<Xpetra::Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node>> regX,        ///< left-hand side (solution)
                     const RCP<Xpetra::Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node>> regB,        ///< right-hand side (forcing term)
                     const RCP<Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>> regionMats,  ///< matrix in region format
                     const Teuchos::ParameterList& params                                              ///< parameter with fast MatVec parameters and pre-computed communication patterns
) {
#include "Xpetra_UseShortNames.hpp"
  using TST = Teuchos::ScalarTraits<Scalar>;
  using Teuchos::TimeMonitor;

  RCP<TimeMonitor> tm = rcp(new TimeMonitor(*TimeMonitor::getNewTimer("computeResidual: use fast MatVec")));

  // Get pre-communicated communication patterns for the fast MatVec
  const ArrayRCP<LocalOrdinal> regionInterfaceLIDs = params.get<ArrayRCP<LO>>("Fast MatVec: interface LIDs");
  const RCP<Import> regionInterfaceImporter        = params.get<RCP<Import>>("Fast MatVec: interface importer");

  // Step 1: Compute region version of y = Ax and store it in regRes
  regionMats->apply(*regX, *regRes, Teuchos::NO_TRANS, TST::one(), TST::zero(), true, regionInterfaceImporter, regionInterfaceLIDs);

  // Step 2: Compute region version of r = b - y
  regRes->update(TST::one(), *regB, -TST::one(), *regRes, TST::zero());

  tm = Teuchos::null;
}  // computeResidual

#endif  // MUELU_SETUPREGIONMATRIX_DEF_HPP
