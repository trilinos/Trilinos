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
//                    Ray Tuminaro      (rstumin@sandia.gov)
//
// ***********************************************************************
//
// @HEADER
#ifndef MUELU_SETUPREGIONMATRIX_DEF_HPP
#define MUELU_SETUPREGIONMATRIX_DEF_HPP

#include <vector>
#include <iostream>

#define RegionsSpanProcs  1
#define MultipleRegionsPerProc  2

#include <Kokkos_DefaultNode.hpp>
#include <KokkosSparse_spmv.hpp>

#include <Xpetra_ConfigDefs.hpp>
#include <Xpetra_Export.hpp>
#include <Xpetra_Import.hpp>
#include <Xpetra_Map.hpp>
#include <Xpetra_MultiVector.hpp>
#include <Xpetra_Vector.hpp>
#include <Xpetra_CrsMatrixWrap.hpp>

using Teuchos::RCP;
using Teuchos::ArrayRCP;
using Teuchos::Array;

/*! \brief Find common regions of two nodes
 *
 */
template<class LocalOrdinal, class GlobalOrdinal, class Node>
Teuchos::Array<int> findCommonRegions(const GlobalOrdinal nodeA, ///< GID of first node
                                      const GlobalOrdinal nodeB, ///< GID of second node
                                      const Array<ArrayRCP<const LocalOrdinal> > nodesToRegions, ///< mapping of nodes to regions
                                      RCP<const Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node> > nodesToRegionsMap
                                      )
{
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
void MakeQuasiregionMatrices(const RCP<Xpetra::CrsMatrixWrap<Scalar, LocalOrdinal, GlobalOrdinal, Node> > AComp,
                             RCP<Xpetra::MultiVector<LocalOrdinal, LocalOrdinal, GlobalOrdinal, Node> > regionsPerGIDWithGhosts,
                             RCP<const Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node> > rowMap,
                             RCP<const Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node> > colMap,
                             RCP<Xpetra::Import<LocalOrdinal, GlobalOrdinal, Node> >& rowImport,
                             RCP<Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> >& quasiRegionMats) {
#include "Xpetra_UseShortNames.hpp"
  using Teuchos::RCP;
  using Teuchos::TimeMonitor;

  Array<ArrayRCP<const LO> > regionPerGIDWithGhostsData(regionsPerGIDWithGhosts->getNumVectors());
  for(size_t vecIdx = 0; vecIdx < regionsPerGIDWithGhosts->getNumVectors(); ++vecIdx) {
    regionPerGIDWithGhostsData[vecIdx] = regionsPerGIDWithGhosts->getData(vecIdx);
  }

  /* We use the edge-based splitting, i.e. we first modify off-diagonal
   * entries in the composite matrix, then decompose it into region matrices
   * and finally take care of diagonal entries by enforcing the nullspace
   * preservation constraint.
   */

  // Import data from AComp into the quasiRegion matrices
  // Since the stencil size does not grow between composite and region format
  // use AComp->getCrsGraph()->getNodeMaxNumRowEntries() to get an upper
  // bound of the number of nonzeros per row in quasiRegionMat
  RCP<TimeMonitor> tm = rcp(new TimeMonitor(*TimeMonitor::getNewTimer("MakeQuasiregionMatrices: 1 - Create Matrix")));

  quasiRegionMats = MatrixFactory::Build(rowMap,
                                         colMap,
                                         AComp->getCrsGraph()->getNodeMaxNumRowEntries());

  tm = Teuchos::null;
  tm = rcp(new TimeMonitor(*TimeMonitor::getNewTimer("MakeQuasiregionMatrices: 2 - Import data")));

  quasiRegionMats->doImport(*AComp,
                            *(rowImport),
                            Xpetra::INSERT);

  tm = Teuchos::null;
  tm = rcp(new TimeMonitor(*TimeMonitor::getNewTimer("MakeQuasiregionMatrices: 3 - Scale interface entries")));

  RCP<CrsMatrixWrap> quasiRegionCrsWrap = Teuchos::rcp_dynamic_cast<CrsMatrixWrap>(quasiRegionMats);
  RCP<CrsMatrix> quasiRegionCrs = quasiRegionCrsWrap->getCrsMatrix();
  const LO numRows = Teuchos::as<LO>(quasiRegionCrs->getNodeNumRows());

  for (LO row = 0; row < numRows; ++row) { // loop over local rows of composite matrix
    GO rowGID = rowMap->getGlobalElement(row);
    std::size_t numEntries = quasiRegionMats->getNumEntriesInLocalRow(row); // number of entries in this row
    Teuchos::Array<SC> vals(numEntries); // non-zeros in this row
    Teuchos::Array<LO> inds(numEntries); // local column indices
    quasiRegionMats->getLocalRowCopy(row, inds, vals, numEntries);

    for (std::size_t c = 0; c < Teuchos::as<std::size_t>(inds.size()); ++c) { // loop over all entries in this row
      LocalOrdinal col = inds[c];
      GlobalOrdinal colGID = colMap->getGlobalElement(col);
      Array<int> commonRegions;
      if (rowGID != colGID) { // Skip the diagonal entry. It will be processed later.
        // commonRegions = findCommonRegions(rowGID, colGID, *regionsPerGIDWithGhosts);
        commonRegions = findCommonRegions(rowGID, colGID, regionPerGIDWithGhostsData, regionsPerGIDWithGhosts->getMap());
      }

      std::size_t sizeOfCommonRegions = commonRegions.size();
      if (sizeOfCommonRegions > 1) {
        vals[c] /= Teuchos::as<double>(sizeOfCommonRegions);
      }
    }

    quasiRegionMats->replaceLocalValues(row, inds, vals);
  }

  tm = Teuchos::null;
  tm = rcp(new TimeMonitor(*TimeMonitor::getNewTimer("MakeQuasiregionMatrices: 4 - fillComplete")));

  quasiRegionMats->fillComplete();

  tm = Teuchos::null;
} // MakeQuasiregionMatrices


template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void MakeRegionMatrices(const RCP<Xpetra::CrsMatrixWrap<Scalar, LocalOrdinal, GlobalOrdinal, Node> > AComp,
                        const RCP<const Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node> > mapComp,
                        RCP<const Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node> > rowMap,
                        RCP<const Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node> > revisedRowMap,
                        RCP<const Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node> > revisedColMap,
                        RCP<Xpetra::Import<LocalOrdinal, GlobalOrdinal, Node> >& rowImport,
                        RCP<Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> >& quasiRegionMats,
                        RCP<Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> >& regionMats) {
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

      for(LocalOrdinal idx = 0; idx < static_cast<LocalOrdinal>(rowptrRegion.size()); ++idx) {
        rowptrRegion[idx] = rowptrQuasiRegion[idx];
      }

      for(LocalOrdinal idx = 0; idx < static_cast<LocalOrdinal>(colindRegion.size()); ++idx) {
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
    regNspViolation = VectorFactory::Build(revisedRowMap, true);
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
  regDiag = VectorFactory::Build(revisedRowMap, true);
  regionMats->getLocalDiagCopy(*regDiag);
  regDiag->update(-SC_ONE, *regCorrection, SC_ONE, *regNspViolation, SC_ONE);

  // Extract current region matrix in as CrsMatrix
  RCP<CrsMatrix> regionCrsMat = Teuchos::rcp_dynamic_cast<CrsMatrixWrap>(regionMats)->getCrsMatrix();
  regionCrsMat->replaceDiag(*regDiag);

  tm = Teuchos::null;
} // MakeRegionMatrices


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
void regionalToComposite(const RCP<Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> >& regMat, ///< Matrix in region layout [in]
                         const RCP<const Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node> > rowMap, ///< row maps in quasiRegion layout [in]
                         const RCP<const Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node> > colMap, ///< col maps in quasiRegion layout [in]
                         const RCP<Xpetra::Import<LocalOrdinal, GlobalOrdinal, Node> > rowImport, ///< row importer in region layout [in]
                         const Xpetra::CombineMode combineMode, ///< Combine mode for import/export [in]
                         RCP<Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> >& compMat ///< Matrix in composite layout [in/out]
                         )
{
#include "Xpetra_UseShortNames.hpp"
  using Teuchos::TimeMonitor;
  using Teuchos::rcp;
  using std::size_t;

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
                                           regMat->getCrsGraph()->getNodeMaxNumRowEntries()));

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

    for(LocalOrdinal idx = 0; idx < static_cast<LocalOrdinal>(rowptrQuasiRegion.size()); ++idx) {
      rowptrQuasiRegion[idx] = rowptrRegion[idx];
    }

    for(LocalOrdinal idx = 0; idx < static_cast<LocalOrdinal>(colindQuasiRegion.size()); ++idx) {
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
                                           8*regMat->getCrsGraph()->getNodeMaxNumRowEntries());
  partialCompMat->doExport(*(quasiRegMat), *(rowImport), Xpetra::INSERT);
  partialCompMat->fillComplete();

  // Add all partialCompMat together
  MatrixMatrix::TwoMatrixAdd(*partialCompMat, false, SC_ONE, *compMat, SC_ONE);

  compMat->fillComplete();

  return;
} // regionalToComposite


/*! \brief Compute local data needed to perform a MatVec in region format

 */
template <class LocalOrdinal, class GlobalOrdinal, class Node>
void SetupMatVec(const Teuchos::RCP<Xpetra::MultiVector<GlobalOrdinal, LocalOrdinal, GlobalOrdinal, Node> >& interfaceGIDsMV,
                 const Teuchos::RCP<Xpetra::MultiVector<LocalOrdinal, LocalOrdinal, GlobalOrdinal, Node> >& regionsPerGIDWithGhosts,
                 const Teuchos::RCP<const Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node> >& regionRowMap,
                 const Teuchos::RCP<Xpetra::Import<LocalOrdinal, GlobalOrdinal, Node> >& rowImport,
                 Teuchos::ArrayRCP<LocalOrdinal>& regionMatVecLIDs,
                 Teuchos::RCP<Xpetra::Import<LocalOrdinal, GlobalOrdinal, Node> >& regionInterfaceImporter) {
#include "Xpetra_UseShortNamesOrdinal.hpp"
  using Teuchos::TimeMonitor;

  RCP<TimeMonitor> tm;
  tm = rcp(new TimeMonitor(*TimeMonitor::getNewTimer("SetupMatVec: 1 - sumInterfaceValues")));

  const LO  maxRegPerGID  = static_cast<LO>(regionsPerGIDWithGhosts->getNumVectors());
  const int myRank        = regionRowMap->getComm()->getRank();
  interfaceGIDsMV->replaceMap(regionRowMap);
  RCP<Xpetra::MultiVector<GO, LO, GO, NO> > interfaceGIDs;
  interfaceGIDs = interfaceGIDsMV;
  sumInterfaceValues(interfaceGIDs, regionRowMap, rowImport);

  tm = Teuchos::null;
  tm = rcp(new TimeMonitor(*TimeMonitor::getNewTimer("SetupMatVec: 2 - build regionInterfaceMap")));

  Teuchos::Array<LO> regionMatVecLIDstmp;
  Teuchos::Array<GO> regionMatVecGIDs;
  Array<ArrayRCP<const LO> > regionsPerGIDWithGhostsData(maxRegPerGID);
  Array<ArrayRCP<const GO> > interfaceGIDsData(maxRegPerGID);
  for(LO regionIdx = 0; regionIdx < maxRegPerGID; ++regionIdx) {
    regionsPerGIDWithGhostsData[regionIdx] = regionsPerGIDWithGhosts->getData(regionIdx);
    interfaceGIDsData[regionIdx] = interfaceGIDs->getData(regionIdx);
    for(LO idx = 0; idx < static_cast<LO>(regionsPerGIDWithGhostsData[regionIdx].size()); ++idx) {
      if((regionsPerGIDWithGhostsData[regionIdx][idx] != -1)
         && (regionsPerGIDWithGhostsData[regionIdx][idx] != myRank)) {
        regionMatVecLIDstmp.push_back(idx);
        regionMatVecGIDs.push_back(interfaceGIDsData[regionIdx][idx]);
      }
    }
  }

  // Copy the temporary regionMatVecLIDstmp into an ArrayRCP
  // so we can store it and retrieve it easily later on.
  regionMatVecLIDs.deepCopy(regionMatVecLIDstmp());

  RCP<Map> regionInterfaceMap = Xpetra::MapFactory<LO,GO,Node>::Build(regionRowMap->lib(),
                                                                      Teuchos::OrdinalTraits<GO>::invalid(),
                                                                      regionMatVecGIDs(),
                                                                      regionRowMap->getIndexBase(),
                                                                      regionRowMap->getComm());

  tm = Teuchos::null;
  tm = rcp(new TimeMonitor(*TimeMonitor::getNewTimer("SetupMatVec: 3 - Build importer")));

  regionInterfaceImporter = ImportFactory::Build(regionRowMap, regionInterfaceMap);

  tm = Teuchos::null;
} // SetupMatVec

/*! \brief Compute a matrix vector product \f$Y = beta*Y + alpha*Ax\f$
 *
 *  The residual is computed based on matrices and vectors in a regional layout.
 *  1. Compute y = A*x in regional layout.
 *  2. Sum interface values of y to account for duplication of interface DOFs.
 */
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void ApplyMatVec(const Scalar alpha,
                 const RCP<Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> >& regionMatrix,
                 const RCP<Xpetra::Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node> >& X,
                 const Scalar beta,
                 const RCP<Xpetra::Import<LocalOrdinal, GlobalOrdinal, Node> >& regionInterfaceImporter,
                 const Teuchos::ArrayRCP<LocalOrdinal>& regionInterfaceLIDs,
                 RCP<Xpetra::Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node> >& Y,
                 const Teuchos::ETransp transposeMode,
                 const bool sumInterfaceValues) {
#include "Xpetra_UseShortNames.hpp"
  using Teuchos::TimeMonitor;
  using local_matrix_type = typename Xpetra::Matrix<SC,LO,GO,Node>::local_matrix_type;
  using device_type       = typename local_matrix_type::device_type;

  RCP<TimeMonitor> tm = rcp(new TimeMonitor(*TimeMonitor::getNewTimer("ApplyMatVec: 1 - local apply")));
  RCP<const Map> regionInterfaceMap = regionInterfaceImporter->getTargetMap();

  // Step 1: apply the local operator
  // since in region formate the matrix is block diagonal
  // regionMatrix->apply(*X, *Y, Teuchos::NO_TRANS, alpha, beta);
  local_matrix_type localA = regionMatrix->getLocalMatrix();
  auto localX = X->getDeviceLocalView();
  auto localY = Y->getDeviceLocalView();
  char spmvMode = KokkosSparse::NoTranspose[0];
  if (transposeMode == Teuchos::TRANS)
    spmvMode = KokkosSparse::Transpose[0];
  else
    TEUCHOS_TEST_FOR_EXCEPT_MSG(false, "Unsupported mode.");
  KokkosSparse::spmv(&spmvMode, alpha, localA, localX, beta, localY);

  if (sumInterfaceValues)
  {
    tm = Teuchos::null;
    tm = rcp(new TimeMonitor(*TimeMonitor::getNewTimer("ApplyMatVec: 2 - communicate data")));

    // Step 2: preform communication to propagate local interface
    // values to all the processor that share interfaces.
    RCP<MultiVector> matvecInterfaceTmp = MultiVectorFactory::Build(regionInterfaceMap, 1);
    matvecInterfaceTmp->doImport(*Y, *regionInterfaceImporter, Xpetra::INSERT);

    tm = Teuchos::null;
    tm = rcp(new TimeMonitor(*TimeMonitor::getNewTimer("ApplyMatVec: 3 - sum interface contributions")));

    // Step 3: sum all contributions to interface values
    // on all ranks
    ArrayRCP<Scalar> YData = Y->getDataNonConst(0);
    ArrayRCP<Scalar> interfaceData = matvecInterfaceTmp->getDataNonConst(0);
    for(LO interfaceIdx = 0; interfaceIdx < static_cast<LO>(interfaceData.size()); ++interfaceIdx) {
      YData[regionInterfaceLIDs[interfaceIdx]] += interfaceData[interfaceIdx];
    }
  }

  tm = Teuchos::null;
} // ApplyMatVec


/*! \brief Compute the residual \f$r = b - Ax\f$ with pre-computed communication patterns
 *
 *  The residual is computed based on matrices and vectors in a regional layout.
 *  1. Compute and sum interface values of y = A*x in regional layout using the fast MatVev
 *  2. Compute residual via r = b - y
 */
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void
computeResidual(RCP<Xpetra::Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node> >& regRes, ///< residual (to be evaluated)
                const RCP<Xpetra::Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node> > regX, ///< left-hand side (solution)
                const RCP<Xpetra::Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node> > regB, ///< right-hand side (forcing term)
                const RCP<Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> > regionMats, ///< matrix in region format
                const Teuchos::ParameterList& params ///< parameter with fast MatVec parameters and pre-computed communication patterns
    )
{
#include "Xpetra_UseShortNames.hpp"
  using TST = Teuchos::ScalarTraits<Scalar>;
  using Teuchos::TimeMonitor;

  RCP<TimeMonitor> tm = rcp(new TimeMonitor(*TimeMonitor::getNewTimer("computeResidual: use fast MatVec")));

  // Get pre-communicated communication patterns for the fast MatVec
  const ArrayRCP<LocalOrdinal> regionInterfaceLIDs = params.get<ArrayRCP<LO>>("Fast MatVec: interface LIDs");
  const RCP<Import> regionInterfaceImporter = params.get<RCP<Import>>("Fast MatVec: interface importer");

  // Step 1: Compute region version of y = Ax
  RCP<Vector> aTimesX = VectorFactory::Build(regionMats->getRangeMap(), true);
  ApplyMatVec(TST::one(), regionMats, regX,
      TST::zero(), regionInterfaceImporter, regionInterfaceLIDs, aTimesX, Teuchos::NO_TRANS, true);

  // Step 2: Compute region version of r = b - y
  regRes->update(TST::one(), *regB, -TST::one(), *aTimesX, TST::zero());

  tm = Teuchos::null;
} // computeResidual

#endif // MUELU_SETUPREGIONMATRIX_DEF_HPP
