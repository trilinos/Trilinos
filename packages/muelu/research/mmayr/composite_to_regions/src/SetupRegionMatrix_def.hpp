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


template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class widget>
void MakeGroupRegionRowMaps(const int myRank,
                            const RCP<Xpetra::CrsMatrixWrap<Scalar, LocalOrdinal, GlobalOrdinal, Node> > AComp,
                            LocalOrdinal (*LIDregion)(void*, const LocalOrdinal, int),
                            widget appData,
                            const std::vector<LocalOrdinal> myRegions,
                            std::vector<Teuchos::RCP<Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node> > >& rowMapPerGrp) {
#include "Xpetra_UseShortNames.hpp"
  // Make group region row maps
  // std::cout << myRank << " | Creating region group row maps ..." << std::endl;

  Teuchos::Array<GlobalOrdinal> rowGIDsReg;
  const Teuchos::ArrayView<const GlobalOrdinal> colGIDsComp = AComp->getColMap()->getNodeElementList();

  for (int k = 0; k < static_cast<int>(myRegions.size()); k++) {
    rowGIDsReg.resize(0);
    std::vector<int> tempRegIDs(AComp->getColMap()->getNodeNumElements());
    for (int i = 0; i < static_cast<int>(AComp->getColMap()->getNodeNumElements()); i++) {
      tempRegIDs[i] = (*LIDregion)(&appData, i, k);
    }

    std::vector<int> idx(tempRegIDs.size());
    std::iota(idx.begin(), idx.end(), 0);
    std::sort(idx.begin(), idx.end(), [tempRegIDs](int i1,int i2) { return tempRegIDs[i1] < tempRegIDs[i2];});

    for (int i = 0; i < static_cast<int>(AComp->getColMap()->getNodeNumElements()); i++) {
      if (tempRegIDs[idx[i]] != -1)
        rowGIDsReg.push_back(colGIDsComp[idx[i]]);
    }
    rowMapPerGrp[k] = MapFactory::Build(AComp->getMap()->lib(),
                                        Teuchos::OrdinalTraits<GlobalOrdinal>::invalid(),
                                        rowGIDsReg,
                                        Teuchos::OrdinalTraits<GlobalOrdinal>::zero(),
                                        AComp->getMap()->getComm());
  }

  for (int k = static_cast<int>(myRegions.size()); k < appData.maxRegPerProc; k++) {
    rowMapPerGrp[k] = MapFactory::Build(AComp->getMap()->lib(),
                                        Teuchos::OrdinalTraits<GlobalOrdinal>::invalid(),
                                        Teuchos::OrdinalTraits<GlobalOrdinal>::zero(),
                                        AComp->getMap()->getComm());
  }
}


template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class widget>
void MakeGroupRegionColumnMaps(const int myRank,
                               const int whichCase,
                               const RCP<Xpetra::CrsMatrixWrap<Scalar, LocalOrdinal, GlobalOrdinal, Node> > AComp,
                               LocalOrdinal (*LIDregion)(void*, const LocalOrdinal, int),
                               widget appData,
                               const std::vector<LocalOrdinal> myRegions,
                               std::vector<Teuchos::RCP<Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node> > >& rowMapPerGrp,
                               std::vector<Teuchos::RCP<Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node> > >& colMapPerGrp) {
#include "Xpetra_UseShortNames.hpp"
  // std::cout << myRank << " | Creating region group column maps ..." << std::endl;

  if (whichCase == MultipleRegionsPerProc) {//so maxRegPerProc > 1
    // clone rowMap
    for (int j = 0; j < appData.maxRegPerProc; j++) {
      if (j < static_cast<int>(myRegions.size())) {
        colMapPerGrp[j] = MapFactory::Build(AComp->getMap()->lib(),
                                            Teuchos::OrdinalTraits<int>::invalid(),
                                            rowMapPerGrp[j]->getNodeElementList(),
                                            Teuchos::OrdinalTraits<int>::zero(),
                                            AComp->getMap()->getComm());
      } else {
        colMapPerGrp[j] = MapFactory::Build(AComp->getMap()->lib(),
                                            Teuchos::OrdinalTraits<GlobalOrdinal>::invalid(),
                                            Teuchos::OrdinalTraits<GlobalOrdinal>::zero(),
                                            AComp->getMap()->getComm());
      }
    }
  } else if (whichCase == RegionsSpanProcs) {//so maxRegPerProc = 1
    Teuchos::Array<GlobalOrdinal> colIDsReg;

    // copy the rowmap
    Teuchos::ArrayView<const GlobalOrdinal> rowGIDsReg = rowMapPerGrp[0]->getNodeElementList();
    for (LocalOrdinal i = 0; i < static_cast<LocalOrdinal>(rowMapPerGrp[0]->getNodeNumElements()); i++) {
      colIDsReg.push_back(rowGIDsReg[i]);
    }

    // append additional ghosts who are in my region and
    // for whom I have a LID
    LocalOrdinal LID;
    Teuchos::ArrayView<const GlobalOrdinal> colGIDsComp =  AComp->getColMap()->getNodeElementList();
    for (std::size_t i = 0; i < AComp->getColMap()->getNodeNumElements(); i++) {
      LID = LIDregion(&appData, i, 0);
      if (LID == -1) {
        for (int j = 0; j < appData.maxRegPerGID; j++) {
          Teuchos::ArrayRCP<const LocalOrdinal> jthRegions = appData.regionsPerGIDWithGhosts->getData(j);
          if  (static_cast<int>(jthRegions[i]) == myRegions[0]) {
            colIDsReg.push_back(colGIDsComp[i]);
            break;
          }
        }
      }
    }
    if (static_cast<int>(myRegions.size()) > 0) {
      colMapPerGrp[0] = MapFactory::Build(AComp->getMap()->lib(),
                                          Teuchos::OrdinalTraits<int>::invalid(),
                                          colIDsReg,
                                          Teuchos::OrdinalTraits<int>::zero(),
                                          AComp->getMap()->getComm());
    } else {
      colMapPerGrp[0] = MapFactory::Build(AComp->getMap()->lib(),
                                          Teuchos::OrdinalTraits<GlobalOrdinal>::invalid(),
                                          Teuchos::OrdinalTraits<GlobalOrdinal>::zero(),
                                          AComp->getMap()->getComm());
    }
  } else {
    fprintf(stderr,"whichCase not set properly\n");
    exit(1);
  }
} // MakeGroupRegionColumnMaps

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class widget>
void MakeGroupExtendedMaps(const int myRank,
                           const int whichCase,
                           const RCP<Xpetra::CrsMatrixWrap<Scalar, LocalOrdinal, GlobalOrdinal, Node> > AComp,
                           const RCP<Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node> > mapComp,
                           LocalOrdinal (*LIDregion)(void*, const LocalOrdinal, int),
                           widget appData,
                           const std::vector<LocalOrdinal> myRegions,
                           std::vector<RCP<Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node> > >& rowMapPerGrp,
                           std::vector<RCP<Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node> > >& colMapPerGrp,
                           std::vector<RCP<Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node> > >& revisedRowMapPerGrp,
                           std::vector<RCP<Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node> > >& revisedColMapPerGrp) {
#include "Xpetra_UseShortNames.hpp"
  // std::cout << myRank << " | Creating extended group region maps ..." << std::endl;

  const Scalar SC_ONE = Teuchos::ScalarTraits<Scalar>::one();

  LocalOrdinal nLocal = AComp->getRowMap()->getNodeNumElements();
  LocalOrdinal nExtended = AComp->getColMap()->getNodeNumElements();
  LocalOrdinal nTotal = 0;
  Teuchos::reduceAll<LocalOrdinal, LocalOrdinal>(*(AComp->getMap()->getComm()),
                                                 Teuchos::REDUCE_SUM,
                                                 nLocal,
                                                 Teuchos::outArg(nTotal));

  // first step of NewGID calculation just counts the number of NewGIDs
  // and sets firstNewGID[k] such that it is equal to the number of
  // NewGIDs that have already been counted in (0:k-1) that myRank
  // is responsible for.

  RCP<Vector> firstNewGID = VectorFactory::Build(mapComp, true);
  ArrayRCP<Scalar> firstNewGIDData = firstNewGID->getDataNonConst(0);
  for (int k = 0; k < nLocal-1; k++) {
    firstNewGIDData[k + 1] = firstNewGIDData[k] - 1;
    // firstNewGID->replaceLocalValue(k + 1, (firstNewGID->getData(0))[k]-1);
    for (int j = 0; j < appData.maxRegPerGID; j++) {
      Teuchos::ArrayRCP<const LocalOrdinal> jthRegions = appData.regionsPerGIDWithGhosts->getData(j);
      if (jthRegions[k] != -1) firstNewGIDData[k+1] = SC_ONE;
      // if (jthRegions[k] != -1) firstNewGID->sumIntoLocalValue(k+1, SC_ONE);
    }
  }
  // So firstNewGID[nLocal-1] is number of NewGIDs up to nLocal-2
  // To account for newGIDs associated with nLocal-1, we'll just
  // use an upper bound (to avoid the little loop above).
  // By adding maxRegPerGID-1 we account for the maximum
  // number of possible newGIDs due to last composite id
  GlobalOrdinal upperBndNumNewGIDs = Teuchos::as<GlobalOrdinal>(firstNewGIDData[nLocal-1]) + appData.maxRegPerGID-1;
  GlobalOrdinal upperBndNumNewGIDsAllProcs;
  Teuchos::reduceAll<GlobalOrdinal,GlobalOrdinal>(*(AComp->getMap()->getComm()),
                                                  Teuchos::REDUCE_MAX,
                                                  upperBndNumNewGIDs,
                                                  Teuchos::outArg(upperBndNumNewGIDsAllProcs));

  //    std::cout << "upperBndNumNewGIDsAllProcs: " << upperBndNumNewGIDsAllProcs << std::endl;

  // Now that we have an upper bound on the maximum number of
  // NewGIDs over all procs, we sweep through firstNewGID again
  // to assign ids to the first NewGID associated with each row of
  // regionsPerGIDWithGhosts (by adding an offset)
  for (LocalOrdinal k = 0; k < nLocal; k++) {
    firstNewGIDData[k] += upperBndNumNewGIDsAllProcs*myRank+nTotal;
    // firstNewGID->sumIntoLocalValue(k, upperBndNumNewGIDsAllProcs*myRank+nTotal);
  }

  RCP<Import> Importer = ImportFactory::Build(mapComp, AComp->getColMap());
  RCP<Vector> firstNewGIDWithGhost = VectorFactory::Build(Importer->getTargetMap()); //AComp->getColMap());
  firstNewGIDWithGhost->doImport(*firstNewGID, *Importer, Xpetra::INSERT);

  Teuchos::Array<GlobalOrdinal> revisedGIDs;

  for (std::size_t k = 0; k < myRegions.size(); k++) {
    revisedGIDs.resize(0);
    int curRegion = myRegions[k];
    Teuchos::ArrayView<const GlobalOrdinal> colGIDsComp = AComp->getColMap()->getNodeElementList();
    Teuchos::Array<LocalOrdinal> tempRegIDs(nExtended);

    // must put revisedGIDs in application-provided order by
    // invoking LIDregion() and sorting
    for (LocalOrdinal i = 0; i < nExtended; i++) {
      tempRegIDs[i] = LIDregion(&appData, i, k);
    }
    Teuchos::Array<LocalOrdinal> idx(tempRegIDs.size());
    std::iota(idx.begin(),idx.end(),0);
    std::sort(idx.begin(),idx.end(),[tempRegIDs](int i1,int i2){return tempRegIDs[i1] < tempRegIDs[i2];});

    // Now sweep through regionsPerGIDWithGhosts looking for those
    // associated with curRegion and record the revisedGID
    int j;
    for (LocalOrdinal i = 0; i < nExtended; i++) {

      if (tempRegIDs[idx[i]] != -1) {// if a valid LID for this region
        for (j = 0; j < appData.maxRegPerGID; j++) {
          Teuchos::ArrayRCP<const LocalOrdinal> jthRegions = appData.regionsPerGIDWithGhosts->getData(j);
          if (jthRegions[idx[i]] == curRegion) break;
        }

        // (regionsPerGIDWithGhosts)[0] entries keep original GID
        // while others use firstNewGID to determine NewGID

        if (j == 0) revisedGIDs.push_back(colGIDsComp[idx[i]]);
        else if (j < appData.maxRegPerGID) {
          revisedGIDs.push_back((GlobalOrdinal) firstNewGIDWithGhost->getData(0)[idx[i]] + j - 1);
        }

        // add entry to listDulicatedGIDs
      }
    }

    revisedRowMapPerGrp[k] = MapFactory::Build(AComp->getMap()->lib(),
                                               Teuchos::OrdinalTraits<GlobalOrdinal>::invalid(),
                                               revisedGIDs,
                                               Teuchos::OrdinalTraits<GlobalOrdinal>::zero(),
                                               AComp->getMap()->getComm());

    // now append more stuff to handle ghosts ... needed for
    // revised version of column map

    for (LocalOrdinal i = 0; i < nExtended; i++) {
      if (tempRegIDs[i] == -1) {// only in revised col map
        // note: since sorting not used when making the
        // original regional colmap, we can't use
        // it here either ... so no idx[]'s.
        for (j = 0; j < appData.maxRegPerGID; j++) {
          Teuchos::ArrayRCP<const LocalOrdinal> jthRegions = appData.regionsPerGIDWithGhosts->getData(j);
          if  (jthRegions[i] == curRegion) break;
        }
        // (*regionsPerGIDWithGhosts)[0] entries keep original GID
        // while others use firstNewGID to determine NewGID

        if (j == 0) revisedGIDs.push_back(colGIDsComp[i]);
        else if (j < appData.maxRegPerGID) {
          revisedGIDs.push_back((int) firstNewGIDWithGhost->getData(0)[i] + j - 1);
        }
      }
    }
    revisedColMapPerGrp[k] = MapFactory::Build(AComp->getMap()->lib(),
                                               Teuchos::OrdinalTraits<GlobalOrdinal>::invalid(),
                                               revisedGIDs,
                                               Teuchos::OrdinalTraits<GlobalOrdinal>::zero(),
                                               AComp->getMap()->getComm());
  }
  for (std::size_t k = myRegions.size(); k < Teuchos::as<std::size_t>(appData.maxRegPerProc); k++) {
    revisedRowMapPerGrp[k] = MapFactory::Build(AComp->getMap()->lib(),
                                               Teuchos::OrdinalTraits<GlobalOrdinal>::invalid(),
                                               Teuchos::OrdinalTraits<GlobalOrdinal>::zero(),
                                               AComp->getMap()->getComm());
    revisedColMapPerGrp[k] = MapFactory::Build(AComp->getMap()->lib(),
                                               Teuchos::OrdinalTraits<GlobalOrdinal>::invalid(),
                                               Teuchos::OrdinalTraits<GlobalOrdinal>::zero(),
                                               AComp->getMap()->getComm());
  }
} // MakeGroupExtendedMaps

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void MakeQuasiregionMatrices(const RCP<Xpetra::CrsMatrixWrap<Scalar, LocalOrdinal, GlobalOrdinal, Node> > AComp,
                             const int maxRegPerProc,
                             RCP<Xpetra::MultiVector<LocalOrdinal, LocalOrdinal, GlobalOrdinal, Node> > regionsPerGIDWithGhosts,
                             std::vector<RCP<Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node> > >& rowMapPerGrp,
                             std::vector<RCP<Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node> > >& colMapPerGrp,
                             std::vector<RCP<Xpetra::Import<LocalOrdinal, GlobalOrdinal, Node> > >& rowImportPerGrp,
                             std::vector<RCP<Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> > >& quasiRegionGrpMats) {
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
  for (int j = 0; j < maxRegPerProc; j++) {
    RCP<TimeMonitor> tm = rcp(new TimeMonitor(*TimeMonitor::getNewTimer("MakeQuasiregionMatrices: 1 - Create Matrix")));

    quasiRegionGrpMats[j] = MatrixFactory::Build(rowMapPerGrp[j],
                                                 colMapPerGrp[j],
                                                 AComp->getCrsGraph()->getNodeMaxNumRowEntries());

    tm = Teuchos::null;
    tm = rcp(new TimeMonitor(*TimeMonitor::getNewTimer("MakeQuasiregionMatrices: 2 - Import data")));

    quasiRegionGrpMats[j]->doImport(*AComp,
                                    *(rowImportPerGrp[j]),
                                    Xpetra::INSERT);

    tm = Teuchos::null;
    tm = rcp(new TimeMonitor(*TimeMonitor::getNewTimer("MakeQuasiregionMatrices: 3 - Scale interface entries")));

    RCP<CrsMatrixWrap> quasiRegionCrsWrap = Teuchos::rcp_dynamic_cast<CrsMatrixWrap>(quasiRegionGrpMats[j]);
    RCP<CrsMatrix> quasiRegionCrs = quasiRegionCrsWrap->getCrsMatrix();
    const LO numRows = Teuchos::as<LO>(quasiRegionCrs->getNodeNumRows());

    for (LO row = 0; row < numRows; ++row) { // loop over local rows of composite matrix
      GO rowGID = rowMapPerGrp[j]->getGlobalElement(row);
      std::size_t numEntries = quasiRegionGrpMats[j]->getNumEntriesInLocalRow(row); // number of entries in this row
      Teuchos::Array<SC> vals(numEntries); // non-zeros in this row
      Teuchos::Array<LO> inds(numEntries); // local column indices
      quasiRegionGrpMats[j]->getLocalRowCopy(row, inds, vals, numEntries);

      for (std::size_t c = 0; c < Teuchos::as<std::size_t>(inds.size()); ++c) { // loop over all entries in this row
        LocalOrdinal col = inds[c];
        GlobalOrdinal colGID = colMapPerGrp[j]->getGlobalElement(col);
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

      quasiRegionGrpMats[j]->replaceLocalValues(row, inds, vals);
    }

    tm = Teuchos::null;
    tm = rcp(new TimeMonitor(*TimeMonitor::getNewTimer("MakeQuasiregionMatrices: 4 - fillComplete")));

    quasiRegionGrpMats[j]->fillComplete();

    tm = Teuchos::null;
  }
} // MakeQuasiregionMatrices

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void MakeRegionMatrices(const RCP<Xpetra::CrsMatrixWrap<Scalar, LocalOrdinal, GlobalOrdinal, Node> > AComp,
                        const RCP<const Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node> > mapComp,
                        std::vector<RCP<Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node> > >& rowMapPerGrp,
                        std::vector<RCP<Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node> > >& revisedRowMapPerGrp,
                        std::vector<RCP<Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node> > >& revisedColMapPerGrp,
                        std::vector<RCP<Xpetra::Import<LocalOrdinal, GlobalOrdinal, Node> > >& rowImportPerGrp,
                        const LocalOrdinal maxRegPerProc,
                        std::vector<RCP<Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> > >& quasiRegionGrpMats,
                        std::vector<RCP<Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> > >& regionGrpMats) {
#include "Xpetra_UseShortNames.hpp"
  using Teuchos::RCP;
  using Teuchos::TimeMonitor;
  const SC SC_ONE  = Teuchos::ScalarTraits<SC>::one();
  const SC SC_ZERO = Teuchos::ScalarTraits<SC>::zero();

  RCP<TimeMonitor> tm = rcp(new TimeMonitor(*TimeMonitor::getNewTimer("MakeRegionMatrices: 1 - Create Matrix")));

  // Copy data from quasiRegionGrpMats, but into new map layout
  {
    for (int j = 0; j < maxRegPerProc; j++) {
      regionGrpMats[j] = rcp(new CrsMatrixWrap(revisedRowMapPerGrp[j], revisedColMapPerGrp[j], 9));

      // Extract current region CrsMatrix
      RCP<CrsMatrix> regionCrsMat = Teuchos::rcp_dynamic_cast<CrsMatrixWrap>(regionGrpMats[j])->getCrsMatrix();

      // Extract current quasi-region CrsMatrix
      RCP<CrsMatrix> quasiRegionCrsMat = Teuchos::rcp_dynamic_cast<CrsMatrixWrap>(quasiRegionGrpMats[j])->getCrsMatrix();

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
      regionCrsMat->expertStaticFillComplete(revisedRowMapPerGrp[j], revisedRowMapPerGrp[j]);
    }
  }

  tm = Teuchos::null;
  tm = rcp(new TimeMonitor(*TimeMonitor::getNewTimer("MakeRegionMatrices: 2 - Enforce nullspace constraint")));

  // enforce nullspace constraint
  Array<RCP<Vector> > regNspViolation(maxRegPerProc);
  {
    // compute violation of nullspace property close to DBCs
    RCP<Vector> nspVec = VectorFactory::Build(mapComp);
    nspVec->putScalar(SC_ONE);
    RCP<Vector> nspViolation = VectorFactory::Build(mapComp, true);
    AComp->apply(*nspVec, *nspViolation);

    // move to regional layout
    Array<RCP<Vector> > quasiRegNspViolation(maxRegPerProc);
    createRegionalVector(quasiRegNspViolation, rowMapPerGrp);
    createRegionalVector(regNspViolation, revisedRowMapPerGrp);
    compositeToRegional(nspViolation, quasiRegNspViolation, regNspViolation,
                        revisedRowMapPerGrp, rowImportPerGrp);

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
      Array<RCP<Vector> > interfaceScaling(maxRegPerProc);
      for (int j = 0; j < maxRegPerProc; j++) {
        interfaceScaling[j] = VectorFactory::Build(revisedRowMapPerGrp[j]);
        interfaceScaling[j]->putScalar(SC_ONE);
      }

      // transform to composite layout while adding interface values via the Export() combine mode
      RCP<Vector> compInterfaceScalingSum = VectorFactory::Build(mapComp, true);
      regionalToComposite(interfaceScaling, compInterfaceScalingSum, rowImportPerGrp);

      /* transform composite layout back to regional layout. Now, GIDs associated
       * with region interface should carry a scaling factor (!= 1).
       */
      Array<RCP<Vector> > quasiRegInterfaceScaling(maxRegPerProc);
      compositeToRegional(compInterfaceScalingSum, quasiRegInterfaceScaling,
                          interfaceScaling,
                          revisedRowMapPerGrp, rowImportPerGrp);

      // modify its interface entries
      for (int j = 0; j < maxRegPerProc; j++) {
        RCP<Vector> inverseInterfaceScaling = VectorFactory::Build(interfaceScaling[j]->getMap(), true);
        inverseInterfaceScaling->reciprocal(*interfaceScaling[j]);
        regNspViolation[j]->elementWiseMultiply(SC_ONE, *regNspViolation[j], *inverseInterfaceScaling, SC_ZERO);
      }
    }
  }

  tm = Teuchos::null;
  tm = rcp(new TimeMonitor(*TimeMonitor::getNewTimer("MakeRegionMatrices: 3 - Replace diagonal")));

  Array<RCP<Vector> > regNsp(maxRegPerProc);
  Array<RCP<Vector> > regCorrection(maxRegPerProc);
  for (int j = 0; j < maxRegPerProc; j++) {
    regNsp[j] = VectorFactory::Build(revisedRowMapPerGrp[j]);
    regNsp[j]->putScalar(SC_ONE);

    regCorrection[j] = VectorFactory::Build(revisedRowMapPerGrp[j], true);
    regionGrpMats[j]->apply(*regNsp[j], *regCorrection[j]);
  }

  RCP<Vector> regDiag = Teuchos::null;
  for (int j = 0; j < maxRegPerProc; j++) {
    regDiag = VectorFactory::Build(revisedRowMapPerGrp[j], true);
    regionGrpMats[j]->getLocalDiagCopy(*regDiag);
    regDiag->update(-SC_ONE, *regCorrection[j], SC_ONE, *regNspViolation[j], SC_ONE);

    // Extract current region matrix in as CrsMatrix
    RCP<CrsMatrix> regionCrsMat = Teuchos::rcp_dynamic_cast<CrsMatrixWrap>(regionGrpMats[j])->getCrsMatrix();
    regionCrsMat->replaceDiag(*regDiag);
  }

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
void regionalToComposite(const std::vector<RCP<Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> > >& regMat, ///< Matrix in region layout [in]
                         const std::vector<RCP<Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node> > > rowMapPerGrp, ///< row maps in quasiRegion layout [in]
                         const std::vector<RCP<Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node> > > colMapPerGrp, ///< col maps in quasiRegion layout [in]
                         const std::vector<RCP<Xpetra::Import<LocalOrdinal, GlobalOrdinal, Node> > > rowImportPerGrp, ///< row importer in region layout [in]
                         const Xpetra::CombineMode combineMode, ///< Combine mode for import/export [in]
                         RCP<Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> >& compMat ///< Matrix in composite layout [in/out]
                         )
{
#include "Xpetra_UseShortNames.hpp"
  using Teuchos::TimeMonitor;
  using Teuchos::rcp;
  using std::size_t;

  // Get max number of regions per proc
  const int maxRegPerProc = regMat.size();

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
  std::vector<RCP<Matrix> > quasiRegMat(maxRegPerProc);
  {
    for (int j = 0; j < maxRegPerProc; j++) {
      quasiRegMat[j] = rcp(new CrsMatrixWrap(rowMapPerGrp[j],
                                             colMapPerGrp[j],
                                             regMat[j]->getCrsGraph()->getNodeMaxNumRowEntries()));

      // Extract current quasi-region CrsMatrix
      RCP<CrsMatrix> quasiRegionCrsMat = Teuchos::rcp_dynamic_cast<CrsMatrixWrap>(quasiRegMat[j])->getCrsMatrix();

      // Extract current region CrsMatrix
      RCP<CrsMatrix> regionCrsMat = Teuchos::rcp_dynamic_cast<CrsMatrixWrap>(regMat[j])->getCrsMatrix();

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
      quasiRegionCrsMat->expertStaticFillComplete(rowMapPerGrp[j], rowMapPerGrp[j]);
    }
  }

  // Export from quasiRegional format to composite layout
  std::vector<RCP<Matrix> > partialCompMat(maxRegPerProc);
  for (int j = 0; j < maxRegPerProc; j++) {
    partialCompMat[j] = MatrixFactory::Build(compMat->getRowMap(),
                                             8*regMat[0]->getCrsGraph()->getNodeMaxNumRowEntries());
    partialCompMat[j]->doExport(*(quasiRegMat[j]), *(rowImportPerGrp[j]), Xpetra::INSERT);
    partialCompMat[j]->fillComplete();
  }

  // Add all partialCompMat together
  for (int j = 0; j < maxRegPerProc; j++) {
    MatrixMatrix::TwoMatrixAdd(*partialCompMat[j], false, SC_ONE, *compMat, SC_ONE);
  }

  compMat->fillComplete();

  return;
} // regionalToComposite

/*! \brief Compute the residual \f$r = b - Ax\f$ using two rounds of communication
 *
 *  The residual is computed based on matrices and vectors in a regional layout.
 *  1. Compute y = A*x in regional layout.
 *  2. Sum interface values of y to account for duplication of interface DOFs.
 *  3. Compute r = b - y
 */
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
Array<RCP<Xpetra::Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node> > >
computeResidual(Array<RCP<Xpetra::Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node> > >& regRes, ///< residual (to be evaluated)
                const Array<RCP<Xpetra::Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node> > > regX, ///< left-hand side (solution)
                const Array<RCP<Xpetra::Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node> > > regB, ///< right-hand side (forcing term)
                const std::vector<RCP<Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> > > regionGrpMats,
                const std::vector<RCP<Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node> > > revisedRowMapPerGrp, ///< revised row maps in region layout [in] (actually extracted from regionGrpMats)
                const std::vector<RCP<Xpetra::Import<LocalOrdinal, GlobalOrdinal, Node> > > rowImportPerGrp ///< row importer in region layout [in]
    )
{
#include "Xpetra_UseShortNames.hpp"
  using Teuchos::TimeMonitor;
  const int maxRegPerProc = regX.size();

  RCP<TimeMonitor> tm = rcp(new TimeMonitor(*TimeMonitor::getNewTimer("computeResidual: 1 - Rreg = Areg*Xreg")));

  /* Update the residual vector
   * 1. Compute tmp = A * regX in each region
   * 2. Sum interface values in tmp due to duplication (We fake this by scaling to reverse the basic splitting)
   * 3. Compute r = B - tmp
   */
  for (int j = 0; j < maxRegPerProc; j++) { // step 1
    regionGrpMats[j]->apply(*regX[j], *regRes[j]);
    //    TEUCHOS_ASSERT(regionGrpMats[j]->getDomainMap()->isSameAs(*regX[j]->getMap()));
    //    TEUCHOS_ASSERT(regionGrpMats[j]->getRangeMap()->isSameAs(*regRes[j]->getMap()));
  }

  tm = Teuchos::null;
  tm = rcp(new TimeMonitor(*TimeMonitor::getNewTimer("computeResidual: 2 - sumInterfaceValues")));

  sumInterfaceValues(regRes, revisedRowMapPerGrp, rowImportPerGrp);

  tm = Teuchos::null;
  tm = rcp(new TimeMonitor(*TimeMonitor::getNewTimer("computeResidual: 3 - Rreg = Breg - Rreg")));

  for (int j = 0; j < maxRegPerProc; j++) { // step 3
    regRes[j]->update(1.0, *regB[j], -1.0);
    //    TEUCHOS_ASSERT(regRes[j]->getMap()->isSameAs(*regB[j]->getMap()));
  }

  tm = Teuchos::null;

  return regRes;
} // computeResidual

/*! \brief Compute local data needed to perform a MatVec in region format

 */
template <class LocalOrdinal, class GlobalOrdinal, class Node>
void SetupMatVec(const Teuchos::RCP<Xpetra::MultiVector<GlobalOrdinal, LocalOrdinal, GlobalOrdinal, Node> >& interfaceGIDsMV,
                 const Teuchos::RCP<Xpetra::MultiVector<LocalOrdinal, LocalOrdinal, GlobalOrdinal, Node> >& regionsPerGIDWithGhosts,
                 const std::vector<Teuchos::RCP<Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node> > >& regionRowMap,
                 const std::vector<Teuchos::RCP<Xpetra::Import<LocalOrdinal, GlobalOrdinal, Node> > >& rowImportPerGrp,
                 Teuchos::ArrayRCP<LocalOrdinal>& regionMatVecLIDs,
                 Teuchos::RCP<Xpetra::Import<LocalOrdinal, GlobalOrdinal, Node> >& regionInterfaceImporter) {
#include "Xpetra_UseShortNamesOrdinal.hpp"
  using Teuchos::TimeMonitor;

  RCP<TimeMonitor> tm;
  tm = rcp(new TimeMonitor(*TimeMonitor::getNewTimer("SetupMatVec: 1 - sumInterfaceValues")));

  const LO  maxRegPerProc = static_cast<LO>(regionRowMap.size());
  const LO  maxRegPerGID  = static_cast<LO>(regionsPerGIDWithGhosts->getNumVectors());
  const int myRank        = regionRowMap[0]->getComm()->getRank();
  interfaceGIDsMV->replaceMap(regionRowMap[0]);
  Array<RCP<Xpetra::MultiVector<GO, LO, GO, NO> > > interfaceGIDsPerGrp(maxRegPerProc);
  interfaceGIDsPerGrp[0] = interfaceGIDsMV;
  sumInterfaceValues(interfaceGIDsPerGrp, regionRowMap, rowImportPerGrp);

  tm = Teuchos::null;
  tm = rcp(new TimeMonitor(*TimeMonitor::getNewTimer("SetupMatVec: 2 - build regionInterfaceMap")));

  Teuchos::Array<LO> regionMatVecLIDstmp;
  Teuchos::Array<GO> regionMatVecGIDs;
  Array<ArrayRCP<const LO> > regionsPerGIDWithGhostsData(maxRegPerGID);
  Array<ArrayRCP<const GO> > interfaceGIDsData(maxRegPerGID);
  for(LO regionIdx = 0; regionIdx < maxRegPerGID; ++regionIdx) {
    regionsPerGIDWithGhostsData[regionIdx] = regionsPerGIDWithGhosts->getData(regionIdx);
    interfaceGIDsData[regionIdx] = interfaceGIDsPerGrp[0]->getData(regionIdx);
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

  RCP<Map> regionInterfaceMap = Xpetra::MapFactory<LO,GO,Node>::Build(regionRowMap[0]->lib(),
                                                                      Teuchos::OrdinalTraits<GO>::invalid(),
                                                                      regionMatVecGIDs(),
                                                                      regionRowMap[0]->getIndexBase(),
                                                                      regionRowMap[0]->getComm());

  tm = Teuchos::null;
  tm = rcp(new TimeMonitor(*TimeMonitor::getNewTimer("SetupMatVec: 3 - Build importer")));

  regionInterfaceImporter = ImportFactory::Build(regionRowMap[0], regionInterfaceMap);

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
  auto localX = X->template getLocalView<device_type>();
  auto localY = Y->template getLocalView<device_type>();
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
computeResidual(Array<RCP<Xpetra::Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node> > >& regRes, ///< residual (to be evaluated)
                const Array<RCP<Xpetra::Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node> > > regX, ///< left-hand side (solution)
                const Array<RCP<Xpetra::Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node> > > regB, ///< right-hand side (forcing term)
                const std::vector<RCP<Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> > > regionGrpMats, ///< matrix in region format
                const Teuchos::ParameterList& params ///< parameter with fast MatVec parameters and pre-computed communication patterns
    )
{
#include "Xpetra_UseShortNames.hpp"
  using TST = Teuchos::ScalarTraits<Scalar>;
  using Teuchos::TimeMonitor;

  TEUCHOS_TEST_FOR_EXCEPT_MSG(regX.size()>1, "Fast MatVec only implemented for one region per group.");

  RCP<TimeMonitor> tm = rcp(new TimeMonitor(*TimeMonitor::getNewTimer("computeResidual: use fast MatVec")));

  // Get pre-communicated communication patterns for the fast MatVec
  const ArrayRCP<LocalOrdinal> regionInterfaceLIDs = params.get<ArrayRCP<LO>>("Fast MatVec: interface LIDs");
  const RCP<Import> regionInterfaceImporter = params.get<RCP<Import>>("Fast MatVec: interface importer");

  // Step 1: Compute region version of y = Ax
  RCP<Vector> aTimesX = VectorFactory::Build(regionGrpMats[0]->getRangeMap(), true);
  ApplyMatVec(TST::one(), regionGrpMats[0], regX[0],
      TST::zero(), regionInterfaceImporter, regionInterfaceLIDs, aTimesX, Teuchos::NO_TRANS, true);

  // Step 2: Compute region version of r = b - y
  regRes[0]->update(TST::one(), *regB[0], -TST::one(), *aTimesX, TST::zero());

  tm = Teuchos::null;
} // computeResidual

#endif // MUELU_SETUPREGIONMATRIX_DEF_HPP
