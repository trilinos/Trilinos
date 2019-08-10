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
#include <numeric>

#define RegionsSpanProcs  1
#define MultipleRegionsPerProc  2
#ifdef HAVE_MPI
#include "mpi.h"
#endif

#include <Kokkos_DefaultNode.hpp>

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
                                      const Xpetra::MultiVector<LocalOrdinal, LocalOrdinal, GlobalOrdinal, Node>& nodesToRegions ///< mapping of nodes to regions
                                      )
{
#include "Xpetra_UseShortNamesOrdinal.hpp"

  // extract node-to-regions mapping for both nodes A and B
  Array<int> regionsA;
  Array<int> regionsB;
  {
    RCP<const Map> map = nodesToRegions.getMap();
    for (std::size_t i = 0; i < nodesToRegions.getNumVectors(); ++i) {
      regionsA.push_back(nodesToRegions.getData(i)[map->getLocalElement(nodeA)]);
      regionsB.push_back(nodesToRegions.getData(i)[map->getLocalElement(nodeB)]);
    }
  }

  // identify common regions
  std::vector<int> commonRegions(nodesToRegions.getNumVectors());
  std::sort(regionsA.begin(), regionsA.end());
  std::sort(regionsB.begin(), regionsB.end());

  std::vector<int>::iterator it = std::set_intersection(regionsA.begin(),
      regionsA.end(), regionsB.begin(), regionsB.end(), commonRegions.begin());
  commonRegions.resize(it - commonRegions.begin());

  // remove '-1' entries
  Teuchos::Array<int> finalCommonRegions;
  for (std::size_t i = 0; i < commonRegions.size(); ++i) {
    if (commonRegions[i] != -1)
      finalCommonRegions.push_back(commonRegions[i]);
  }

  return finalCommonRegions;
}

//! Create an empty vector in the regional layout
template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void createRegionalVector(std::vector<RCP<Xpetra::Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node> > >& regVecs, ///< regional vector to be filled
                          const int maxRegPerProc, ///< max number of regions per process
                          const std::vector<RCP<Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node> > > revisedRowMapPerGrp ///< regional map
                          )
{
#include "Xpetra_UseShortNames.hpp"
  regVecs.resize(maxRegPerProc);
  for (int j = 0; j < maxRegPerProc; j++)
    regVecs[j] = VectorFactory::Build(revisedRowMapPerGrp[j], true);

  return;
} // createRegionalVector

/*! \brief Transform composite vector to regional layout
 *
 *  Starting from a vector in composite layout, we
 *  1. import it into an auxiliary vector in the quasiRegional layout
 *  2. replace the quasiRegional map of the auxiliary vector with the regional map
 */
template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void compositeToRegional(RCP<Xpetra::Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node> > compVec, ///< Vector in composite layout [in]
                         std::vector<RCP<Xpetra::Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node> > >& quasiRegVecs, ///< Vector in quasiRegional layout [in/out]
                         std::vector<RCP<Xpetra::Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node> > >& regVecs, ///< Vector in regional layout [in/out]
                         const int maxRegPerProc, ///< max number of regions per proc [in]
                         const std::vector<RCP<Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node> > > rowMapPerGrp, ///< row maps in region layout [in]
                         const std::vector<RCP<Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node> > > revisedRowMapPerGrp, ///< revised row maps in region layout [in]
                         const std::vector<RCP<Xpetra::Import<LocalOrdinal, GlobalOrdinal, Node> > > rowImportPerGrp ///< row importer in region layout [in]
                         )
{
#include "Xpetra_UseShortNames.hpp"
  // quasiRegional layout
  for (int j = 0; j < maxRegPerProc; j++) {
    // create empty vectors and fill it by extracting data from composite vector
    quasiRegVecs[j] = VectorFactory::Build(rowMapPerGrp[j], true);
    TEUCHOS_ASSERT(!quasiRegVecs[j].is_null());
    quasiRegVecs[j]->doImport(*compVec, *(rowImportPerGrp[j]), Xpetra::INSERT);
  }

  // regional layout
  for (int j = 0; j < maxRegPerProc; j++) {
    // create regVecs vector (copy from quasiRegVecs and swap the map)
    regVecs[j] = quasiRegVecs[j]; // assignment operator= does deep copy in Xpetra
    TEUCHOS_ASSERT(!regVecs[j].is_null());
    regVecs[j]->replaceMap(revisedRowMapPerGrp[j]);
  }

  return;
} // compositeToRegional

/*! \brief Transform regional vector to composite layout
 *
 *  Starting from a \c Epetra_Vector in regional layout, we
 *  1. replace the regional map with the quasiRegional map
 *  2. export it into a vector with composite layout using the \c CombineMode \c Add.
 *     Note: on-process values also have to be added to account for region interfaces inside a process.
 *
 *  \note We also need the capability to add processor-local values. This is not supported by
 *  available CombineMode options in Xpetra/Tpetra, so we use a manual implementation here.
 */
template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void regionalToComposite(const std::vector<RCP<Xpetra::Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node> > >& regVec, ///< Vector in region layout [in]
                         RCP<Xpetra::Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node> > compVec, ///< Vector in composite layout [in/out]
                         const int maxRegPerProc, ///< max number of regions per proc
                         const std::vector<RCP<Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node> > > rowMapPerGrp, ///< row maps in quasiRegion layout [in]
                         const std::vector<RCP<Xpetra::Import<LocalOrdinal, GlobalOrdinal, Node> > > rowImportPerGrp, ///< row importer in region layout [in]
                         const Xpetra::CombineMode combineMode ///< Combine mode for import/export [in]
                         )
{
  /* Let's fake an ADD combine mode that also adds local values by
   * 1. exporting quasiRegional vectors to auxiliary composite vectors (1 per group)
   * 2. add all auxiliary vectors together
   */
#include "Xpetra_UseShortNames.hpp"

  const SC SC_ZERO = Teuchos::ScalarTraits<Scalar>::zero();

  const size_t compVecLocalLength = compVec->getLocalLength();
  compVec->putScalar(SC_ZERO);

  {
    RCP<Vector> quasiRegVec;
    for(int grpIdx = 0; grpIdx < maxRegPerProc; ++grpIdx) {
      quasiRegVec = regVec[grpIdx];
      TEUCHOS_ASSERT(Teuchos::nonnull(quasiRegVec));
      quasiRegVec->replaceMap(rowImportPerGrp[grpIdx]->getTargetMap());

      RCP<Vector> partialCompVec = VectorFactory::Build(rowImportPerGrp[0]->getSourceMap(), true);
      TEUCHOS_ASSERT(Teuchos::nonnull(partialCompVec));
      TEUCHOS_ASSERT(partialCompVec->getLocalLength() == compVecLocalLength);
      // ToDo (mayr.mt) Use input variable 'combineMode'
      partialCompVec->doExport(*quasiRegVec, *(rowImportPerGrp[grpIdx]), Xpetra::ADD);

      Teuchos::ArrayRCP<const SC> partialCompVecData = partialCompVec->getData(0);
      for(size_t entryIdx = 0; entryIdx < compVecLocalLength; ++entryIdx) {
        compVec->sumIntoLocalValue(static_cast<LO>(entryIdx), partialCompVecData[entryIdx]);
      }
    }
  }

  return;
} // regionalToComposite

/*! \brief Sum region interface values
 *
 *  Sum values of interface GIDs using the underlying Export() routines. Technically, we perform the
 *  exchange/summation of interface data by exporting a regional vector to the composite layout and
 *  then immediately importing it back to the regional layout. The Export() involved when going to the
 *  composite layout takes care of the summation of interface values.
 */
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void sumInterfaceValues(std::vector<RCP<Xpetra::Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node> > >& regVec,
                        const RCP<Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node> > compMap,
                        const int maxRegPerProc, ///< max number of regions per proc [in]
                        const std::vector<RCP<Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node> > > rowMapPerGrp,///< row maps in region layout [in]
                        const std::vector<RCP<Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node> > > revisedRowMapPerGrp,///< revised row maps in region layout [in]
                        const std::vector<RCP<Xpetra::Import<LocalOrdinal, GlobalOrdinal, Node> > > rowImportPerGrp ///< row importer in region layout [in])
                        )
{
#include "Xpetra_UseShortNames.hpp"
  Teuchos::RCP<Vector> compVec = VectorFactory::Build(compMap, true);
  TEUCHOS_ASSERT(!compVec.is_null());

  std::vector<Teuchos::RCP<Vector> > quasiRegVec(maxRegPerProc);
  regionalToComposite(regVec, compVec, maxRegPerProc, rowMapPerGrp,
                      rowImportPerGrp, Xpetra::ADD);

  compositeToRegional(compVec, quasiRegVec, regVec, maxRegPerProc,
                      rowMapPerGrp, revisedRowMapPerGrp, rowImportPerGrp);

  return;
} // sumInterfaceValues

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class widget>
void MakeGroupRegionRowMaps(const int myRank,
                            const RCP<Xpetra::CrsMatrixWrap<Scalar, LocalOrdinal, GlobalOrdinal, Node> > AComp,
                            LocalOrdinal (*LIDregion)(void*, const LocalOrdinal, int),
                            widget appData,
                            const std::vector<LocalOrdinal> myRegions,
                            std::vector<Teuchos::RCP<Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node> > >& rowMapPerGrp) {
#include "Xpetra_UseShortNames.hpp"
  // Make group region row maps
  std::cout << myRank << " | Creating region group row maps ..." << std::endl;

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
  std::cout << myRank << " | Creating region group column maps ..." << std::endl;

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
}

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
  std::cout << myRank << " | Creating extended group region maps ..." << std::endl;

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
  for (int k = 0; k < nLocal-1; k++) {
    firstNewGID->replaceLocalValue(k + 1, (firstNewGID->getData(0))[k]-1);
    for (int j = 0; j < appData.maxRegPerGID; j++) {
      Teuchos::ArrayRCP<const LocalOrdinal> jthRegions = appData.regionsPerGIDWithGhosts->getData(j);
      if (jthRegions[k] != -1) firstNewGID->sumIntoLocalValue(k+1, SC_ONE);
    }
  }
  // So firstNewGID[nLocal-1] is number of NewGIDs up to nLocal-2
  // To account for newGIDs associated with nLocal-1, we'll just
  // use an upper bound (to avoid the little loop above).
  // By adding maxRegPerGID-1 we account for the maximum
  // number of possible newGIDs due to last composite id
  Teuchos::ArrayRCP<const Scalar> firstNewGIDData = firstNewGID->getData(0);
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
  ArrayRCP<Scalar> firstNewGIDDataNonConst = firstNewGID->getDataNonConst(0);
  for (LocalOrdinal k = 0; k < nLocal; k++) {
    //      const GlobalOrdinal tmpGID = firstNewGIDDataNonConst[k];
    //      firstNewGIDDataNonConst[k] = tmpGID + upperBndNumNewGIDsAllProcs*myRank+nTotal;
    firstNewGID->sumIntoLocalValue(k, upperBndNumNewGIDsAllProcs*myRank+nTotal);
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
}

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

  std::cout << AComp->getMap()->getComm()->getRank() << " | Forming quasiRegion matrices ..." << std::endl;

  /* We use the edge-based splitting, i.e. we first modify off-diagonal
   * entries in the composite matrix, then decompose it into region matrices
   * and finally take care of diagonal entries by enforcing the nullspace
   * preservation constraint.
   */


  // // copy and modify the composite matrix
  // RCP<Matrix> ACompSplit = MatrixFactory::BuildCopy(AComp);
  // ACompSplit->resumeFill();

  // for (LocalOrdinal row = 0; row < Teuchos::as<LocalOrdinal>(ACompSplit->getNodeNumRows()); row++) { // loop over local rows of composite matrix
  //   GlobalOrdinal rowGID = ACompSplit->getRowMap()->getGlobalElement(row);
  //   std::size_t numEntries = ACompSplit->getNumEntriesInLocalRow(row); // number of entries in this row
  //   Teuchos::Array<Scalar> vals(numEntries); // non-zeros in this row
  //   Teuchos::Array<LocalOrdinal> inds(numEntries); // local column indices
  //   ACompSplit->getLocalRowCopy(row, inds, vals, numEntries);

  //   for (std::size_t c = 0; c < Teuchos::as<std::size_t>(inds.size()); ++c) { // loop over all entries in this row
  //     LocalOrdinal col = inds[c];
  //     GlobalOrdinal colGID = ACompSplit->getColMap()->getGlobalElement(col);
  //     Array<int> commonRegions;
  //     if (rowGID != colGID) { // Skip the diagonal entry. It will be processed later.
  //       commonRegions = findCommonRegions(rowGID, colGID, *regionsPerGIDWithGhosts);
  //     }

  //     if (commonRegions.size() > 1) {
  //       vals[c] *= (1.0 / static_cast<double>(commonRegions.size()));
  //     }
  //   }

  //   ACompSplit->replaceLocalValues(row, inds, vals);
  // }
  // ACompSplit->fillComplete(AComp->getDomainMap(), AComp->getRangeMap());

  // Import data from ACompSplit into the quasiRegion matrices
  // Since the stencil size does not grow between composite and region format
  // use ACompSplit->getCrsGraph()->getNodeMaxNumRowEntries() to get an upper
  // bound of the number of nonzeros per row in quasiRegionMat
  for (int j = 0; j < maxRegPerProc; j++) {
    quasiRegionGrpMats[j] = MatrixFactory::Build(rowMapPerGrp[j],
                                                 colMapPerGrp[j],
                                                 AComp->getCrsGraph()->getNodeMaxNumRowEntries(),
                                                 Xpetra::DynamicProfile);
    quasiRegionGrpMats[j]->doImport(*AComp,
                                    *(rowImportPerGrp[j]),
                                    Xpetra::INSERT);

    for (LO row = 0; row < Teuchos::as<LO>(quasiRegionGrpMats[j]->getNodeNumRows()); ++row) { // loop over local rows of composite matrix
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
          commonRegions = findCommonRegions(rowGID, colGID, *regionsPerGIDWithGhosts);
        }

        if (commonRegions.size() > 1) {
          vals[c] *= (1.0 / static_cast<double>(commonRegions.size()));
        }
      }

      quasiRegionGrpMats[j]->replaceLocalValues(row, inds, vals);
    }

    quasiRegionGrpMats[j]->fillComplete();
  }
}

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
  const SC SC_ONE  = Teuchos::ScalarTraits<SC>::one();
  const SC SC_ZERO = Teuchos::ScalarTraits<SC>::zero();

  std::cout << mapComp->getComm()->getRank() << " | Forming region matrices ..." << std::endl;

  // Copy data from quasiRegionGrpMats, but into new map layout
  {
    for (int j = 0; j < maxRegPerProc; j++) {
      regionGrpMats[j] = rcp(new CrsMatrixWrap(revisedRowMapPerGrp[j], revisedColMapPerGrp[j], 9, Xpetra::DynamicProfile));

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

  // enforce nullspace constraint
  std::vector<RCP<Vector> > regNspViolation(maxRegPerProc);
  {
    // compute violation of nullspace property close to DBCs
    RCP<Vector> nspVec = VectorFactory::Build(mapComp);
    nspVec->putScalar(SC_ONE);
    RCP<Vector> nspViolation = VectorFactory::Build(mapComp, true);
    AComp->apply(*nspVec, *nspViolation);

    // move to regional layout
    std::vector<RCP<Vector> > quasiRegNspViolation(maxRegPerProc);
    createRegionalVector(quasiRegNspViolation, maxRegPerProc, rowMapPerGrp);
    createRegionalVector(regNspViolation, maxRegPerProc, revisedRowMapPerGrp);
    compositeToRegional(nspViolation, quasiRegNspViolation, regNspViolation,
                        maxRegPerProc, rowMapPerGrp, revisedRowMapPerGrp, rowImportPerGrp);

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
      std::vector<RCP<Vector> > interfaceScaling(maxRegPerProc);
      for (int j = 0; j < maxRegPerProc; j++) {
        interfaceScaling[j] = VectorFactory::Build(revisedRowMapPerGrp[j]);
        interfaceScaling[j]->putScalar(SC_ONE);
      }

      // transform to composite layout while adding interface values via the Export() combine mode
      RCP<Vector> compInterfaceScalingSum = VectorFactory::Build(mapComp, true);
      regionalToComposite(interfaceScaling, compInterfaceScalingSum,
                          maxRegPerProc, rowMapPerGrp, rowImportPerGrp, Xpetra::ADD);

      /* transform composite layout back to regional layout. Now, GIDs associated
       * with region interface should carry a scaling factor (!= 1).
       */
      std::vector<RCP<Vector> > quasiRegInterfaceScaling(maxRegPerProc);
      compositeToRegional(compInterfaceScalingSum, quasiRegInterfaceScaling,
                          interfaceScaling, maxRegPerProc, rowMapPerGrp,
                          revisedRowMapPerGrp, rowImportPerGrp);

      // modify its interface entries
      for (int j = 0; j < maxRegPerProc; j++) {
        RCP<Vector> inverseInterfaceScaling = VectorFactory::Build(interfaceScaling[j]->getMap(), true);
        inverseInterfaceScaling->reciprocal(*interfaceScaling[j]);
        regNspViolation[j]->elementWiseMultiply(SC_ONE, *regNspViolation[j], *inverseInterfaceScaling, SC_ZERO);
      }
    }
  }

  std::vector<RCP<Vector> > regNsp(maxRegPerProc);
  std::vector<RCP<Vector> > regCorrection(maxRegPerProc);
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
}

#endif // MUELU_SETUPREGIONMATRIX_DEF_HPP
