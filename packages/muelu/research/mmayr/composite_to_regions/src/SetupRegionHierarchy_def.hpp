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
#ifndef MUELU_SETUPREGIONHIERARCHY_DEF_HPP
#define MUELU_SETUPREGIONHIERARCHY_DEF_HPP

#include <vector>
#include <iostream>
#include <numeric>

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

#include <MueLu_CreateXpetraPreconditioner.hpp>

#include "SetupRegionMatrix_def.hpp"

using Teuchos::RCP;
using Teuchos::ArrayRCP;
using Teuchos::Array;
using Teuchos::ArrayView;
using Teuchos::ParameterList;

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
                         const int maxRegPerProc, ///< max number of regions per proc
                         const std::vector<RCP<Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node> > > rowMapPerGrp, ///< row maps in quasiRegion layout [in]
                         const std::vector<RCP<Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node> > > colMapPerGrp, ///< col maps in quasiRegion layout [in]
                         const std::vector<RCP<Xpetra::Import<LocalOrdinal, GlobalOrdinal, Node> > > rowImportPerGrp, ///< row importer in region layout [in]
                         const Xpetra::CombineMode combineMode, ///< Combine mode for import/export [in]
                         RCP<Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> >& compMat ///< Matrix in composite layout [in/out]
                         )
{
#include "Xpetra_UseShortNames.hpp"
  using Teuchos::rcp;
  using std::size_t;

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
      quasiRegMat[j] = rcp(new CrsMatrixWrap(rowMapPerGrp[j], colMapPerGrp[j], 9, Xpetra::DynamicProfile));

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
    partialCompMat[j] = MatrixFactory::Build(compMat->getRowMap(), 3, Xpetra::DynamicProfile);
    partialCompMat[j]->doExport(*(quasiRegMat[j]), *(rowImportPerGrp[j]), Xpetra::INSERT);
    partialCompMat[j]->fillComplete(compMat->getDomainMap(), compMat->getRangeMap());
  }

  // Add all partialCompMat together
  for (int j = 0; j < maxRegPerProc; j++) {
    MatrixMatrix::TwoMatrixAdd(*partialCompMat[j], false, SC_ONE, *compMat, SC_ONE);
  }

  compMat->fillComplete();

  return;
}

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
void createContinuousCoarseLevelMaps(RCP<const Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node> > rowMap, ///< row map
                                     RCP<const Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node> > colMap, ///< column map
                                     RCP<Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node> >& contRowMap, ///< row map with continuous GIDs
                                     RCP<Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node> >& contColMap ///< column map with continuous GIDs
                                     )
{
#include "Xpetra_UseShortNamesOrdinal.hpp"
  //   /!\ This function is pure ordinal, no scalar type is passed as input
  //       This means that use only three template paramters and that we
  //       do not use the Scalar dependent short names!

  // Create row map with continuous GIDs
  contRowMap = MapFactory::Build(rowMap->lib(),
                                 Teuchos::OrdinalTraits<GO>::invalid(),
                                 rowMap->getNodeNumElements(),
                                 Teuchos::OrdinalTraits<GO>::zero(),
                                 rowMap->getComm());

  /* Create column map based on row map with continuous GIDs
   *
   * We use an Importer to create an auxiliary vector in colMap format containing
   * the GIDs of the contRowMap as its entries. By looping over its LIDs, we can
   * then form the contColMap.
   */
  RCP<Import> rowColImport = ImportFactory::Build(rowMap, colMap);
  RCP<Xpetra::Vector<GO, LO, GO, Node> > colGIDVec = Xpetra::VectorFactory<GO, LO, GO, Node>::Build(rowMap, true);
  ArrayRCP<GO> colGIDVecData = colGIDVec->getDataNonConst(0);
  RCP<Xpetra::Vector<GO, LO, GO, Node> > contColGIDVec = Xpetra::VectorFactory<GO, LO, GO, Node>::Build(colMap, true);
  contColGIDVec->doImport(*colGIDVec, *rowColImport, Xpetra::INSERT);

  ArrayRCP<const GO> constColGIDVecData = colGIDVec->getData(0);
  std::vector<GO> contColGIDs;
  for (size_t i = 0; i < contColGIDVec->getLocalLength(); ++i) {
    contColGIDs.push_back(constColGIDVecData[i]);
  }
  contColMap = MapFactory::Build(rowMap->lib(),
                                 Teuchos::OrdinalTraits<GO>::invalid(),
                                 contColGIDs,
                                 Teuchos::OrdinalTraits<GO>::zero(),
                                 rowMap->getComm());

  return;
}



/* Reconstruct coarse-level maps
 *
 * We know the regional map on the coarse levels since they are just the
 * row maps of the coarse level operators. Though, we need to re-construct
 * the quasiRegional and composite maps ourselves.
 */
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void MakeCoarseLevelMaps(const int maxRegPerProc, const int numLevels,
                         Array<RCP<Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node> > >& compRowMaps,
                         Array<std::vector<RCP<Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node> > > >& regRowMaps,
                         Array<std::vector<RCP<Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node> > > >& quasiRegRowMaps,
                         Array<std::vector<RCP<Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node> > > >& regColMaps,
                         Array<std::vector<RCP<Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node> > > >& quasiRegColMaps,
                         Array<std::vector<RCP<Xpetra::Import<LocalOrdinal, GlobalOrdinal, Node> > > >& regRowImporters,
                         Array<std::vector<RCP<Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> > > >& regProlong) {
  /* General strategy to create coarse level maps:
   * =============================================
   *
   * We need three row maps on the next coarser level:
   * - regional row map: just extract the RowMap() from the prolongator
   * - quasiRegional row map: needs to be formed manually
   * - composite row map: needs to be formed manually
   *
   * After extracting the RowMap() from the prolongator to be used as the
   * coarse level's regional map, we need to identify pairs/sets of GIDs
   * that represent the same (duplicated) DOF at a region interface in
   * order to form the quasiRegional and composite map for the next
   * coarser level.
   *
   * Identifying this pairs consists of several steps:
   * 1. Identify GIDs that represent the same (duplicated) DOF at a region
   *   interface on the current level
   * 2. Identify GIDs that represent the same (duplicated) DOF at a region
   *   interface on the next coarser level
   *
   * To form the quasiRegional map, we loop over the regional map and
   * - copy each GID that's associated with an interior DOF away from a
   *   region interface.
   * - replace each GID that's associated with an interface DOF by its
   *   quasiRegional counterpart.
   *
   * Note: We choose the quasiRegional (and composite) GID of an interface
   * DOF to be the std::min of all GIDs that represent that same
   * duplicated DOF.
   *
   * To form the composite map, we loop over the quasiRegional map and
   * - every proc accepts every GID that it owns.
   * - every proc ignores every GID that it doesn't own.
   *
   * To form a composite operator on the coarsest level, we also need the
   * column map on the coarsest level. Hence, we have to recursively
   * re-construct all column maps on each level, namely
   * - regional col map: just extract the ColMap() from the prolongator
   * - quasiRegional col map: needs to be formed manually
   * - composite col map: needs to be formed manually
   *
   * To form column maps, lots of information that we used to generate row
   * maps can be reused. The central piece of information is the mapping
   * of duplicated interface GIDs.
   */

#include "Xpetra_UseShortNames.hpp"

  RCP<Teuchos::FancyOStream> fos = Teuchos::fancyOStream(Teuchos::rcpFromRef(std::cout));
  const SC SC_ONE = Teuchos::ScalarTraits<SC>::one();
  const Xpetra::UnderlyingLib lib = regRowMaps[0][0]->lib();
  const Teuchos::RCP<const Teuchos::Comm<int> > Comm = regRowMaps[0][0]->getComm();

  ////////////////////////////////////////////////////////////////////////
  // IDENTIFY DUPLICATED GIDS ON FINE LEVEL
  ////////////////////////////////////////////////////////////////////////

  RCP<Matrix> summedDuplicateMatrix = Teuchos::null;

  Array<std::vector<std::vector<LocalOrdinal> > > interfaceLIDs; // local IDs of interface nodes on each level in each group
  Array<std::vector<std::vector<std::vector<GlobalOrdinal> > > > interfaceGIDPairs; // pairs of GIDs of interface nodes on each level in each group

  // Find lists of duplicated GIDs locally
  interfaceLIDs.resize(numLevels);
  interfaceGIDPairs.resize(numLevels);
  for (int l = 0; l < numLevels; ++l) {
    interfaceLIDs[l].resize(maxRegPerProc);
    interfaceGIDPairs[l].resize(maxRegPerProc);
  }

  for (int l = 0; l < numLevels - 1; ++l) {
    //        Comm->barrier();
    //        if (Comm->getRank() == 0) {
    //          std::cout << std::endl << std::endl
    //              << "Processing GID pairs on level " << l
    //              << std::endl << std::endl;
    //        }
    //        Comm->barrier();

    //        sleep(1);
    //        std::cout << "Prolongator" << std::endl;
    //        regProlong[l+1][0]->describe(*fos, Teuchos::VERB_HIGH);

    // create list of LIDs per group
    for (int j = 0; j < maxRegPerProc; j++) {
      for (size_t i = 0; i < regRowMaps[l][j]->getNodeNumElements(); ++i) {
        if (regRowMaps[l][j]->getGlobalElement(i) != quasiRegRowMaps[l][j]->getGlobalElement(i)) {
          // This is an interface node
          interfaceLIDs[l][j].push_back(i);
        }
      }
    }

    //        for (int j = 0; j < maxRegPerProc; j++) {
    //          std::cout << myRank << " | group = " << j << " | LIDs = ";
    //          for (std::size_t i = 0; i < interfaceLIDs[l][j].size(); ++i)
    //            std::cout << interfaceLIDs[l][j][i] << ", ";
    //          std::cout << std::endl;
    //        }

    for (int j = 0; j < maxRegPerProc; j++)
      interfaceGIDPairs[l][j].resize(interfaceLIDs[l][j].size());

    // create list of LIDPairs per group
    for (int j = 0; j < maxRegPerProc; j++)
      for (std::size_t i = 0; i < interfaceLIDs[l][j].size(); ++i)
        interfaceGIDPairs[l][j][i].push_back(quasiRegRowMaps[l][j]->getGlobalElement(interfaceLIDs[l][j][i]));
    for (int j = 0; j < maxRegPerProc; j++)
      for (std::size_t i = 0; i < interfaceLIDs[l][j].size(); ++i)
        if (regRowMaps[l][j]->getGlobalElement(interfaceLIDs[l][j][i]) != quasiRegRowMaps[l][j]->getGlobalElement(interfaceLIDs[l][j][i]))
          interfaceGIDPairs[l][j][i].push_back(regRowMaps[l][j]->getGlobalElement(interfaceLIDs[l][j][i]));

    //        std::cout << myRank << " | Print interfaceGIDPairs:" << std::endl;
    //        for (int j = 0; j < maxRegPerProc; j++) {
    //          std::cout << myRank << " | " << "Level " << l <<" | Group " << j << ":" << std::endl;
    //          for (std::size_t i = 0; i < interfaceGIDPairs[l][j].size(); ++i) {
    //            std::cout << "   " << myRank << " | ";
    //            for (std::size_t k = 0; k < interfaceGIDPairs[l][j][i].size(); ++k)
    //              std::cout << interfaceGIDPairs[l][j][i][k] << ", ";
    //            std::cout << std::endl;
    //          }
    //        }

    std::vector<RCP<Vector> > regDupGIDVec(maxRegPerProc); // Vector with zeros everywhere, but ones at duplicated GID spots
    std::vector<RCP<Vector> > quasiRegDupGIDVec(maxRegPerProc); // Vector with zeros everywhere, but ones at duplicated GID spots
    for (int j = 0; j < maxRegPerProc; ++j) {
      regDupGIDVec[j] = VectorFactory::Build(regRowMaps[l][j], true);
      quasiRegDupGIDVec[j] = VectorFactory::Build(quasiRegRowMaps[l][j], true);

      for (std::size_t i = 0; i < interfaceLIDs[l][j].size(); ++i)
        regDupGIDVec[j]->replaceLocalValue(interfaceLIDs[l][j][i], SC_ONE);
    }

    RCP<Vector> compDupGIDVec = VectorFactory::Build(compRowMaps[l], true);
    regionalToComposite(regDupGIDVec, compDupGIDVec, maxRegPerProc,
                        quasiRegRowMaps[l], regRowImporters[l], Xpetra::ADD);

    compositeToRegional(compDupGIDVec, quasiRegDupGIDVec, regDupGIDVec,
                        maxRegPerProc, quasiRegRowMaps[l], regRowMaps[l], regRowImporters[l]);

    // create row/range/domain map for fine level duplicates mapping matrix
    RCP<Map> duplicateMap = Teuchos::null;
    {
      Teuchos::Array<GlobalOrdinal> myIntGIDs;
      for (int j = 0; j < maxRegPerProc; j++) {
        for (std::size_t i = 0; i < interfaceGIDPairs[l][j].size(); ++i) {
          for (std::size_t k = 0; k < interfaceGIDPairs[l][j][i].size(); ++k) {
            if (regRowMaps[l][j]->isNodeGlobalElement(interfaceGIDPairs[l][j][i][k]))
              myIntGIDs.push_back(interfaceGIDPairs[l][j][i][k]);
          }
        }
      }

      duplicateMap = MapFactory::Build(lib, Teuchos::OrdinalTraits<GlobalOrdinal>::invalid(),
                                       myIntGIDs,Teuchos::OrdinalTraits<GlobalOrdinal>::zero(), Comm);

      //          sleep(1);
      //          Comm->barrier();
      //          std::cout << myRank << " | duplicateMap:" << std::endl;
      //          duplicateMap->describe(*fos, Teuchos::VERB_HIGH);
    }

    // create row/range/domain map for the transpose of the fine level duplicates mapping matrix
    RCP<Map> fullDuplicateMap = Teuchos::null;
    {
      Array<GlobalOrdinal> myIntGIDs;
      ArrayRCP<const Scalar> regDupGIDVecData = regDupGIDVec[0]->getData(0);
      for (size_t i = 0; i < regRowMaps[l][0]->getNodeNumElements(); ++i) {
        if (regDupGIDVecData[i] != 0)
          myIntGIDs.push_back(Teuchos::as<GlobalOrdinal>(regRowMaps[l][0]->getGlobalElement(i)));
      }

      //          std::cout << myRank << " | myIntGIDs = ";
      //          for (auto gid : myIntGIDs)
      //            std::cout << gid << ", ";
      //          std::cout << std::endl;

      fullDuplicateMap = MapFactory::Build(lib, Teuchos::OrdinalTraits<GlobalOrdinal>::invalid(),
                                           myIntGIDs, Teuchos::OrdinalTraits<GlobalOrdinal>::zero(), Comm);

      //          sleep(1);
      //          Comm->barrier();
      //          std::cout << myRank << " | fullDuplicateMap:" << std::endl;
      //          fullDuplicateMap->describe(*fos, TEUCHOS::VERB_HIGH);
    }

    // Create and fill matrix
    RCP<Matrix> duplicateMatrix = MatrixFactory::Build(fullDuplicateMap, 2, Xpetra::DynamicProfile);
    {
      // Fill diagonal
      {
        Array<Scalar> vals(1);
        Array<GlobalOrdinal> cols(1);
        for (size_t i = 0; i < fullDuplicateMap->getNodeNumElements(); ++i) {
          GlobalOrdinal globalRow = fullDuplicateMap->getGlobalElement(i);

          vals[0] = static_cast<Scalar>(globalRow);
          cols[0] = globalRow;

          duplicateMatrix->insertGlobalValues(globalRow, cols, vals);
        }
      }

      // Fill off-diagonals (if known)
      Array<Scalar> vals(2);
      Array<GlobalOrdinal> cols(2);
      for (size_t i = 0; i < duplicateMap->getNodeNumElements(); ++i) {
        GlobalOrdinal globalRow = duplicateMap->getGlobalElement(i);

        for (std::size_t k = 0; k < interfaceGIDPairs[l][0][i].size(); ++k) {
          vals[k] = static_cast<Scalar>(globalRow);
          cols[k] = interfaceGIDPairs[l][0][i][k];
        }

        duplicateMatrix->insertGlobalValues(globalRow, cols, vals);
      }
    }

    duplicateMatrix->fillComplete();

    //        sleep(1);
    //        Comm->barrier();
    //        std::cout << myRank << " | duplicateMatrix:" << std::endl;
    //        duplicateMatrix->describe(*fos, Teuchos::VERB_HIGH);
    //
    //        sleep(1);
    //        Comm->barrier();
    //        std::cout << myRank << " | duplicateMatrix->RangeMap():" << std::endl;
    //        duplicateMatrix->getRangeMap()->describe(*fos, Teuchos::VERB_HIGH);
    //
    //        sleep(1);
    //        Comm->barrier();
    //        std::cout << myRank << " | duplicateMatrix->DomainMap():" << std::endl;
    //        duplicateMatrix->getDomainMap()->describe(*fos, Teuchos::VERB_HIGH);
    //
    //        sleep(1);
    //        Comm->barrier();
    //        std::cout << myRank << " | duplicateMatrix->ColMap():" << std::endl;
    //        duplicateMatrix->getColMap()->describe(*fos, Teuchos::VERB_HIGH);

    {
      RCP<Matrix> transDuplicateMatrix = MueLu::Utilities<Scalar,LocalOrdinal,GlobalOrdinal,Node>::Transpose(*duplicateMatrix);
      TEUCHOS_ASSERT(!transDuplicateMatrix.is_null());

      transDuplicateMatrix->fillComplete();
      TEUCHOS_ASSERT(transDuplicateMatrix->isFillComplete());

      //          sleep(1);
      //          Comm->barrier();
      //          std::cout << myRank << " | transDuplicateMatrix:" << std::endl;
      //          transDuplicateMatrix->describe(*fos, Teuchos::VERB_HIGH);

      summedDuplicateMatrix = MatrixFactory::Build(fullDuplicateMap, 2, Xpetra::DynamicProfile);
      MatrixMatrix::TwoMatrixAdd(*transDuplicateMatrix, false, SC_ONE, *duplicateMatrix, false, SC_ONE, summedDuplicateMatrix, *fos);
    }

    summedDuplicateMatrix->fillComplete(fullDuplicateMap, fullDuplicateMap);
    //        summedDuplicateMatrix->fillComplete();

    //        sleep(1);
    //        Comm->barrier();
    //        std::cout << myRank << " | summedDuplicateMatrix:" << std::endl;
    //        summedDuplicateMatrix->describe(*fos, Teuchos::VERB_HIGH);

    std::vector<std::vector<GlobalOrdinal> > myDuplicates; // pairs of duplicates with locally owned GID listed first.
    myDuplicates.resize(summedDuplicateMatrix->getNodeNumRows());
    for (size_t lRow = 0; lRow < summedDuplicateMatrix->getNodeNumRows(); ++lRow) {

      ArrayView<const LocalOrdinal> inds;
      ArrayView<const Scalar> vals;

      summedDuplicateMatrix->getLocalRowView(lRow, inds, vals);

      std::vector<GlobalOrdinal> gidsToSort;
      for (typename ArrayView<const LocalOrdinal>::size_type k = 0; k < inds.size(); ++k) {
        //            std::cout << myRank << " | inds[" << k << "] = " << inds[k] << std::endl;
        gidsToSort.push_back(summedDuplicateMatrix->getColMap()->getGlobalElement(inds[k]));
      }

      /* Note: At intersections of more than two regions, i.e. at interior
       * vertices of the regional interface, the list of gidsToSort is not
       * identical on all procs. In particular, only the proc owning the
       * original GID of this vertex knows the full list of all duplicated
       * GIDs. Those procs, that own a duplicated GID, only know about the
       * two-member pair of their duplicate and the original GID, but not
       * about the duplicated GIDs of all other procs.
       *
       * However, this does not matter since we are only looking for the
       * mapping of each proc's duplicated GID to the original GID, which
       * every proc knows about.
       *
       * Bottom line, this should be sufficient.
       */

      //          sleep(1);
      //          std::cout << myRank << " | gidsToSort:" << std::endl << "  ";
      //          for (auto gid : gidsToSort)
      //             std::cout << gid << ", ";
      //           std::cout << std::endl;
      //           sleep(1);

      /* sort s.t. my GID is first
       *
       * 1. Find index of the one GID that I own
       * 2. Insert this GID at the beginning of the vector
       * 3. Erase this GID from its initial position in the vector
       *
       * ToDo (mayr.mt) Is this really necessary?
       */
      {
        //            int indOfMyGID = -1;
        //            for (std::size_t k = 0; k < gidsToSort.size(); ++k) {
        //              if (regRowMaps[l][0]->MyGID(gidsToSort[k])) {
        //                indOfMyGID = k;
        //                break;
        //              }
        //            }
        //            TEUCHOS_ASSERT(indOfMyGID >= 0);
        //
        //            int tmpIndex = gidsToSort[indOfMyGID];
        //            gidsToSort.erase(gidsToSort.begin() + indOfMyGID);
        //            gidsToSort.insert(gidsToSort.begin(), tmpIndex);

        //            for (std::size_t k = 0; i < gidsToSort.size(); ++k)
        //              myDuplicates[i].push_back(gidsToSort[i]);

        myDuplicates[lRow] = gidsToSort;
      }
    }

    //        sleep(myRank);
    //        std::cout << std::endl << std::endl << std::endl << std::endl;
    //        for (std::size_t i = 0; i < myDuplicates.size(); ++i) {
    //          std::cout << myRank << " | myDuplicates:" << std::endl << "  ";
    //          for (std::size_t k = 0; k < myDuplicates[i].size(); ++k)
    //            std::cout << myDuplicates[i][k] << ", ";
    //          std::cout << std::endl;
    //        }
    //        std::cout << std::endl << std::endl << std::endl << std::endl;

    ////////////////////////////////////////////////////////////////////////
    // IDENTIFY DUPLICATED GIDS ON COARSE LEVEL
    ////////////////////////////////////////////////////////////////////////

    std::vector<std::vector<GlobalOrdinal> > myCoarseDuplicates; // pairs of duplicates with locally owned GID listed first on coarse level.
    myCoarseDuplicates.resize(myDuplicates.size());

    std::vector<GlobalOrdinal> myCoarseInterfaceGIDs;
    for (size_t i = 0; i < myDuplicates.size(); ++i) {
      const GlobalOrdinal rowGID = myDuplicates[i][0];
      const LocalOrdinal rowLID = regProlong[l+1][0]->getRowMap()->getLocalElement(rowGID);

      ArrayView<const Scalar> vals;
      ArrayView<const LocalOrdinal> inds;
      regProlong[l+1][0]->getLocalRowView(rowLID, inds, vals);
      //          std::cout << myRank << " | ExtractMyRowView err = " << err << std::endl;
      TEUCHOS_ASSERT(inds.size() == 1); // tentative P: only one entry per row
      TEUCHOS_ASSERT(vals.size() == 1); // tentative P: only one entry per row

      myCoarseInterfaceGIDs.push_back(regProlong[l+1][0]->getColMap()->getGlobalElement(inds[0]));
    }

    //        std::cout << "myCoarseInterfaceGIDs on proc " << myRank << ": ";
    //        for (auto gid : myCoarseInterfaceGIDs) {
    //          std::cout << gid << ", ";
    //        }
    //        std::cout << std::endl;

    // Build row/range/domain map of duplicate mapping matrix
    RCP<Map> interfaceMap = Teuchos::null;
    {
      std::vector<GlobalOrdinal> interfaceMapGIDs(myDuplicates.size());
      for (std::size_t i = 0; i < myDuplicates.size(); ++i) {
        interfaceMapGIDs[i] = myDuplicates[i][0];
      }

      //          std::cout << "interfaceMapGIDs on proc " << myRank << ":      ";
      //          for (auto gid : interfaceMapGIDs) {
      //            std::cout << gid << ", ";
      //          }
      //          std::cout << std::endl;

      interfaceMap = MapFactory::Build(lib, Teuchos::OrdinalTraits<GlobalOrdinal>::invalid(),
                                       interfaceMapGIDs, Teuchos::OrdinalTraits<GlobalOrdinal>::zero(), Comm);
    }

    RCP<Matrix> duplicateMapping = MatrixFactory::Build(interfaceMap, 2, Xpetra::DynamicProfile);

    // Fill the matrix
    RCP<Matrix> transDuplicateMapping = Teuchos::null;
    {
      size_t numRows = myDuplicates.size();
      std::vector<GlobalOrdinal> rowPtr(numRows);
      std::vector<Array<Scalar> > vals(numRows);
      std::vector<Array<LocalOrdinal> > colInds(numRows);

      for (size_t rowIdx = 0; rowIdx < numRows; ++rowIdx) {
        rowPtr[rowIdx] = myDuplicates[rowIdx][0];
        for (std::size_t k = 0; k < myDuplicates[rowIdx].size(); ++k) {
          vals[rowIdx].push_back(myCoarseInterfaceGIDs[rowIdx]);
          colInds[rowIdx].push_back(myDuplicates[rowIdx][k]);
        }
      }

      //          for (int rowIdx = 0; rowIdx < numRows; ++rowIdx) {
      //            std::cout << myRank << " | rowIdx = " << rowIdx << ", rowGID = " << rowPtr[rowIdx] << ":";
      //            for (std::size_t i = 0; i < vals[rowIdx].size(); ++i)
      //              std::cout << "(" << colInds[rowIdx][i] << "|" << vals[rowIdx][i] << "), ";
      //            std::cout << std::endl;
      //          }

      // local dummy insertion
      for (size_t rowIdx = 0; rowIdx < numRows; ++rowIdx) {
        duplicateMapping->insertGlobalValues(rowPtr[rowIdx], colInds[rowIdx], vals[rowIdx]);
      }

      duplicateMapping->fillComplete();
      TEUCHOS_ASSERT(duplicateMapping->isFillComplete());

      //          sleep(1);
      //          Comm.Barrier();
      //          duplicateMapping->Print(std::cout);

      transDuplicateMapping = MueLu::Utilities<Scalar,LocalOrdinal,GlobalOrdinal,Node>::Transpose(*duplicateMapping);

      Comm->barrier();

      TEUCHOS_ASSERT(!transDuplicateMapping.is_null());

      transDuplicateMapping->fillComplete();
      TEUCHOS_ASSERT(transDuplicateMapping->isFillComplete());

      //          sleep(1);
      //          std::cout << myRank << " | Printing the transpose ..." << std::endl;
      //          Comm.Barrier();
      //          transDuplicateMapping->Print(std::cout);
    }

    sleep(1);

    /* Extract coarse level duplicates from transDuplicateMapping
     *
     * Note: In 2D, there will be duplicates which need to be removed.
     * One could maybe use something along the lines of:
     *    sort( vec.begin(), vec.end() );
     *    vec.erase( unique( vec.begin(), vec.end() ), vec.end() );
     *
     */
    std::vector<std::vector<GlobalOrdinal> > myCoarseInterfaceDuplicates(transDuplicateMapping->getNodeNumRows());
    {
      for (size_t lRow = 0; lRow < transDuplicateMapping->getNodeNumRows(); ++lRow) {
        ArrayView<const Scalar> vals;
        ArrayView<const LocalOrdinal> inds;
        transDuplicateMapping->getLocalRowView(lRow, inds, vals);

        const size_t numEntries = inds.size();
        myCoarseInterfaceDuplicates[lRow].resize(numEntries);
        myCoarseInterfaceDuplicates[lRow].assign(vals.begin(), vals.end());

        //            std::cout << myRank << " | myCoarseInterfaceDuplicates[" << i << "] = ";
        //            for (auto id : myCoarseInterfaceDuplicates[i])
        //              std::cout << id << ", ";
        //            std::cout << std::endl;
      }

    }

    ////////////////////////////////////////////////////////////////////////
    // CREATE COARSE LEVEL MAPS
    ////////////////////////////////////////////////////////////////////////
    {
      //          sleep(1);
      //          std::cout << myRank << " | Printing regRowMaps[" << l+1 << "][0] ..." << std::endl;
      //          Comm->barrier();
      //          regRowMaps[l+1][0]->describe(*fos, Teuchos::VERB_HIGH);
      //          sleep(2);

      // create quasiRegional row map
      {
        std::vector<GlobalOrdinal> myQuasiRegGIDs;

        for (size_t i = 0; i < regRowMaps[l+1][0]->getNodeNumElements(); ++i) {
          // grab current regional GID to be processed
          GlobalOrdinal currGID = regRowMaps[l+1][0]->getGlobalElement(i);
          GlobalOrdinal quasiGID = currGID; // assign dummy value

          // find quasiRegional counterpart
          {
            // Is this an interface GID?
            bool isInterface = false;
            for (std::size_t k = 0; k < myCoarseInterfaceDuplicates.size(); ++k) {
              for (std::size_t kk = 0; kk < myCoarseInterfaceDuplicates[k].size(); ++kk) {
                if (currGID == myCoarseInterfaceDuplicates[k][kk])
                  isInterface = true;
              }
            }

            if (isInterface) {
              for (std::size_t k = 0; k < myCoarseInterfaceDuplicates.size(); ++k) {
                bool found = false;
                for (std::size_t kk = 0; kk < myCoarseInterfaceDuplicates[k].size(); ++kk) {
                  if (currGID == myCoarseInterfaceDuplicates[k][kk])
                    found = true;
                }

                if (found) {
                  for (std::size_t kk = 0; kk < myCoarseInterfaceDuplicates[k].size(); ++kk)
                    quasiGID = std::min(quasiGID, Teuchos::as<const GlobalOrdinal>(myCoarseInterfaceDuplicates[k][kk]));

                  //                      std::cout << myRank << " | Interface GID " << currGID << " is replaced by quasiGID " << quasiGID << std::endl;
                  break;
                }
              }
            }
            else { // interior node --> take GID from regional map
              quasiGID = currGID;
            }
          }

          TEUCHOS_ASSERT(quasiGID>=0);
          myQuasiRegGIDs.push_back(quasiGID);
        }

        quasiRegRowMaps[l+1][0] = MapFactory::Build(lib, Teuchos::OrdinalTraits<GlobalOrdinal>::invalid(),
                                                    myQuasiRegGIDs, Teuchos::OrdinalTraits<GlobalOrdinal>::zero(), Comm);
        TEUCHOS_ASSERT(!quasiRegRowMaps[l+1][0].is_null());

        //            sleep(1);
        //            std::cout << myRank << " | Printing quasiRegRowMaps[" << l+1 << "][0] ..." << std::endl;
        //            Comm->barrier();
        //            quasiRegRowMaps[l+1][0]->describe(*fos, Teuchos::VERB_HIGH);
      }

      // create composite row map
      {
        std::vector<GlobalOrdinal> myCompGIDs;
        for (int j = 0; j < maxRegPerProc; j++) {
          for (size_t i = 0; i < quasiRegRowMaps[l+1][j]->getNodeNumElements(); ++i) {
            const GlobalOrdinal trialGID = quasiRegRowMaps[l+1][j]->getGlobalElement(i);

            if (regRowMaps[l+1][j]->isNodeGlobalElement(trialGID))
              myCompGIDs.push_back(trialGID);
          }
        }

        compRowMaps[l+1] = MapFactory::Build(lib, Teuchos::OrdinalTraits<GlobalOrdinal>::invalid(),
                                             myCompGIDs, Teuchos::OrdinalTraits<GlobalOrdinal>::zero(), Comm);
        TEUCHOS_ASSERT(!compRowMaps[l+1].is_null());

        //            sleep(1);
        //            std::cout << myRank << " | Printing compRowMaps["<< l+1 << "] ..." << std::endl;
        //            Comm->barrier();
        //            compRowMaps[l+1]->describe(*fos, Teuchos::VERB_HIGH);
      }

      // create regRowImporter
      for (int j = 0; j < maxRegPerProc; ++j) {
        regRowImporters[l+1][j] = ImportFactory::Build(compRowMaps[l+1], quasiRegRowMaps[l+1][j]);
        TEUCHOS_ASSERT(!regRowImporters[l+1][j].is_null());
      }

      // Create quasiRegional column map
      {
        std::vector<GlobalOrdinal> myQuasiRegGIDs;

        for (size_t i = 0; i < regColMaps[l+1][0]->getNodeNumElements(); ++i) {
          // grab current regional GID to be processed
          GlobalOrdinal currGID = regColMaps[l+1][0]->getGlobalElement(i);
          GlobalOrdinal quasiGID = currGID; // assign dummy value

          /* Find quasiRegional counterpart
           *
           * We should be able to use the same procedure as for the
           * quasiRegRowMap since duplicated GIDs only live on the interface.
           * Hence, even if the column map reaches across an interface,
           * GIDs from 'the other side of the interface' do not have to be
           * modified as they have not been duplicated (only those on the
           * interface).
           */
          {
            // Is this an interface GID?
            bool isInterface = false;
            for (std::size_t k = 0; k < myCoarseInterfaceDuplicates.size(); ++k) {
              for (std::size_t kk = 0; kk < myCoarseInterfaceDuplicates[k].size(); ++kk) {
                if (currGID == myCoarseInterfaceDuplicates[k][kk])
                  isInterface = true;
              }
            }

            if (isInterface) {
              for (std::size_t k = 0; k < myCoarseInterfaceDuplicates.size(); ++k) {
                bool found = false;
                for (std::size_t kk = 0; kk < myCoarseInterfaceDuplicates[k].size(); ++kk) {
                  if (currGID == myCoarseInterfaceDuplicates[k][kk])
                    found = true;
                }

                if (found) {
                  for (std::size_t kk = 0; kk < myCoarseInterfaceDuplicates[k].size(); ++kk)
                    quasiGID = std::min(quasiGID, Teuchos::as<GlobalOrdinal>(myCoarseInterfaceDuplicates[k][kk]));

                  //                      std::cout << myRank << " | Interface GID " << currGID << " is replaced by quasiGID " << quasiGID << std::endl;
                  break;
                }
              }
            }
            else { // interior node --> take GID from regional map
              quasiGID = currGID;
            }
          }

          TEUCHOS_ASSERT(quasiGID>=0);
          myQuasiRegGIDs.push_back(quasiGID);
        }

        quasiRegColMaps[l+1][0] = MapFactory::Build(lib, Teuchos::OrdinalTraits<GlobalOrdinal>::invalid(),
                                                    myQuasiRegGIDs, Teuchos::OrdinalTraits<GlobalOrdinal>::zero(), Comm);
        TEUCHOS_ASSERT(!quasiRegColMaps[l+1][0].is_null());

        //            sleep(1);
        //            std::cout << myRank << " | Printing quasiRegColMaps[" << l+1 << "][0] ..." << std::endl;
        //            Comm->barrier();
        //            quasiRegColMaps[l+1][0]->describe(*fos, Teuchos::VERB_HIGH);
      }
    }
  }
}


// Form the composite coarse level operator
template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void MakeCoarseCompositeOperator(const int maxRegPerProc, const int numLevels,
                                 Array<RCP<Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node> > >& compRowMaps,
                                 Array<std::vector<RCP<Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node> > > >& quasiRegRowMaps,
                                 Array<std::vector<RCP<Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node> > > >& quasiRegColMaps,
                                 Array<std::vector<RCP<Xpetra::Import<LocalOrdinal, GlobalOrdinal, Node> > > >& regRowImporters,
                                 Array<std::vector<RCP<Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> > > >& regMatrices,
                                 RCP<Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> >& coarseCompOp)
{
#include "Xpetra_UseShortNames.hpp"
  const int maxLevel = numLevels - 1;
  coarseCompOp = MatrixFactory::Build(compRowMaps[maxLevel], 3, Xpetra::DynamicProfile);
  //      coarseCompOp->setAllToScalar(SC_ZERO);
  //      coarseCompOp->describe(*fos, Teuchos::VERB_EXTREME);

  regionalToComposite(regMatrices[maxLevel], maxRegPerProc,
                      quasiRegRowMaps[maxLevel], quasiRegColMaps[maxLevel],
                      regRowImporters[maxLevel], Xpetra::ADD,
                      coarseCompOp);

  //      coarseCompOp->fillComplete(compRowMaps[maxLevel], compRowMaps[maxLevel]);
  //      TEUCHOS_ASSERT(coarseCompOp->isFillComplete());
  //
  //      sleep(1);
  //      std::cout << myRank << " | Printing coarseCompOp ..." << std::endl;
  //      Comm->barrier();
  //      coarseCompOp->describe(*fos, Teuchos::VERB_HIGH);
}


  // Make interface scaling factors recursively
template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void MakeInterfaceScalingFactors(const int numLevels, const int maxRegPerProc,
                                 Array<RCP<Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node> > >& compRowMaps,
                                 Array<std::vector<RCP<Xpetra::Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node> > > >& regInterfaceScalings,
                                 Array<std::vector<RCP<Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node> > > >& regRowMaps,
                                 Array<std::vector<RCP<Xpetra::Import<LocalOrdinal, GlobalOrdinal, Node> > > >& regRowImporters,
                                 Array<std::vector<RCP<Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node> > > >& quasiRegRowMaps)
{
#include "Xpetra_UseShortNames.hpp"
  std::cout << compRowMaps[0]->getComm()->getRank() << " | Computing interface scaling factors ..." << std::endl;

  const SC SC_ONE = Teuchos::ScalarTraits<SC>::one();

  TEUCHOS_TEST_FOR_EXCEPT_MSG(!(numLevels>0), "We require numLevel > 0. Probably, numLevel has not been set, yet.");

  for (int l = 0; l < numLevels; l++) {
    // initialize region vector with all ones.
    for (int j = 0; j < maxRegPerProc; j++) {
      regInterfaceScalings[l][j] = VectorFactory::Build(regRowMaps[l][j]);
      regInterfaceScalings[l][j]->putScalar(SC_ONE);
    }

    // transform to composite layout while adding interface values via the Export() combine mode
    RCP<Vector> compInterfaceScalingSum = VectorFactory::Build(compRowMaps[l], true);
    regionalToComposite(regInterfaceScalings[l], compInterfaceScalingSum, maxRegPerProc, quasiRegRowMaps[l], regRowImporters[l], Xpetra::ADD);

    /* transform composite layout back to regional layout. Now, GIDs associated
     * with region interface should carry a scaling factor (!= 1).
     */
    std::vector<RCP<Vector> > quasiRegInterfaceScaling(maxRegPerProc); // Is that vector really needed?
    compositeToRegional(compInterfaceScalingSum, quasiRegInterfaceScaling,
                        regInterfaceScalings[l], maxRegPerProc, quasiRegRowMaps[l],
                        regRowMaps[l], regRowImporters[l]);
  }
}


template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void createRegionHierarchy(const int maxRegPerProc,
                           const int numDimensions,
                           const Array<Array<int> > lNodesPerDim,
                           const Array<std::string> aggregationRegionType,
                           const std::string xmlFileName,
                           Array<RCP<Xpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node> > >& nullspace,
                           Array<RCP<Xpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node> > >& coordinates,
                           std::vector<RCP<Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> > >& regionGrpMats,
                           const RCP<Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node> > mapComp,
                           const std::vector<RCP<Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node> > > rowMapPerGrp,
                           const std::vector<RCP<Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node> > > colMapPerGrp,
                           const std::vector<RCP<Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node> > > revisedRowMapPerGrp,
                           const std::vector<RCP<Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node> > > revisedColMapPerGrp,
                           const std::vector<RCP<Xpetra::Import<LocalOrdinal, GlobalOrdinal, Node> > > rowImportPerGrp,
                           Array<RCP<Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node> > >& compRowMaps,
                           Array<RCP<Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node> > >& compColMaps,
                           Array<std::vector<RCP<Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node> > > >& regRowMaps,
                           Array<std::vector<RCP<Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node> > > >& regColMaps,
                           Array<std::vector<RCP<Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node> > > >& quasiRegRowMaps,
                           Array<std::vector<RCP<Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node> > > >& quasiRegColMaps,
                           Array<std::vector<RCP<Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> > > >& regMatrices,
                           Array<std::vector<RCP<Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> > > >& regProlong,
                           Array<std::vector<RCP<Xpetra::Import<LocalOrdinal, GlobalOrdinal, Node> > > >& regRowImporters,
                           Array<std::vector<RCP<Xpetra::Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node> > > >& regInterfaceScalings,
                           RCP<Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> >& coarseCompOp) {
#include "Xpetra_UseShortNames.hpp"

  typedef MueLu::Hierarchy<Scalar, LocalOrdinal, GlobalOrdinal, Xpetra::EpetraNode> Hierarchy;
  typedef MueLu::Utilities<Scalar, LocalOrdinal, GlobalOrdinal, Xpetra::EpetraNode> Utilities;

  std::cout << mapComp->getComm()->getRank() << " | Setting up MueLu hierarchies ..." << std::endl;
  int numLevels = 0;

  // A hierarchy for each group
  std::vector<RCP<Hierarchy> > regGrpHierarchy(maxRegPerProc);

  for (int j = 0; j < maxRegPerProc; j++) {

    /* Set number of nodes per processor per dimension
     *
     * We don't use the number of owned nodes provided on input.
     * Use the region dimensions instead. This is the right thing to do
     * since duplication of interface nodes has added duplicated nodes to those regions
     * where inpData_ownedX/inpData_ownedY and inpData_regionX/inpData_regionY have been different on input.
     */

    // Read MueLu parameter list form xml file
    RCP<ParameterList> mueluParams = Teuchos::getParametersFromXmlFile(xmlFileName);

    // Insert region-specific data into parameter list
    const std::string userName = "user data";
    Teuchos::ParameterList& userParamList = mueluParams->sublist(userName);
    userParamList.set<int>("int numDimensions", numDimensions);
    userParamList.set<Array<int> >("Array<LO> lNodesPerDim", lNodesPerDim[j]);
    userParamList.set<std::string>("string aggregationRegionType", aggregationRegionType[j]);

    // Setup hierarchy
    regGrpHierarchy[j] = MueLu::CreateXpetraPreconditioner(regionGrpMats[j], *mueluParams, coordinates[j], nullspace[j]);
  }

  // resize Arrays and vectors
  {
    // resize level containers
    numLevels = regGrpHierarchy[0]->GetNumLevels();
    compRowMaps.resize(numLevels);
    compColMaps.resize(numLevels);
    regRowMaps.resize(numLevels);
    regColMaps.resize(numLevels);
    quasiRegRowMaps.resize(numLevels);
    quasiRegColMaps.resize(numLevels);
    regMatrices.resize(numLevels);
    regProlong.resize(numLevels);
    regRowImporters.resize(numLevels);
    regInterfaceScalings.resize(numLevels);

    // resize group containers on each level
    for (int l = 0; l < numLevels; ++l) {
      regRowMaps[l].resize(maxRegPerProc);
      regColMaps[l].resize(maxRegPerProc);
      quasiRegRowMaps[l].resize(maxRegPerProc);
      quasiRegColMaps[l].resize(maxRegPerProc);
      regMatrices[l].resize(maxRegPerProc);
      regProlong[l].resize(maxRegPerProc);
      regRowImporters[l].resize(maxRegPerProc);
      regInterfaceScalings[l].resize(maxRegPerProc);
    }
  }

  // Fill fine level with our data
  {
    compRowMaps[0] = mapComp;
    quasiRegRowMaps[0] = rowMapPerGrp;
    quasiRegColMaps[0] = colMapPerGrp;
    regRowMaps[0] = revisedRowMapPerGrp;
    regColMaps[0] = revisedColMapPerGrp;
    regRowImporters[0] = rowImportPerGrp;
    regMatrices[0] = regionGrpMats;

    /* MueLu stores prolongator on coarse level, so there is no prolongator
     * on the fine level. To have level containers of the same size, let's
     * just put in dummy data
     */
    std::vector<RCP<Matrix> > fineLevelProlong(maxRegPerProc);
    for (int j = 0; j < maxRegPerProc; ++j)
      fineLevelProlong[j] = Teuchos::null;
    regProlong[0] = fineLevelProlong;
  }

  /* Get coarse level matrices and prolongators from MueLu hierarchy
   * Note: fine level has been dealt with previously, so we start at level 1 here.
   */
  for (int l = 1; l < numLevels; ++l) { // Note: we start at level 1 (which is the first coarse level)
    for (int j = 0; j < maxRegPerProc; ++j) {
      RCP<MueLu::Level> level = regGrpHierarchy[j]->GetLevel(l);

      regProlong[l][j] = level->Get<RCP<Matrix> >("P", MueLu::NoFactory::get());
      regMatrices[l][j] = level->Get<RCP<Matrix> >("A", MueLu::NoFactory::get());

      regRowMaps[l][j] = Teuchos::rcp_const_cast<Map>(regMatrices[l][j]->getRowMap()); // ToDo (mayr.mt) Should we rather copy?
      regColMaps[l][j] = Teuchos::rcp_const_cast<Map>(regMatrices[l][j]->getColMap()); // ToDo (mayr.mt) Should we rather copy?
    }
  }

  MakeCoarseLevelMaps(maxRegPerProc, numLevels,
                      compRowMaps,
                      regRowMaps,
                      quasiRegRowMaps,
                      regColMaps,
                      quasiRegColMaps,
                      regRowImporters,
                      regProlong);

  MakeCoarseCompositeOperator(maxRegPerProc, numLevels,
                              compRowMaps,
                              quasiRegRowMaps,
                              quasiRegColMaps,
                              regRowImporters,
                              regMatrices,
                              coarseCompOp);

  MakeInterfaceScalingFactors(numLevels, maxRegPerProc,
                              compRowMaps,
                              regInterfaceScalings,
                              regRowMaps,
                              regRowImporters,
                              quasiRegRowMaps);
}





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
}


/*! \brief Compute the residual \f$r = b - Ax\f$
 *
 *  The residual is computed based on matrices and vectors in a regional layout.
 *  1. Compute y = A*x in regional layout.
 *  2. Sum interface values of y to account for duplication of interface DOFs.
 *  3. Compute r = b - y
 */
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
std::vector<RCP<Xpetra::Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node> > >
computeResidual(std::vector<RCP<Xpetra::Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node> > >& regRes, ///< residual (to be evaluated)
                const std::vector<RCP<Xpetra::Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node> > > regX, ///< left-hand side (solution)
                const std::vector<RCP<Xpetra::Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node> > > regB, ///< right-hand side (forcing term)
                const std::vector<RCP<Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> > > regionGrpMats,
                const RCP<Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node> > mapComp, ///< composite map, computed by removing GIDs > numDofs in revisedRowMapPerGrp
                const std::vector<RCP<Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node> > > rowMapPerGrp, ///< row maps in region layout [in] requires the mapping of GIDs on fine mesh to "filter GIDs"
                const std::vector<RCP<Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node> > > revisedRowMapPerGrp, ///< revised row maps in region layout [in] (actually extracted from regionGrpMats)
                const std::vector<RCP<Xpetra::Import<LocalOrdinal, GlobalOrdinal, Node> > > rowImportPerGrp ///< row importer in region layout [in]
    )
{
#include "Xpetra_UseShortNames.hpp"
  const int maxRegPerProc = regX.size();

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

  sumInterfaceValues(regRes, mapComp, maxRegPerProc, rowMapPerGrp,
                     revisedRowMapPerGrp, rowImportPerGrp);

  for (int j = 0; j < maxRegPerProc; j++) { // step 3
    regRes[j]->update(1.0, *regB[j], -1.0);
    //    TEUCHOS_ASSERT(regRes[j]->getMap()->isSameAs(*regB[j]->getMap()));
  }

  return regRes;
}


/*! \brief Do Jacobi smoothing
 *
 *  Perform Jacobi smoothing in the region layout using the true diagonal value
 *  recovered from the splitted matrix.
 */
template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void jacobiIterate(const int maxIter,
                   const double omega,
                   std::vector<RCP<Xpetra::Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node> > >& regX, // left-hand side (or solution)
                   const std::vector<RCP<Xpetra::Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node> > > regB, // right-hand side (or residual)
                   const std::vector<RCP<Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> > > regionGrpMats, // matrices in true region layout
                   const std::vector<RCP<Xpetra::Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node> > > regionInterfaceScaling, // recreate on coarse grid by import Add on region vector of ones
                   const int maxRegPerProc, ///< max number of regions per proc [in]
                   const RCP<Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node> > mapComp, ///< composite map
                   const std::vector<RCP<Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node> > > rowMapPerGrp, ///< row maps in region layout [in] requires the mapping of GIDs on fine mesh to "filter GIDs"
                   const std::vector<RCP<Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node> > > revisedRowMapPerGrp, ///< revised row maps in region layout [in] (actually extracted from regionGrpMats)
                   const std::vector<RCP<Xpetra::Import<LocalOrdinal, GlobalOrdinal, Node> > > rowImportPerGrp ///< row importer in region layout [in]
    )
{
#include "Xpetra_UseShortNames.hpp"
  const Scalar SC_ZERO = Teuchos::ScalarTraits<Scalar>::zero();
  const Scalar SC_ONE = Teuchos::ScalarTraits<Scalar>::one();

  std::vector<RCP<Vector> > regRes(maxRegPerProc);
  createRegionalVector(regRes, maxRegPerProc, revisedRowMapPerGrp);

  // extract diagonal from region matrices, recover true diagonal values, invert diagonal
  std::vector<RCP<Vector> > diag(maxRegPerProc);
  for (int j = 0; j < maxRegPerProc; j++) {
    // extract inverse of diagonal from matrix
    diag[j] = VectorFactory::Build(regionGrpMats[j]->getRowMap(), true);
    regionGrpMats[j]->getLocalDiagCopy(*diag[j]);
    diag[j]->elementWiseMultiply(SC_ONE, *diag[j], *regionInterfaceScaling[j], SC_ZERO); // ToDo Does it work to pass in diag[j], but also return into the same variable?
    diag[j]->reciprocal(*diag[j]);
  }

  for (int iter = 0; iter < maxIter; ++iter) {

    /* Update the residual vector
     * 1. Compute tmp = A * regX in each region
     * 2. Sum interface values in tmp due to duplication (We fake this by scaling to reverse the basic splitting)
     * 3. Compute r = B - tmp
     */
    for (int j = 0; j < maxRegPerProc; j++) { // step 1

//      Teuchos::RCP<Teuchos::FancyOStream> fos = Teuchos::fancyOStream(Teuchos::rcpFromRef(std::cout));
//      regRes[j]->getMap()->describe(*fos, Teuchos::VERB_EXTREME);
//      regionGrpMats[j]->getRangeMap()->describe(*fos, Teuchos::VERB_EXTREME);

//      TEUCHOS_ASSERT(regionGrpMats[j]->getDomainMap()->isSameAs(*regX[j]->getMap()));
//      TEUCHOS_ASSERT(regionGrpMats[j]->getRangeMap()->isSameAs(*regRes[j]->getMap()));

      regionGrpMats[j]->apply(*regX[j], *regRes[j]);
    }

    sumInterfaceValues(regRes, mapComp, maxRegPerProc, rowMapPerGrp,
        revisedRowMapPerGrp, rowImportPerGrp);

    for (int j = 0; j < maxRegPerProc; j++) { // step 3
      regRes[j]->update(1.0, *regB[j], -1.0);
    }

    // check for convergence
    {
      RCP<Vector> compRes = VectorFactory::Build(mapComp, true);
      regionalToComposite(regRes, compRes, maxRegPerProc, rowMapPerGrp,
          rowImportPerGrp, Xpetra::ADD);
      typename Teuchos::ScalarTraits<Scalar>::magnitudeType normRes = compRes->norm2();

      if (normRes < 1.0e-12)
        return;
    }

    for (int j = 0; j < maxRegPerProc; j++) {
      // update solution according to Jacobi's method
      regX[j]->elementWiseMultiply(omega, *diag[j], *regRes[j], SC_ONE);
    }
  }

  return;
}


//! Recursive V-cycle in region fashion
template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void vCycle(const int l, ///< ID of current level
            const int numLevels, ///< Total number of levels
            const int maxFineIter, ///< max. sweeps on fine and intermediate levels
            const int maxCoarseIter, ///< max. sweeps on coarse level
            const double omega, ///< damping parameter for Jacobi smoother
            const int maxRegPerProc, ///< Max number of regions per process
            std::vector<RCP<Xpetra::Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node> > >& fineRegX, ///< solution
            std::vector<RCP<Xpetra::Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node> > > fineRegB, ///< right hand side
            Array<std::vector<RCP<Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> > > > regMatrices, ///< Matrices in region layout
            Array<std::vector<RCP<Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> > > > regProlong, ///< Prolongators in region layout
            Array<RCP<Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node> > > compRowMaps, ///< composite maps
            Array<std::vector<RCP<Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node> > > > quasiRegRowMaps, ///< quasiRegional row maps
            Array<std::vector<RCP<Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node> > > > regRowMaps, ///< regional row maps
            Array<std::vector<RCP<Xpetra::Import<LocalOrdinal, GlobalOrdinal, Node> > > > regRowImporters, ///< regional row importers
            Array<std::vector<RCP<Xpetra::Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node> > > > regInterfaceScalings, ///< regional interface scaling factors
            RCP<Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> > coarseCompMat ///< Coarsest level composite operator
            )
{
#include "Xpetra_UseShortNames.hpp"
  const Scalar SC_ZERO = Teuchos::ScalarTraits<Scalar>::zero();
  const Scalar SC_ONE = Teuchos::ScalarTraits<Scalar>::one();

  if (l < numLevels - 1) { // fine or intermediate levels

//    std::cout << "level: " << l << std::endl;

    // pre-smoothing
    jacobiIterate(maxFineIter, omega, fineRegX, fineRegB, regMatrices[l],
                  regInterfaceScalings[l], maxRegPerProc, compRowMaps[l],
                  quasiRegRowMaps[l], regRowMaps[l], regRowImporters[l]);

    std::vector<RCP<Vector> > regRes(maxRegPerProc);
    createRegionalVector(regRes, maxRegPerProc, regRowMaps[l]);
    computeResidual(regRes, fineRegX, fineRegB, regMatrices[l], compRowMaps[l],
                    quasiRegRowMaps[l], regRowMaps[l], regRowImporters[l]);

    // Transfer to coarse level
    std::vector<RCP<Vector> > coarseRegX(maxRegPerProc);
    std::vector<RCP<Vector> > coarseRegB(maxRegPerProc);
    for (int j = 0; j < maxRegPerProc; j++) {
      coarseRegX[j] = VectorFactory::Build(regRowMaps[l+1][j], true);
      coarseRegB[j] = VectorFactory::Build(regRowMaps[l+1][j], true);

      regRes[j]->elementWiseMultiply(SC_ONE, *regRes[j], *((regInterfaceScalings[l])[j]), SC_ZERO);

      regProlong[l+1][j]->apply(*regRes[j], *coarseRegB[j], Teuchos::TRANS);
      TEUCHOS_ASSERT(regProlong[l+1][j]->getRangeMap()->isSameAs(*regRes[j]->getMap()));
      TEUCHOS_ASSERT(regProlong[l+1][j]->getDomainMap()->isSameAs(*coarseRegB[j]->getMap()));
    }
    sumInterfaceValues(coarseRegB, compRowMaps[l+1], maxRegPerProc,
                       quasiRegRowMaps[l+1], regRowMaps[l+1], regRowImporters[l+1]);

    // Call V-cycle recursively
    vCycle(l+1, numLevels, maxFineIter, maxCoarseIter, omega, maxRegPerProc,
           coarseRegX, coarseRegB, regMatrices, regProlong, compRowMaps,
           quasiRegRowMaps, regRowMaps, regRowImporters, regInterfaceScalings, coarseCompMat);

    // Transfer coarse level correction to fine level
    std::vector<RCP<Vector> > regCorrection(maxRegPerProc);
    for (int j = 0; j < maxRegPerProc; j++) {
      regCorrection[j] = VectorFactory::Build(regRowMaps[l][j], true);
      regProlong[l+1][j]->apply(*coarseRegX[j], *regCorrection[j]);
      TEUCHOS_ASSERT(regProlong[l+1][j]->getDomainMap()->isSameAs(*coarseRegX[j]->getMap()));
      TEUCHOS_ASSERT(regProlong[l+1][j]->getRangeMap()->isSameAs(*regCorrection[j]->getMap()));
    }

    // apply coarse grid correction
    for (int j = 0; j < maxRegPerProc; j++) {
      fineRegX[j]->update(SC_ONE, *regCorrection[j], SC_ONE);
    }

//    std::cout << "level: " << l << std::endl;

    // post-smoothing
    jacobiIterate(maxFineIter, omega, fineRegX, fineRegB, regMatrices[l],
                  regInterfaceScalings[l], maxRegPerProc, compRowMaps[l], quasiRegRowMaps[l],
                  regRowMaps[l], regRowImporters[l]);
  } else {

//    // create row and column maps with continuous GIDs
//    RCP<Map> contigRowMap = Teuchos::null; //= MapFactory::Build(coarseCompMat->getRowMap());
//    RCP<Map> contigColMap = Teuchos::null; //MapFactory::Build(coarseCompMat->getColMap());
//    createContinuousCoarseLevelMaps(coarseCompMat->getRowMap(),
//                                    coarseCompMat->getColMap(),
//                                    contigRowMap, contigColMap);
//
//    // Store non contiguous maps for later
//    RCP<const Map> noncontigRowMap = coarseCompMat->getRowMap();
//    RCP<const Map> noncontigColMap = coarseCompMat->getColMap();
//
//    // Create composite error vector (zero initial guess)
//    RCP<Vector> compX = VectorFactory::Build(coarseCompMat->getRowMap(), true);
//
//    // Create composite right-hand side vector
//    RCP<Vector> compRhs = VectorFactory::Build(coarseCompMat->getRowMap(), true);
//    {
//      for (int j = 0; j < maxRegPerProc; j++) {
//        RCP<Vector> inverseInterfaceScaling = VectorFactory::Build(regInterfaceScalings[l][j]->getMap());
//        inverseInterfaceScaling->reciprocal(*regInterfaceScalings[l][j]);
//        fineRegB[j]->elementWiseMultiply(SC_ONE, *fineRegB[j], *inverseInterfaceScaling, SC_ZERO);
//      }
//
//      regionalToComposite(fineRegB, compRhs, maxRegPerProc, quasiRegRowMaps[l],
//        regRowImporters[l], Xpetra::ADD);
//    }

    RCP<Teuchos::FancyOStream> fos = Teuchos::fancyOStream(Teuchos::rcpFromRef(std::cout));
    fos->setOutputToRootOnly(0);
    *fos << "+++++++++++++++++++++++++++ WARNING +++++++++++++++++++++++++\n"
         << "+ Coarse level solver has not been migrated to Xpetra, yet. +\n"
         << "+ We skip it for now.                                       +\n"
         << "+++++++++++++++++++++++++++ WARNING +++++++++++++++++++++++++"
        << std::endl;

//    TEUCHOS_TEST_FOR_EXCEPT_MSG(true, "Coarse level solver not migrated to Xpetra, yet.");

//    {
//      Level level;
//      RCP<FactoryManager> factoryHandler = rcp(new FactoryManager());
//      factoryHandler->SetKokkosRefactor(false);
//      level.SetFactoryManager(factoryHandler);
//      level.SetLevelID(0);
//      level.Set("A", coarseCompMat);
//      level.setlib(compX->getMap()->UnderlyingLib());
//    }


    {
      /*

       1. Create DirectSolver by calling its constructor
       2. Call DirectSolver::Copy() to obtain Amesos/Amesos2 object wrapped into a SmootherPrototype
       3. Call Setup() and Apply() on the SmootherPrototype

       */
    }

//    {
//      using DirectSolver = MueLu::DirectSolver<Scalar, LocalOrdinal, GlobalOrdinal, Node>;
//      using FactoryManager = MueLu::FactoryManager<Scalar, LocalOrdinal, GlobalOrdinal, Node>;
//      using Hierarchy = MueLu::Hierarchy<Scalar, LocalOrdinal, GlobalOrdinal, Node>;
//      using SmootherPrototype = MueLu::SmootherPrototype<Scalar, LocalOrdinal, GlobalOrdinal, Node>;
//      using SmootherFactory = MueLu::SmootherFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>;
//
//      RCP<Hierarchy> H = rcp(new Hierarchy(coarseCompMat));
//
//      RCP<SmootherPrototype> coarseProto = rcp(new DirectSolver("Klu"));
//      RCP<SmootherFactory> coarseSolveFact = rcp(new SmootherFactory(coarseProto));
//
//      FactoryManager M;
//      M.SetFactory("CoarseSolver", coarseSolveFact);
//
//      H->Setup(M, 0, 1);
//      H->Iterate(*compRhs, *compX, 1);
//    }

    // Replace non-continuos maps by continuous maps
//    compX->replaceMap(contigRowMap);
//    compRhs->replaceMap(contigRowMap);



//    coarseCompMat->replaceRowMap(*contigRowMap);

    //    TEUCHOS_ASSERT(err==0);
//    err = coarseCompMat->ReplaceColMap(*contigColMap);
//    TEUCHOS_ASSERT(err==0);
//    err = coarseCompMat->ExpertStaticFillComplete(*contigRowMap, *contigRowMap);
//    TEUCHOS_ASSERT(err==0);
//
//    // create a linear problem object
//    Epetra_LinearProblem problem(coarseCompMat.get(), &(*compX), &(*compRhs));
//
//    // Direct solver
//    {
//      Teuchos::ParameterList pList;
//      pList.set("PrintTiming",true);
//      pList.set("PrintStatus",true);
//      pList.set("MaxProcs", coarseCompMat->Comm().NumProc());
//
//      Amesos Factory;
//      Amesos_BaseSolver* solver = Factory.Create("Amesos_Umfpack", problem);
//      TEUCHOS_ASSERT(solver!=NULL);
//
//      solver->SetParameters(pList);
//      solver->SetUseTranspose(false);
//
//      solver->SymbolicFactorization();
//      solver->NumericFactorization();
//      solver->Solve();
//    }
//
//    // Replace maps with the original non-continuous maps
//    compX->ReplaceMap(*noncontigColMap);
//    compRhs->ReplaceMap(*noncontigRowMap);
  }

  return;
}

#endif // MUELU_SETUPREGIONHIERARCHY_DEF_HPP
