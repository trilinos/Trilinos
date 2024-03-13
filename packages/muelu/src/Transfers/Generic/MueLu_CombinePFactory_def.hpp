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
//                    Andrey Prokopenko (aprokop@sandia.gov)
//                    Ray Tuminaro      (rstumin@sandia.gov)
//
// ***********************************************************************
//
// @HEADER
#ifndef MUELU_COMBINEPFACTORY_DEF_HPP
#define MUELU_COMBINEPFACTORY_DEF_HPP

#include <stdlib.h>
#include <iomanip>

// #include <Teuchos_LAPACK.hpp>
#include <Teuchos_SerialDenseMatrix.hpp>
#include <Teuchos_SerialDenseVector.hpp>
#include <Teuchos_SerialDenseSolver.hpp>

#include <Xpetra_CrsMatrixWrap.hpp>
#include <Xpetra_ImportFactory.hpp>
#include <Xpetra_Matrix.hpp>
#include <Xpetra_MapFactory.hpp>
#include <Xpetra_MultiVectorFactory.hpp>
#include <Xpetra_VectorFactory.hpp>

#include <Xpetra_IO.hpp>

#include "MueLu_FactoryManagerBase.hpp"
#include "MueLu_CombinePFactory_decl.hpp"

#include "MueLu_MasterList.hpp"
#include "MueLu_Monitor.hpp"

namespace MueLu {

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
RCP<const ParameterList> CombinePFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::GetValidParameterList() const {
  RCP<ParameterList> validParamList = rcp(new ParameterList());
  validParamList->setEntry("combine: numBlks", ParameterEntry(1));
  validParamList->set<RCP<const FactoryBase>>("A", Teuchos::null, "Generating factory of the matrix A");

  return validParamList;
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void CombinePFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    DeclareInput(Level& fineLevel, Level& /* coarseLevel */) const {
  //    Input(fineLevel, "subblock");
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void CombinePFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Build(Level& fineLevel,
                                                                       Level& coarseLevel) const {
  return BuildP(fineLevel, coarseLevel);
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void CombinePFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::BuildP(Level& fineLevel,
                                                                        Level& coarseLevel) const {
  FactoryMonitor m(*this, "Build", coarseLevel);

  const ParameterList& pL = GetParameterList();
  const LO nBlks          = as<LO>(pL.get<int>("combine: numBlks"));

  RCP<Matrix> A = Get<RCP<Matrix>>(fineLevel, "A");

  // Record all matrices that each define a block in block diagonal comboP
  // matrix used for PDE/multiblock interpolation.  Additionally, count
  // total number of local rows, nonzeros, coarseDofs, and colDofs.

  Teuchos::ArrayRCP<RCP<Matrix>> arrayOfMatrices(nBlks);
  size_t nComboRowMap = 0, nnzCombo = 0, nComboColMap = 0, nComboDomMap = 0;

  LO nTotalNumberLocalColMapEntries = 0;
  Teuchos::ArrayRCP<size_t> DomMapSizePerBlk(nBlks);
  Teuchos::ArrayRCP<size_t> ColMapSizePerBlk(nBlks);
  Teuchos::ArrayRCP<size_t> ColMapLocalSizePerBlk(nBlks);
  Teuchos::ArrayRCP<size_t> ColMapRemoteSizePerBlk(nBlks);
  Teuchos::ArrayRCP<size_t> ColMapLocalCumulativePerBlk(nBlks + 1);   // hardwire 0th entry so that it has the value of 0
  Teuchos::ArrayRCP<size_t> ColMapRemoteCumulativePerBlk(nBlks + 1);  // hardwire 0th entry so that it has the value of 0
  for (int j = 0; j < nBlks; j++) {
    std::string blockName = "Psubblock" + Teuchos::toString(j);
    if (coarseLevel.IsAvailable(blockName, NoFactory::get())) {
      arrayOfMatrices[j] = coarseLevel.Get<RCP<Matrix>>(blockName, NoFactory::get());
      nComboRowMap += Teuchos::as<size_t>((arrayOfMatrices[j])->getRowMap()->getLocalNumElements());
      DomMapSizePerBlk[j] = Teuchos::as<size_t>((arrayOfMatrices[j])->getDomainMap()->getLocalNumElements());
      ColMapSizePerBlk[j] = Teuchos::as<size_t>((arrayOfMatrices[j])->getColMap()->getLocalNumElements());
      nComboDomMap += DomMapSizePerBlk[j];
      nComboColMap += ColMapSizePerBlk[j];
      nnzCombo += Teuchos::as<size_t>((arrayOfMatrices[j])->getLocalNumEntries());
      TEUCHOS_TEST_FOR_EXCEPTION((arrayOfMatrices[j])->getDomainMap()->getIndexBase() != 0, Exceptions::RuntimeError, "interpolation subblocks must use 0 indexbase");

      // figure out how many empty entries in each column map
      int tempii = 0;
      for (int i = 0; i < (int)DomMapSizePerBlk[j]; i++) {
        //          if ( (arrayOfMatrices[j])->getDomainMap()->getGlobalElement(i) == (arrayOfMatrices[j])->getColMap()->getGlobalElement(tempii) )  nTotalNumberLocalColMapEntries++;
        if ((arrayOfMatrices[j])->getDomainMap()->getGlobalElement(i) == (arrayOfMatrices[j])->getColMap()->getGlobalElement(tempii)) tempii++;
      }
      nTotalNumberLocalColMapEntries += tempii;
      ColMapLocalSizePerBlk[j]  = tempii;
      ColMapRemoteSizePerBlk[j] = ColMapSizePerBlk[j] - ColMapLocalSizePerBlk[j];
    } else {
      arrayOfMatrices[j]        = Teuchos::null;
      ColMapLocalSizePerBlk[j]  = 0;
      ColMapRemoteSizePerBlk[j] = 0;
    }
    ColMapLocalCumulativePerBlk[j + 1]  = ColMapLocalSizePerBlk[j] + ColMapLocalCumulativePerBlk[j];
    ColMapRemoteCumulativePerBlk[j + 1] = ColMapRemoteSizePerBlk[j] + ColMapRemoteCumulativePerBlk[j];
  }
  TEUCHOS_TEST_FOR_EXCEPTION(nComboRowMap != A->getRowMap()->getLocalNumElements(), Exceptions::RuntimeError, "sum of subblock rows != #row's Afine");

  // build up csr arrays for combo block diagonal P
  Teuchos::ArrayRCP<size_t> comboPRowPtr(nComboRowMap + 1);
  Teuchos::ArrayRCP<LocalOrdinal> comboPCols(nnzCombo);
  Teuchos::ArrayRCP<Scalar> comboPVals(nnzCombo);

  size_t nnzCnt = 0, nrowCntFromPrevBlks = 0, ncolCntFromPrevBlks = 0;
  LO maxNzPerRow = 0;
  for (int j = 0; j < nBlks; j++) {
    // grab csr pointers for individual blocks of P
    if (arrayOfMatrices[j] != Teuchos::null) {
      Teuchos::ArrayRCP<const size_t> subblockRowPtr((arrayOfMatrices[j])->getLocalNumRows());
      Teuchos::ArrayRCP<const LocalOrdinal> subblockCols((arrayOfMatrices[j])->getLocalNumEntries());
      Teuchos::ArrayRCP<const Scalar> subblockVals((arrayOfMatrices[j])->getLocalNumEntries());
      if ((int)(arrayOfMatrices[j])->getLocalMaxNumRowEntries() > maxNzPerRow) maxNzPerRow = (int)(arrayOfMatrices[j])->getLocalMaxNumRowEntries();
      Teuchos::RCP<CrsMatrixWrap> subblockwrap = Teuchos::rcp_dynamic_cast<CrsMatrixWrap>((arrayOfMatrices[j]));
      Teuchos::RCP<CrsMatrix> subblockcrs      = subblockwrap->getCrsMatrix();
      subblockcrs->getAllValues(subblockRowPtr, subblockCols, subblockVals);

      // copy jth block into csr arrays of comboP

      for (decltype(subblockRowPtr.size()) i = 0; i < subblockRowPtr.size() - 1; i++) {
        size_t rowLength                      = subblockRowPtr[i + 1] - subblockRowPtr[i];
        comboPRowPtr[nrowCntFromPrevBlks + i] = nnzCnt;
        for (size_t k = 0; k < rowLength; k++) {
          if ((int)subblockCols[k + subblockRowPtr[i]] < (int)ColMapLocalSizePerBlk[j]) {
            comboPCols[nnzCnt] = subblockCols[k + subblockRowPtr[i]] + ColMapLocalCumulativePerBlk[j];
            if ((int)comboPCols[nnzCnt] >= (int)ColMapLocalCumulativePerBlk[nBlks]) {
              printf("ERROR1\n");
              exit(1);
            }
          } else {
            // Here we subtract off the number of local colmap guys ... so this tell us where we are among ghost unknowns. We then want to stick this ghost after
            // all the Local guys ... so we add ColMapLocalCumulativePerBlk[nBlks] .... finally we need to account for the fact that other ghost blocks may have already
            // been handled ... so we then add  + ColMapRemoteCumulativePerBlk[j];
            comboPCols[nnzCnt] = subblockCols[k + subblockRowPtr[i]] - ColMapLocalSizePerBlk[j] + ColMapLocalCumulativePerBlk[nBlks] + ColMapRemoteCumulativePerBlk[j];
            if ((int)comboPCols[nnzCnt] >= (int)(ColMapLocalCumulativePerBlk[nBlks] + ColMapRemoteCumulativePerBlk[nBlks])) {
              printf("ERROR2\n");
              exit(1);
            }
          }
          comboPVals[nnzCnt++] = subblockVals[k + subblockRowPtr[i]];
        }
      }

      nrowCntFromPrevBlks += Teuchos::as<size_t>((arrayOfMatrices[j])->getRowMap()->getLocalNumElements());
      ncolCntFromPrevBlks += DomMapSizePerBlk[j];  // rst: check this
    }
  }
  comboPRowPtr[nComboRowMap] = nnzCnt;

  // Come up with global IDs for the coarse grid maps. We assume that each xxx
  // block has a minimum GID of 0.  Since MueLu is generally creating these
  // GIDS, this assumption is probably correct, but we'll check it.

  Teuchos::Array<GlobalOrdinal> comboDomainMapGIDs(nComboDomMap);
  Teuchos::Array<GlobalOrdinal> comboColMapGIDs(nComboColMap);

  LO nTotalNumberRemoteColMapEntries = 0;
  GlobalOrdinal offset               = 0;
  size_t domainMapIndex              = 0;
  int nComboColIndex                 = 0;
  for (int j = 0; j < nBlks; j++) {
    int nThisBlkColIndex = 0;

    GlobalOrdinal tempMax = 0, maxGID = 0;
    if (arrayOfMatrices[j] != Teuchos::null) tempMax = (arrayOfMatrices[j])->getDomainMap()->getMaxGlobalIndex();
    Teuchos::reduceAll(*(A->getDomainMap()->getComm()), Teuchos::REDUCE_MAX, tempMax, Teuchos::ptr(&maxGID));

    if (arrayOfMatrices[j] != Teuchos::null) {
      TEUCHOS_TEST_FOR_EXCEPTION(arrayOfMatrices[j]->getDomainMap()->getMinAllGlobalIndex() < 0, Exceptions::RuntimeError, "Global ID assumption for domainMap not met within subblock");

      GO priorDomGID = 0;
      for (size_t c = 0; c < DomMapSizePerBlk[j]; ++c) {  // check this
                                                          //  The global ids of jth block are assumed to go between 0 and maxGID_j.  So the 1th blocks DomainGIDs should start at maxGID_0+1. The 2nd
                                                          //  block DomainDIGS starts at maxGID_0+1 + maxGID_1 + 1. We use offset to keep track of these starting points.
        priorDomGID                          = (arrayOfMatrices[j])->getDomainMap()->getGlobalElement(c);
        comboDomainMapGIDs[domainMapIndex++] = offset + priorDomGID;
        if (priorDomGID == (arrayOfMatrices[j])->getColMap()->getGlobalElement(nThisBlkColIndex)) {
          comboColMapGIDs[nComboColIndex++] = offset + priorDomGID;
          nThisBlkColIndex++;
        }
      }

      for (size_t cc = nThisBlkColIndex; cc < ColMapSizePerBlk[j]; ++cc) {
        comboColMapGIDs[nTotalNumberLocalColMapEntries + nTotalNumberRemoteColMapEntries++] = offset + (arrayOfMatrices[j])->getColMap()->getGlobalElement(cc);
      }
    }
    offset += maxGID + 1;
  }
#ifdef out
  RCP<const Map> coarseDomainMap = Xpetra::MapFactory<LO, GO, NO>::Build(A->getDomainMap()->lib(), Teuchos::OrdinalTraits<Xpetra::global_size_t>::invalid(), comboDomainMapGIDs, 0, A->getDomainMap()->getComm());
  RCP<const Map> coarseColMap    = Xpetra::MapFactory<LO, GO, NO>::Build(A->getDomainMap()->lib(), Teuchos::OrdinalTraits<Xpetra::global_size_t>::invalid(), comboColMapGIDs, 0, A->getDomainMap()->getComm());
#endif

  RCP<const Map> coarseDomainMap = Xpetra::MapFactory<LO, GO, NO>::Build(A->getDomainMap()->lib(), Teuchos::OrdinalTraits<Xpetra::global_size_t>::invalid(), comboDomainMapGIDs, 0, A->getRowMap()->getComm());
  RCP<const Map> coarseColMap    = Xpetra::MapFactory<LO, GO, NO>::Build(A->getDomainMap()->lib(), Teuchos::OrdinalTraits<Xpetra::global_size_t>::invalid(), comboColMapGIDs, 0, A->getRowMap()->getComm());

  Teuchos::RCP<CrsMatrix> comboPCrs = CrsMatrixFactory::Build(A->getRowMap(), coarseColMap, maxNzPerRow);
  // comboPCrs->getCrsGraph(); //.getRowInfo(6122);
  // comboPCrs->getRowInfo(6122);

  //    Teuchos::RCP<CrsMatrix> comboPCrs = CrsMatrixFactory::Build(A->getRowMap(), coarseColMap,nnzCombo+1000);

  //    for (size_t i = 0; i < nComboRowMap; i++) {
  // printf("FIXME\n"); if (nComboRowMap > 6142)  nComboRowMap = 6142;
  for (size_t i = 0; i < nComboRowMap; i++) {
    comboPCrs->insertLocalValues(i, comboPCols.view(comboPRowPtr[i], comboPRowPtr[i + 1] - comboPRowPtr[i]),
                                 comboPVals.view(comboPRowPtr[i], comboPRowPtr[i + 1] - comboPRowPtr[i]));
  }
  comboPCrs->fillComplete(coarseDomainMap, A->getRowMap());

  Teuchos::RCP<Matrix> comboP = Teuchos::rcp(new CrsMatrixWrap(comboPCrs));

  Set(coarseLevel, "P", comboP);
}

}  // namespace MueLu

#define MUELU_COMBINEPFACTORY_SHORT
#endif  // MUELU_COMBINEPFACTORY_DEF_HPP
