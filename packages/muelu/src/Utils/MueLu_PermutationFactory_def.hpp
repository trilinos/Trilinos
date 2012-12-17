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
//                    Jeremie Gaidamour (jngaida@sandia.gov)
//                    Jonathan Hu       (jhu@sandia.gov)
//                    Ray Tuminaro      (rstumin@sandia.gov)
//
// ***********************************************************************
//
// @HEADER
/*
 * MueLu_PermutationFactory_def.hpp
 *
 *  Created on: Nov 28, 2012
 *      Author: wiesner
 */

#ifndef MUELU_PERMUTATIONFACTORY_DEF_HPP_
#define MUELU_PERMUTATIONFACTORY_DEF_HPP_

#include <vector>
//#include <pair>

#include "MueLu_PermutationFactory_decl.hpp"

#include <Xpetra_Map.hpp>
#include <Xpetra_Vector.hpp>
#include <Xpetra_VectorFactory.hpp>
#include <Xpetra_Matrix.hpp>
#include <Xpetra_CrsMatrixWrap.hpp>
#include <Xpetra_Export.hpp>
#include <Xpetra_ExportFactory.hpp>
#include <Xpetra_Import.hpp>
#include <Xpetra_ImportFactory.hpp>

#include "MueLu_Level.hpp"
#include "MueLu_Utilities.hpp"
// #include "MueLu_Monitor.hpp"

// MPI helper
#define sumAll(rcpComm, in, out)                                        \
    Teuchos::reduceAll(*rcpComm, Teuchos::REDUCE_SUM, in, Teuchos::outArg(out));
#define minAll(rcpComm, in, out)                                        \
    Teuchos::reduceAll(*rcpComm, Teuchos::REDUCE_MIN, in, Teuchos::outArg(out));
#define maxAll(rcpComm, in, out)                                        \
    Teuchos::reduceAll(*rcpComm, Teuchos::REDUCE_MAX, in, Teuchos::outArg(out));

namespace MueLu {
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
PermutationFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::PermutationFactory(std::string const & mapName, const RCP<const FactoryBase> & mapFact)
: mapName_(mapName), mapFact_(mapFact)
  { }

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
PermutationFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::~PermutationFactory() {}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
void PermutationFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::DeclareInput(Level &currentLevel) const {
  Input(currentLevel, "A");
  if(mapName_.length() > 0 ) {
    currentLevel.DeclareInput(mapName_,mapFact_.get(),this);
  }
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
void PermutationFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::Build(Level & currentLevel) const {

  Teuchos::RCP<Matrix> A = Get< Teuchos::RCP<Matrix> > (currentLevel, "A");

  Teuchos::RCP<const Map> permRowMap = Teuchos::null;
  if(mapName_.length() > 0 ) {
    permRowMap = currentLevel.Get<RCP<const Map> >(mapName_,mapFact_.get());
  } else {
    permRowMap = A->getRowMap(); // use full row map of A
  }

  GetOStream(Runtime0, 0) << "Perform generation of permutation operators on " << mapName_ << " map with " << permRowMap->getGlobalNumElements() << " elements" << std::endl;

  std::vector<std::pair<GlobalOrdinal, GlobalOrdinal> > permutedDiagCandidates;
  std::vector<std::pair<GlobalOrdinal, GlobalOrdinal> > keepDiagonalEntries;
  std::vector<Scalar> Weights;

  // loop over all local rows in matrix A and keep diagonal entries if corresponding
  // matrix rows are not contained in permRowMap
  for (size_t row = 0; row < A->getRowMap()->getNodeNumElements(); row++) {
    GlobalOrdinal grow = A->getRowMap()->getGlobalElement(row);

    if(permRowMap->isNodeGlobalElement(grow) == true) continue;

    size_t nnz = A->getNumEntriesInLocalRow(row);

    // extract local row information from matrix
    Teuchos::ArrayView<const LocalOrdinal> indices;
    Teuchos::ArrayView<const Scalar> vals;
    A->getLocalRowView(row, indices, vals);

    TEUCHOS_TEST_FOR_EXCEPTION(Teuchos::as<size_t>(indices.size()) != nnz, Exceptions::RuntimeError, "MueLu::PermutationFactory::Build: number of nonzeros not equal to number of indices? Error.");

    // find column entry with max absolute value
    GlobalOrdinal gMaxValIdx = 0;
    Scalar norm1 = 0.0;
    Scalar maxVal = 0.0;
    for (size_t j = 0; j < indices.size(); j++) {
      norm1 += std::abs(vals[j]);
      if(std::abs(vals[j]) > maxVal) {
        maxVal = std::abs(vals[j]);
        gMaxValIdx = A->getColMap()->getGlobalElement(indices[j]);
      }
    }

    if(grow == gMaxValIdx) // only keep row/col pair if it's diagonal dominant!!!
      keepDiagonalEntries.push_back(std::make_pair(grow,grow));

  }

  //////////
  // handle rows that are marked to be relevant for permutations
  for (size_t row = 0; row < permRowMap->getNodeNumElements(); row++) {
    GlobalOrdinal grow = permRowMap->getGlobalElement(row);
    LocalOrdinal lArow = A->getRowMap()->getLocalElement(grow);
    size_t nnz = A->getNumEntriesInLocalRow(lArow);

    // extract local row information from matrix
    Teuchos::ArrayView<const LocalOrdinal> indices;
    Teuchos::ArrayView<const Scalar> vals;
    A->getLocalRowView(lArow, indices, vals);

    TEUCHOS_TEST_FOR_EXCEPTION(Teuchos::as<size_t>(indices.size()) != nnz, Exceptions::RuntimeError, "MueLu::PermutationFactory::Build: number of nonzeros not equal to number of indices? Error.");

    // find column entry with max absolute value
    GlobalOrdinal gMaxValIdx = 0;
    Scalar norm1 = 0.0;
    Scalar maxVal = 0.0;
    for (size_t j = 0; j < indices.size(); j++) {
      norm1 += std::abs(vals[j]);
      if(std::abs(vals[j]) > maxVal) {
        maxVal = std::abs(vals[j]);
        gMaxValIdx = A->getColMap()->getGlobalElement(indices[j]);
      }
    }

    if(std::abs(maxVal) > 0.0) { // keep only max Entries \neq 0.0
      permutedDiagCandidates.push_back(std::make_pair(grow,gMaxValIdx));
      Weights.push_back(maxVal/(norm1*Teuchos::as<Scalar>(nnz)));
    }

  }

  // sort Weights in descending order
  std::vector<int> permutation;
  sortingPermutation(Weights,permutation);

  /*  typedef std::vector<int>::const_iterator I;
  for (I p = permutation.begin(); p != permutation.end(); ++p)
    std::cout << *p << " ";
  std::cout << "\n";*/

  // create new vector with exactly one possible entry for each column

  // todo remove all entries with columns that are already in keepDiagonalEntries
  Teuchos::RCP<Vector> gColVec = VectorFactory::Build(A->getColMap());
  Teuchos::RCP<Vector> gDomVec = VectorFactory::Build(A->getDomainMap());
  gColVec->putScalar(0.0);
  gDomVec->putScalar(0.0);

  // put in all keep diagonal entries
  //typename std::vector<std::pair<GlobalOrdinal, GlobalOrdinal> >::const_iterator it;
  for (typename std::vector<std::pair<GlobalOrdinal, GlobalOrdinal> >::const_iterator p = keepDiagonalEntries.begin(); p != keepDiagonalEntries.end(); ++p) {
    gColVec->sumIntoGlobalValue((*p).second,1.0);
  }

  Teuchos::RCP<Export> exporter = ExportFactory::Build(gColVec->getMap(), gDomVec->getMap());
  gDomVec->doExport(*gColVec,*exporter,Xpetra::ADD);

  std::vector<std::pair<GlobalOrdinal, GlobalOrdinal> > permutedDiagCandidatesFiltered; // TODO reserve memory
  std::map<GlobalOrdinal, Scalar> gColId2Weight;

  for(size_t i = 0; i < permutedDiagCandidates.size(); ++i) {
    // loop over all candidates
    std::pair<GlobalOrdinal, GlobalOrdinal> pp = permutedDiagCandidates[permutation[i]];
    GlobalOrdinal grow = pp.first;
    GlobalOrdinal gcol = pp.second;

    LocalOrdinal lcol = A->getColMap()->getLocalElement(gcol);
    Teuchos::ArrayRCP< const Scalar > ddata = gColVec->getData(0);
    if(ddata[lcol] > 0.0){
      continue; // column already handled by another row
    }

    //gColVec->sumIntoGlobalValue(gcol,1.0); // mark column as already taken
    gColVec->sumIntoLocalValue(A->getColMap()->getLocalElement(gcol),1.0);

    permutedDiagCandidatesFiltered.push_back(std::make_pair(grow,gcol));
    gColId2Weight[gcol] = Weights[permutation[i]];
  }

  // communicate how often each column index is requested by the different procs
  gDomVec->doExport(*gColVec,*exporter,Xpetra::ADD);

  //*****************************************************************************************
  // first communicate ALL global ids of column indices which are requested by more
  // than one proc to all other procs
  // detect which global col indices are requested by more than one proc
  // and store them in multipleColRequests
  std::vector<GlobalOrdinal> multipleColRequests; // store all global column indices from current processor that are alos
  // requested by another processor. This is possible, since they are stored
  // in gDomVec which is based on the nonoverlapping domain map. That is, each
  // global col id is handled by exactly one proc.
  for(size_t sz = 0; sz<gDomVec->getLocalLength(); ++sz) {
    Teuchos::ArrayRCP< const Scalar > arrDomVec = gDomVec->getData(0);
    if(arrDomVec[sz] > 1.0) {
      std::cout << "ERROR: arrDomVec[" << sz << "]=" << arrDomVec[sz] << std::endl;
      multipleColRequests.push_back(gDomVec->getMap()->getGlobalElement(sz));
    }
  }

  // communicate the global number of column indices which are requested by more than one proc
  LocalOrdinal localMultColRequests  = Teuchos::as<LocalOrdinal>(multipleColRequests.size());
  LocalOrdinal globalMultColRequests = 0;

  // sum up all entries in multipleColRequests over all processors
  sumAll(gDomVec->getMap()->getComm(), (LocalOrdinal)localMultColRequests, globalMultColRequests);

  if(globalMultColRequests > 0) {
    const Teuchos::RCP< const Teuchos::Comm< int > > comm = gDomVec->getMap()->getComm();
    // communicate how many interesting column indices are stored on each proc
    int numProcs = comm->getSize();
    int myRank   = comm->getRank();

    // distribute number of multipleColRequests to all processors
    // each processor stores how many column ids for exchange are handled by the cur proc
    std::vector<GlobalOrdinal> numMyMultColRequests(numProcs,0);
    std::vector<GlobalOrdinal> numGlobalMultColRequests(numProcs,0);
    numMyMultColRequests[myRank] = localMultColRequests;
    Teuchos::reduceAll<int,int>(*comm,Teuchos::REDUCE_MAX,numProcs,&numMyMultColRequests[0],&numGlobalMultColRequests[0]);

    // communicate multipleColRequests entries to all processors
    int nMyOffset = 0;
    for (int i=0; i<myRank-1; i++)
      nMyOffset += numGlobalMultColRequests[i]; // calculate offset to store the weights on the corresponding place in procOverlappingWeights

    std::vector<GlobalOrdinal> procMultRequestedColIds(globalMultColRequests,0.0);
    std::vector<GlobalOrdinal> global_procMultRequestedColIds(globalMultColRequests,0.0);

    // loop over all local column GIDs that are also requested by other procs
    for(size_t i = 0; i < multipleColRequests.size(); i++) {
      procMultRequestedColIds[nMyOffset + i] = multipleColRequests[i]; // all weights are > 0 ?
    }

    // template ordinal, package (double)
    Teuchos::reduceAll<int,GlobalOrdinal>(*comm, Teuchos::REDUCE_MAX, (int) globalMultColRequests,&procMultRequestedColIds[0],&global_procMultRequestedColIds[0]);

    // loop over global_procOverlappingWeights and eliminate wrong entries...
    for (size_t k = 0; k<global_procMultRequestedColIds.size(); k++) {
      GlobalOrdinal globColId = global_procMultRequestedColIds[k];

      std::vector<GlobalOrdinal> MyWeightForColId(numProcs,0);
      std::vector<GlobalOrdinal> GlobalWeightForColId(numProcs,0);

      if(gColVec->getMap()->isNodeGlobalElement(globColId)) {
        MyWeightForColId[myRank] = gColId2Weight[globColId];
      } else {
        MyWeightForColId[myRank] = 0.0;
      }

      Teuchos::reduceAll<int,int>(*comm,Teuchos::REDUCE_MAX,numProcs,&MyWeightForColId[0],&GlobalWeightForColId[0]);

      if(gColVec->getMap()->isNodeGlobalElement(globColId)) {
        // note: 2 procs could have the same weight for a column index.
        // pick the first one.
        Scalar winnerValue = 0.0;
        int winnerProcRank = 0;
        for (int proc = 0; proc < numProcs; proc++) {
          if(GlobalWeightForColId[proc] > winnerValue) {
            winnerValue = GlobalWeightForColId[proc];
            winnerProcRank = proc;
          }
        }

        // winnerProcRank is the winner for handling globColId.
        // winnerProcRank is unique (even if two procs have the same weight for a column index)

        if(myRank != winnerProcRank) {
          // remove corresponding entry from permutedDiagCandidatesFiltered
          typename std::vector<std::pair<GlobalOrdinal, GlobalOrdinal> >::iterator p = permutedDiagCandidatesFiltered.begin();
          while(p != permutedDiagCandidatesFiltered.end() )
          {
            GlobalOrdinal jk = (*p).second;
            if((*p).second == globColId)
              p = permutedDiagCandidatesFiltered.erase(p);
            else
              p++;
          }
        }

      } // end if isNodeGlobalElement
    } // end loop over global_procOverlappingWeights and eliminate wrong entries...
  } // end if globalMultColRequests > 0

  // put together all pairs:
  size_t sizeRowColPairs = keepDiagonalEntries.size() + permutedDiagCandidatesFiltered.size();
  std::vector<std::pair<GlobalOrdinal, GlobalOrdinal> > RowColPairs;
  RowColPairs.insert( RowColPairs.end(), keepDiagonalEntries.begin(), keepDiagonalEntries.end());
  RowColPairs.insert( RowColPairs.end(), permutedDiagCandidatesFiltered.begin(), permutedDiagCandidatesFiltered.end());

  //*****************************************************************************************

  //    for (typename std::vector<std::pair<GlobalOrdinal, GlobalOrdinal> >::const_iterator p = RowColPairs.begin(); p != RowColPairs.end(); ++p)
  //    {
  ////      if((*p).first != (*p).second) std::cout << "difference: " << (*p).first << " " << (*p).second << std::endl;
  //      std::cout << (*p).first +1 << " " << (*p).second+1 << std::endl;
  //    }
  //    std::cout << "\n";

  //////////////////////////////////////////////////
  // assumption: on each processor RowColPairs now contains
  // a valid set of (row,column) pairs, where the row entries
  // are a subset of the processor's rows and the column entries
  // are unique throughout all processors.
  // Note: the RowColPairs are only defined for a subset of all rows,
  // so there might be rows without an entry in RowColPairs.
  // It can be, that some rows seem to be missing in RowColPairs, since
  // the entry in that row with maximum absolute value has been reserved
  // by another row already (e.g. as already diagonal dominant row outside
  // of perRowMap).
  // In fact, the RowColPairs vector only defines the (row,column) pairs
  // that will be definitely moved to the diagonal after permutation.

  // build Pperm and Qperm vectors
  Teuchos::RCP<Vector> Pperm = VectorFactory::Build(A->getRowMap());
  Teuchos::RCP<Vector> Qperm = VectorFactory::Build(A->getRowMap());
  Teuchos::RCP<Vector> gRowIndices = VectorFactory::Build(A->getRowMap()); // vector stores, whether a specific row index is already used

  Pperm->putScalar(-1.0);
  Qperm->putScalar(-1.0);
  gRowIndices->putScalar(0.0); // 0: row index can be used, 1: row index is already used

  // loop over local list of future diagonal entries
  Teuchos::ArrayRCP< Scalar > lPpermData = Pperm->getDataNonConst(0);
  Teuchos::ArrayRCP< Scalar > lQpermData = Qperm->getDataNonConst(0);
  Teuchos::ArrayRCP< Scalar > lRowIndicesData = gRowIndices->getDataNonConst(0);

  typename std::vector<std::pair<GlobalOrdinal, GlobalOrdinal> >::iterator p = RowColPairs.begin();
  while(p != RowColPairs.end() )
  {
    GlobalOrdinal ik = (*p).first;
    GlobalOrdinal jk = (*p).second;

    // plausibility check: the jk have to be unique over all procs
    if(lQpermData[A->getColMap()->getLocalElement(jk)] != -1) {
      std::cout << "Error2: Qperm cannot have entries != -1. It seems that the filtered column entries are not unique!" << std::endl;
      std::cout << "ik=" << ik << " jk=" << jk << " lcolidx=" << A->getColMap()->getLocalElement(jk) << " lQpermData " << lQpermData[A->getColMap()->getLocalElement(jk)] << std::endl;
    }


    if(lRowIndicesData[A->getRowMap()->getLocalElement(ik)] == 0.0) {
      lRowIndicesData[A->getRowMap()->getLocalElement(ik)] = 1.0; // use this row id
      lPpermData[A->getRowMap()->getLocalElement(ik)] = ik; // i.e. no row permutation
      lQpermData[A->getColMap()->getLocalElement(jk)] = ik; // only col permutation
      p = RowColPairs.erase(p);
    }
    else
      p++;
  }

  // handle remaining permutations
  LocalOrdinal lNextRowIndex = 0; // remove me???
  LocalOrdinal lCurRowIndex = 0;
  p = RowColPairs.begin();
  while(p != RowColPairs.end() )
  {
    GlobalOrdinal ik = (*p).first;
    GlobalOrdinal jk = (*p).second;

    // plausibility check: the jk have to be unique over all procs
    if(lQpermData[A->getColMap()->getLocalElement(jk)] != -1) {
      std::cout << "Error2: Qperm cannot have entries != -1. It seems that the filtered column entries are not unique!" << std::endl;
      std::cout << "ik=" << ik << " jk=" << jk << " lcolidx=" << A->getColMap()->getLocalElement(jk) << " lQpermData " << lQpermData[A->getColMap()->getLocalElement(jk)] << std::endl;
    }

    // search for first lRowIndicesData neq 0.0
    for(LocalOrdinal l = lNextRowIndex; l<Teuchos::as<LocalOrdinal>(lRowIndicesData.size()); l++) {
      if(lRowIndicesData[l]==0.0) {
        lCurRowIndex  = l;
        lNextRowIndex = l+1;
        break;
      }
    }

    if(lRowIndicesData[lCurRowIndex] == 0.0) {
      lRowIndicesData [lCurRowIndex] = 1.0; // use this row id
      lPpermData[A->getRowMap()->getLocalElement(ik)] = gRowIndices->getMap()->getGlobalElement(lCurRowIndex);
      lQpermData[A->getColMap()->getLocalElement(jk)] = gRowIndices->getMap()->getGlobalElement(lCurRowIndex);
      p = RowColPairs.erase(p);
    }
    else
      p++;
  }

  //  for (typename std::vector<std::pair<GlobalOrdinal, GlobalOrdinal> >::const_iterator p = RowColPairs.begin(); p != RowColPairs.end(); ++p)
  //  {
  //    std::cout << "r/c: " << (*p).first << "/" << (*p).second << std::endl;
  //  }
  //  std::cout << "\n";
  //  std::cout << "Pperm" << *Pperm << std::endl;
  //  std::cout << "Qperm" << *Qperm << std::endl;

  // "copy" RowIndices vector
  // use gRowIndices information for completing Pperm
  // use gColIndices information for completing Qperm
  Teuchos::RCP<Vector> gColIndices = VectorFactory::Build(A->getRowMap()); // vector stores, whether a specific col index is already used
  gColIndices->update(1.0,*gRowIndices,0.0);

  Teuchos::ArrayRCP< Scalar > lColIndicesData = gColIndices->getDataNonConst(0);
  lNextRowIndex = 0; // remove me???
  lCurRowIndex = 0;
  LocalOrdinal lNextColIndex = 0;
  LocalOrdinal lCurColIndex = 0;
  for(LocalOrdinal i=0; i<Teuchos::as<LocalOrdinal>(A->getNodeNumRows()); i++) {
    if(lPpermData[i] < 0.0) {
      // search for first lRowIndicesData neq 0.0
      for(LocalOrdinal l = lNextRowIndex; l<Teuchos::as<LocalOrdinal>(lRowIndicesData.size()); l++) {
        if(lRowIndicesData[l]==0.0) {
          lCurRowIndex  = l;
          lNextRowIndex = l+1;
          break;
        }
      }
      lRowIndicesData [lCurRowIndex] = 1.0; // use this row id
      lPpermData[ i ] = gRowIndices->getMap()->getGlobalElement(lCurRowIndex);
    }
    if(lQpermData[i] < 0.0) {
      // search for first lColIndicesData neq 0.0
      for(LocalOrdinal l = lNextColIndex; l<Teuchos::as<LocalOrdinal>(lColIndicesData.size()); l++) {
        if(lColIndicesData[l]==0.0) {
          lCurColIndex  = l;
          lNextColIndex = l+1;
          break;
        }
      }
      lColIndicesData [lCurColIndex] = 1.0; // use this col id
      lQpermData[ i ] = gColIndices->getMap()->getGlobalElement(lCurColIndex);
    }
  }

  //  std::cout << "Pperm after complete" << *Pperm << std::endl;
  //  std::cout << "Qperm after complete" << *Qperm << std::endl;

  // build permutation matrices

  // create new empty Matrix
  Teuchos::RCP<CrsMatrixWrap> permPTmatrix = Teuchos::rcp(new CrsMatrixWrap(A->getRowMap(),1,Xpetra::StaticProfile));
  Teuchos::RCP<CrsMatrixWrap> permQTmatrix = Teuchos::rcp(new CrsMatrixWrap(A->getRowMap(),1,Xpetra::StaticProfile));

  for(size_t row=0; row<A->getNodeNumRows(); row++) {
    Teuchos::ArrayRCP<GlobalOrdinal> indoutP(1,lPpermData[row]); // column idx for Perm^T
    Teuchos::ArrayRCP<GlobalOrdinal> indoutQ(1,lQpermData[row]); // column idx for Qperm
    Teuchos::ArrayRCP<Scalar> valout(1,1.0);
    permPTmatrix->insertGlobalValues(A->getRowMap()->getGlobalElement(row), indoutP.view(0,indoutP.size()), valout.view(0,valout.size()));
    permQTmatrix->insertGlobalValues (A->getRowMap()->getGlobalElement(row), indoutQ.view(0,indoutQ.size()), valout.view(0,valout.size()));
  }

  permPTmatrix->fillComplete();
  permQTmatrix->fillComplete();

  Teuchos::RCP<Matrix> permPmatrix = Utils2::Transpose(permPTmatrix,true);

  for(size_t row=0; row<permPTmatrix->getNodeNumRows(); row++) {
    if(permPTmatrix->getNumEntriesInLocalRow(row) != 1)
      std::cout<<"#entries in row " << row << " of permPTmatrix is " << permPTmatrix->getNumEntriesInLocalRow(row) << std::endl;
    if(permPmatrix->getNumEntriesInLocalRow(row) != 1)
      std::cout<<"#entries in row " << row << " of permPmatrix is " << permPmatrix->getNumEntriesInLocalRow(row) << std::endl;
    if(permQTmatrix->getNumEntriesInLocalRow(row) != 1)
      std::cout<<"#entries in row " << row << " of permQmatrix is " << permQTmatrix->getNumEntriesInLocalRow(row) << std::endl;
  }

  // build permP * A * permQT
  Teuchos::RCP<Matrix> ApermQt = Utils::TwoMatrixMultiply(A, false, permQTmatrix, false);
  Teuchos::RCP<Matrix> permPApermQt = Utils::TwoMatrixMultiply(Teuchos::rcp_dynamic_cast<Matrix>(permPmatrix), false, ApermQt, false);

  //bool doFillComplete=false;
  //bool optimizeStorage=false;
  //Teuchos::ArrayRCP<SC> diag = Utils::GetMatrixDiagonal(*A);
  //Utils::MyOldScaleMatrix(AP, diag, true, doFillComplete, optimizeStorage); //scale matrix with reciprocal of diag

  /*
  MueLu::Utils<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::Write("A.mat", *A);
  MueLu::Utils<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::Write("permP.mat", *permPmatrix);
  MueLu::Utils<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::Write("permQt.mat", *permQTmatrix);
  MueLu::Utils<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::Write("permPApermQt.mat", *permPApermQt);
   */
  // build scaling matrix
  Teuchos::RCP<Vector> diagVec = VectorFactory::Build(permPApermQt->getRowMap(),true);
  Teuchos::RCP<Vector> invDiagVec = VectorFactory::Build(permPApermQt->getRowMap(),true);
  Teuchos::ArrayRCP< const Scalar > diagVecData = diagVec->getData(0);
  Teuchos::ArrayRCP< Scalar > invDiagVecData = invDiagVec->getDataNonConst(0);

  permPApermQt->getLocalDiagCopy(*diagVec);
  for(size_t i = 0; i<diagVec->getMap()->getNodeNumElements(); ++i) {
    if(diagVecData[i] != 0.0)
      invDiagVecData[i] = 1/diagVecData[i];
    else {
      invDiagVecData[i] = 1.0;
      std::cout << "found zero on diagonal in row " << i << std::endl;
    }
  }

  Teuchos::RCP<CrsMatrixWrap> diagScalingOp = Teuchos::rcp(new CrsMatrixWrap(permPApermQt->getRowMap(),1,Xpetra::StaticProfile));

  for(size_t row=0; row<A->getNodeNumRows(); row++) {
    Teuchos::ArrayRCP<GlobalOrdinal> indout(1,permPApermQt->getRowMap()->getGlobalElement(row)); // column idx for Perm^T
    Teuchos::ArrayRCP<Scalar> valout(1,invDiagVecData[row]);
    diagScalingOp->insertGlobalValues(A->getRowMap()->getGlobalElement(row), indout.view(0,indout.size()), valout.view(0,valout.size()));
  }
  diagScalingOp->fillComplete();

  Teuchos::RCP<Matrix> scaledA = Utils::TwoMatrixMultiply(diagScalingOp, false, permPApermQt, false);
  currentLevel.Set("A", Teuchos::rcp_dynamic_cast<Matrix>(scaledA), this);

  currentLevel.Set("permA", Teuchos::rcp_dynamic_cast<Matrix>(permPApermQt), this);  // TODO careful with this!!!
  currentLevel.Set("permP", Teuchos::rcp_dynamic_cast<Matrix>(permPmatrix), this);
  currentLevel.Set("permQT", Teuchos::rcp_dynamic_cast<Matrix>(permQTmatrix), this);
  currentLevel.Set("permScaling", Teuchos::rcp_dynamic_cast<Matrix>(diagScalingOp), this);
}

} // namespace MueLu


#endif /* MUELU_PERMUTATIONFACTORY_DEF_HPP_ */
