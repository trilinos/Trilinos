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
#include <queue>

#include "MueLu_PermutationFactory_decl.hpp"

#include <Xpetra_Map.hpp>
#include <Xpetra_StridedMap.hpp>    // for nDofsPerNode...
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
#include "MueLu_Monitor.hpp"

#define DEBUG_OUTPUT 1

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
  FactoryMonitor m(*this, "Permutation Factory ", currentLevel);

#ifndef HAVE_MUELU_INST_COMPLEX_INT_INT

  Teuchos::RCP<Matrix> A = Get< Teuchos::RCP<Matrix> > (currentLevel, "A");

  const Teuchos::RCP< const Teuchos::Comm< int > > comm = A->getRowMap()->getComm();
  int numProcs = comm->getSize();
  int myRank   = comm->getRank();


  Teuchos::RCP<const Map> permRowMap = Teuchos::null;
  if(mapName_.length() > 0 ) {
    permRowMap = currentLevel.Get<RCP<const Map> >(mapName_,mapFact_.get());
  } else {
    permRowMap = A->getRowMap(); // use full row map of A
  }

  size_t nDofsPerNode = 1;
  if (A->IsView("stridedMaps")) {
    Teuchos::RCP<const Map> permRowMapStrided = A->getRowMap("stridedMaps");
    nDofsPerNode = Teuchos::rcp_dynamic_cast<const StridedMap>(permRowMapStrided)->getFixedBlockSize();
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
    for (size_t j = 0; j < Teuchos::as<size_t>(indices.size()); j++) {
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
    for (size_t j = 0; j < Teuchos::as<size_t>(indices.size()); j++) {
      norm1 += std::abs(vals[j]);
      if(std::abs(vals[j]) > maxVal) {
        maxVal = std::abs(vals[j]);
        gMaxValIdx = A->getColMap()->getGlobalElement(indices[j]);
      }
    }

    if(std::abs(maxVal) > 0.0) { // keep only max Entries \neq 0.0
      permutedDiagCandidates.push_back(std::make_pair(grow,gMaxValIdx));
      Weights.push_back(maxVal/(norm1*Teuchos::as<Scalar>(nnz)));
    } else {
      std::cout << "ATTENTION: row " << grow << " has only zero entries -> singular matrix!" << std::endl;
    }

  }

  // sort Weights in descending order
  std::vector<int> permutation;
  sortingPermutation(Weights,permutation);

  // create new vector with exactly one possible entry for each column

  // each processor which requests the global column id gcid adds 1 to gColVec
  // gColVec will be summed up over all processors and communicated to gDomVec
  // which is based on the non-overlapping domain map of A.

  Teuchos::RCP<Vector> gColVec = VectorFactory::Build(A->getColMap());
  Teuchos::RCP<Vector> gDomVec = VectorFactory::Build(A->getDomainMap());
  gColVec->putScalar(0.0);
  gDomVec->putScalar(0.0);

  // put in all keep diagonal entries
  for (typename std::vector<std::pair<GlobalOrdinal, GlobalOrdinal> >::const_iterator p = keepDiagonalEntries.begin(); p != keepDiagonalEntries.end(); ++p) {
    gColVec->sumIntoGlobalValue((*p).second,1.0);
  }

  Teuchos::RCP<Export> exporter = ExportFactory::Build(gColVec->getMap(), gDomVec->getMap());
  gDomVec->doExport(*gColVec,*exporter,Xpetra::ADD);  // communicate blocked gcolids to all procs
  gColVec->doImport(*gDomVec,*exporter,Xpetra::INSERT);

  std::vector<std::pair<GlobalOrdinal, GlobalOrdinal> > permutedDiagCandidatesFiltered; // TODO reserve memory
  std::map<GlobalOrdinal, Scalar> gColId2Weight;

  Teuchos::ArrayRCP< Scalar > ddata = gColVec->getDataNonConst(0);
  for(size_t i = 0; i < permutedDiagCandidates.size(); ++i) {
    // loop over all candidates
    std::pair<GlobalOrdinal, GlobalOrdinal> pp = permutedDiagCandidates[permutation[i]];
    GlobalOrdinal grow = pp.first;
    GlobalOrdinal gcol = pp.second;

    LocalOrdinal lcol = A->getColMap()->getLocalElement(gcol);
    //Teuchos::ArrayRCP< Scalar > ddata = gColVec->getDataNonConst(0);
    if(ddata[lcol] > 0.0){
      continue; // skip lcol: column already handled by another row
    }

    // mark column as already taken
    ddata[lcol]++;

    permutedDiagCandidatesFiltered.push_back(std::make_pair(grow,gcol));
    gColId2Weight[gcol] = Weights[permutation[i]];
  }

  // communicate how often each column index is requested by the different procs
  gDomVec->doExport(*gColVec,*exporter,Xpetra::ADD);
  gColVec->doImport(*gDomVec,*exporter,Xpetra::INSERT); // probably not needed // TODO check me

  //*****************************************************************************************
  // first communicate ALL global ids of column indices which are requested by more
  // than one proc to all other procs
  // detect which global col indices are requested by more than one proc
  // and store them in the multipleColRequests vector
  std::vector<GlobalOrdinal> multipleColRequests; // store all global column indices from current processor that are also
                                                  // requested by another processor. This is possible, since they are stored
                                                  // in gDomVec which is based on the nonoverlapping domain map. That is, each
                                                  // global col id is handled by exactly one proc.
  std::queue<GlobalOrdinal> unusedColIdx; // unused column indices on current processor

  for(size_t sz = 0; sz<gDomVec->getLocalLength(); ++sz) {
    Teuchos::ArrayRCP< const Scalar > arrDomVec = gDomVec->getData(0);
    if(arrDomVec[sz] > 1.0) {
      multipleColRequests.push_back(gDomVec->getMap()->getGlobalElement(sz));
    } else if(arrDomVec[sz] == 0.0) {
      unusedColIdx.push(gDomVec->getMap()->getGlobalElement(sz));
    }
  }

  // communicate the global number of column indices which are requested by more than one proc
  LocalOrdinal localMultColRequests  = Teuchos::as<LocalOrdinal>(multipleColRequests.size());
  LocalOrdinal globalMultColRequests = 0;

  // sum up all entries in multipleColRequests over all processors
  sumAll(gDomVec->getMap()->getComm(), (LocalOrdinal)localMultColRequests, globalMultColRequests);

  if(globalMultColRequests > 0) {
    // special handling: two processors request the same global column id.
    // decide which processor gets it

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

    GlobalOrdinal zero=0;
    std::vector<GlobalOrdinal> procMultRequestedColIds(globalMultColRequests,zero);
    std::vector<GlobalOrdinal> global_procMultRequestedColIds(globalMultColRequests,zero);

    // loop over all local column GIDs that are also requested by other procs
    for(size_t i = 0; i < multipleColRequests.size(); i++) {
      procMultRequestedColIds[nMyOffset + i] = multipleColRequests[i]; // all weights are > 0 ?
    }

    // template ordinal, package (double)
    Teuchos::reduceAll<int,GlobalOrdinal>(*comm, Teuchos::REDUCE_MAX, (int) globalMultColRequests,&procMultRequestedColIds[0],&global_procMultRequestedColIds[0]);

    // loop over global_procOverlappingWeights and eliminate wrong entries...
    for (size_t k = 0; k<global_procMultRequestedColIds.size(); k++) {
      GlobalOrdinal globColId = global_procMultRequestedColIds[k];

      std::vector<Scalar> MyWeightForColId(numProcs,0);
      std::vector<Scalar> GlobalWeightForColId(numProcs,0);

      if(gColVec->getMap()->isNodeGlobalElement(globColId)) {
        MyWeightForColId[myRank] = gColId2Weight[globColId];
      } else {
        MyWeightForColId[myRank] = 0.0;
      }

      Teuchos::reduceAll<int,Scalar>(*comm,Teuchos::REDUCE_MAX,numProcs,&MyWeightForColId[0],&GlobalWeightForColId[0]);

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
  //size_t sizeRowColPairs = keepDiagonalEntries.size() + permutedDiagCandidatesFiltered.size();
  std::vector<std::pair<GlobalOrdinal, GlobalOrdinal> > RowColPairs;
  RowColPairs.insert( RowColPairs.end(), keepDiagonalEntries.begin(), keepDiagonalEntries.end());
  RowColPairs.insert( RowColPairs.end(), permutedDiagCandidatesFiltered.begin(), permutedDiagCandidatesFiltered.end());

#ifdef DEBUG_OUTPUT
  //&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
  // plausibility check
  gColVec->putScalar(0.0);
  gDomVec->putScalar(0.0);
  typename std::vector<std::pair<GlobalOrdinal, GlobalOrdinal> >::iterator pl = RowColPairs.begin();
  while(pl != RowColPairs.end() )
  {
    //GlobalOrdinal ik = (*pl).first;
    GlobalOrdinal jk = (*pl).second;

    gColVec->sumIntoGlobalValue(jk,1.0);
    pl++;
  }
  gDomVec->doExport(*gColVec,*exporter,Xpetra::ADD);
  for(size_t sz = 0; sz<gDomVec->getLocalLength(); ++sz) {
    Teuchos::ArrayRCP< const Scalar > arrDomVec = gDomVec->getData(0);
    if(arrDomVec[sz] > 1.0) {
      GetOStream(Runtime0,0) << "RowColPairs has multiple column [" << sz << "]=" << arrDomVec[sz] << std::endl;
    } else if(arrDomVec[sz] == 0.0) {
      GetOStream(Runtime0,0) << "RowColPairs has empty column [" << sz << "]=" << arrDomVec[sz] << std::endl;
    }
  }
  //&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
#endif

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

#ifdef DEBUG_OUTPUT
  //  for (typename std::vector<std::pair<GlobalOrdinal, GlobalOrdinal> >::const_iterator p = RowColPairs.begin(); p != RowColPairs.end(); ++p) {
  //    std::cout << "proc: " << myRank << " r/c: " << (*p).first << "/" << (*p).second << std::endl;
  //  }
  //    for (typename std::vector<std::pair<GlobalOrdinal, GlobalOrdinal> >::const_iterator p = RowColPairs.begin(); p != RowColPairs.end(); ++p)
  //    {
  ////      if((*p).first != (*p).second) std::cout << "difference: " << (*p).first << " " << (*p).second << std::endl;
  //      std::cout << (*p).first +1 << " " << (*p).second+1 << std::endl;
  //    }
  //    std::cout << "\n";
#endif

  // vectors to store permutation information
  Teuchos::RCP<Vector> Pperm  = VectorFactory::Build(A->getRowMap());
  Teuchos::RCP<Vector> Qperm  = VectorFactory::Build(A->getDomainMap()); // global variant (based on domain map)
  Teuchos::RCP<Vector> lQperm = VectorFactory::Build(A->getColMap());  // local variant (based on column map)

  Teuchos::ArrayRCP< Scalar > PpermData = Pperm->getDataNonConst(0);
  Teuchos::ArrayRCP< Scalar > QpermData = Qperm->getDataNonConst(0);

  Pperm->putScalar(0.0);
  Qperm->putScalar(0.0);
  lQperm->putScalar(0.0);

  // setup exporter for Qperm
  Teuchos::RCP<Export> QpermExporter = ExportFactory::Build(lQperm->getMap(), Qperm->getMap());

  Teuchos::RCP<Vector> RowIdStatus  = VectorFactory::Build(A->getRowMap());
  Teuchos::RCP<Vector> ColIdStatus  = VectorFactory::Build(A->getDomainMap()); // global variant (based on domain map)
  Teuchos::RCP<Vector> lColIdStatus = VectorFactory::Build(A->getColMap()); // local variant (based on column map)
  Teuchos::RCP<Vector> ColIdUsed   = VectorFactory::Build(A->getDomainMap()); // mark column ids to be already in use
  Teuchos::ArrayRCP< Scalar > RowIdStatusArray = RowIdStatus->getDataNonConst(0);
  Teuchos::ArrayRCP< Scalar > ColIdStatusArray = ColIdStatus->getDataNonConst(0);
  Teuchos::ArrayRCP< Scalar > lColIdStatusArray = lColIdStatus->getDataNonConst(0);
  Teuchos::ArrayRCP< Scalar > ColIdUsedArray   = ColIdUsed->getDataNonConst(0); // not sure about this
  RowIdStatus->putScalar(0.0);
  ColIdStatus->putScalar(0.0);
  lColIdStatus->putScalar(0.0);
  ColIdUsed->putScalar(0.0);   // no column ids are used

  // count wide-range permutations
  // a wide-range permutation is defined as a permutation of rows/columns which do not
  // belong to the same node
  LocalOrdinal lWideRangeRowPermutations = 0;
  GlobalOrdinal gWideRangeRowPermutations = 0;
  LocalOrdinal lWideRangeColPermutations = 0;
  GlobalOrdinal gWideRangeColPermutations = 0;

  // run 1: mark all "identity" permutations
  typename std::vector<std::pair<GlobalOrdinal, GlobalOrdinal> >::iterator p = RowColPairs.begin();
  while(p != RowColPairs.end() )
  {
    GlobalOrdinal ik = (*p).first;
    GlobalOrdinal jk = (*p).second;

    LocalOrdinal lik = A->getRowMap()->getLocalElement(ik);
    LocalOrdinal ljk = A->getColMap()->getLocalElement(jk);

    if(RowIdStatusArray[lik] == 0.0) {
      RowIdStatusArray[lik] = 1.0; // use this row id
      lColIdStatusArray[ljk] = 1.0; // use this column id
      Pperm->replaceLocalValue(lik, ik);
      lQperm->replaceLocalValue(ljk, ik); // use column map
      ColIdUsed->replaceGlobalValue(ik,1.0); // ik is now used
      p = RowColPairs.erase(p);

      // detect wide range permutations
      if(floor(ik/nDofsPerNode) != floor(jk/nDofsPerNode)) {
        lWideRangeColPermutations++;
      }
    }
    else
      p++;
  }

  // communicate column map -> domain map
  Qperm->doExport(*lQperm,*QpermExporter,Xpetra::ABSMAX);
  ColIdStatus->doExport(*lColIdStatus,*QpermExporter,Xpetra::ABSMAX);

  // plausibility check
  if(RowColPairs.size()>0) GetOStream(Warnings0,0) << "MueLu::PermutationFactory: There are Row/Col pairs left!!!" << std::endl; // TODO fix me

  // close Pperm

  // count, how many row permutations are missing on current proc
  size_t cntFreeRowIdx = 0;
  std::queue<GlobalOrdinal> qFreeGRowIdx;  // store global row ids of "free" rows
  for (size_t lik = 0; lik < RowIdStatus->getLocalLength(); ++lik) {
    if(RowIdStatusArray[lik] == 0.0) {
      cntFreeRowIdx++;
      qFreeGRowIdx.push(RowIdStatus->getMap()->getGlobalElement(lik));
    }
  }

  // fix Pperm
  for (size_t lik = 0; lik < RowIdStatus->getLocalLength(); ++lik) {
    if(RowIdStatusArray[lik] == 0.0) {
      RowIdStatusArray[lik] = 1.0; // use this row id
      Pperm->replaceLocalValue(lik, qFreeGRowIdx.front());
      // detect wide range permutations
      if(floor(qFreeGRowIdx.front()/nDofsPerNode) != floor(RowIdStatus->getMap()->getGlobalElement(lik)/nDofsPerNode)) {
        lWideRangeRowPermutations++;
      }
      qFreeGRowIdx.pop();
    }
  }

  // close Qperm (free permutation entries in Qperm)
  size_t cntFreeColIdx = 0;
  std::queue<GlobalOrdinal> qFreeGColIdx;  // store global column ids of "free" available columns
  for (size_t ljk = 0; ljk < ColIdStatus->getLocalLength(); ++ljk) {
    if(ColIdStatusArray[ljk] == 0.0) {
      cntFreeColIdx++;
      qFreeGColIdx.push(ColIdStatus->getMap()->getGlobalElement(ljk));
    }
  }

  size_t cntUnusedColIdx = 0;
  std::queue<GlobalOrdinal> qUnusedGColIdx;  // store global column ids of "free" available columns
  for (size_t ljk = 0; ljk < ColIdUsed->getLocalLength(); ++ljk) {
    if(ColIdUsedArray[ljk] == 0.0) {
      cntUnusedColIdx++;
      qUnusedGColIdx.push(ColIdUsed->getMap()->getGlobalElement(ljk));
    }
  }

  // fix Qperm with local entries
  for (size_t ljk = 0; ljk < ColIdStatus->getLocalLength(); ++ljk) {
    // stop if no (local) unused column idx are left
    if(cntUnusedColIdx == 0) break;

    if(ColIdStatusArray[ljk] == 0.0) {
      //lColIdStatusArray[ljk] = 1.0; // use this row id
      ColIdStatusArray[ljk] = 1.0; // use this row id
      //ColIdStatus->replaceLocalValue(ljk,1.0); // use this row id
      Qperm->replaceLocalValue(ljk, qUnusedGColIdx.front()); // loop over ColIdStatus (lives on domain map)
      //lQperm->replaceLocalValue(ljk, qUnusedGColIdx.front());
      ColIdUsed->replaceGlobalValue(qUnusedGColIdx.front(),1.0); // ljk is now used, too
      // detect wide range permutations
      if(floor(qUnusedGColIdx.front()/nDofsPerNode) != floor(ColIdStatus->getMap()->getGlobalElement(ljk)/nDofsPerNode)) {
        lWideRangeColPermutations++;
      }
      qUnusedGColIdx.pop();
      cntUnusedColIdx--;
      cntFreeColIdx--;
    }
  }

  //Qperm->doExport(*lQperm,*QpermExporter,Xpetra::ABSMAX); // no export necessary, since changes only locally
  //ColIdStatus->doExport(*lColIdStatus,*QpermExporter,Xpetra::ABSMAX);

  // count, how many unused column idx are needed on current processor
  // to complete Qperm
  cntFreeColIdx = 0;
  for (size_t ljk = 0; ljk < ColIdStatus->getLocalLength(); ++ljk) { // TODO avoid this loop
    if(ColIdStatusArray[ljk] == 0.0) {
      cntFreeColIdx++;
    }
  }

  GlobalOrdinal global_cntFreeColIdx = 0;
  LocalOrdinal  local_cntFreeColIdx = cntFreeColIdx;
  sumAll(comm, (LocalOrdinal)local_cntFreeColIdx, global_cntFreeColIdx);
#ifdef DEBUG_OUTPUT
  std::cout << "global # of empty column idx entries in Qperm: " << global_cntFreeColIdx << std::endl;
#endif

  // avoid global communication if possible
  if(global_cntFreeColIdx > 0) {

    // 1) count how many unused column ids are left
    GlobalOrdinal global_cntUnusedColIdx = 0;
    LocalOrdinal  local_cntUnusedColIdx = cntUnusedColIdx;
    sumAll(comm, (LocalOrdinal)local_cntUnusedColIdx, global_cntUnusedColIdx);
#ifdef DEBUG_OUTPUT
    std::cout << "global # of unused column idx: " << global_cntUnusedColIdx << std::endl;
#endif

    // 2) communicate how many unused column ids are available on procs
    std::vector<LocalOrdinal> local_UnusedColIdxOnProc (numProcs);
    std::vector<LocalOrdinal> global_UnusedColIdxOnProc(numProcs);
    local_UnusedColIdxOnProc[myRank] = local_cntUnusedColIdx;
    Teuchos::reduceAll<int,LocalOrdinal>(*comm,Teuchos::REDUCE_MAX,numProcs,&local_UnusedColIdxOnProc[0],&global_UnusedColIdxOnProc[0]);

#ifdef DEBUG_OUTPUT
    std::cout << "PROC " << myRank << " global num unused indices per proc: ";
    for (size_t ljk = 0; ljk < global_UnusedColIdxOnProc.size(); ++ljk) {
      std::cout << " " << global_UnusedColIdxOnProc[ljk];
    }
    std::cout << std::endl;
#endif

    // 3) build array of length global_cntUnusedColIdx to globally replicate unused column idx
    std::vector<GlobalOrdinal> local_UnusedColIdxVector(Teuchos::as<size_t>(global_cntUnusedColIdx));
    std::vector<GlobalOrdinal> global_UnusedColIdxVector(Teuchos::as<size_t>(global_cntUnusedColIdx));
    GlobalOrdinal global_cntUnusedColIdxStartIter = 0;
    for(int proc=0; proc<myRank; proc++) {
      global_cntUnusedColIdxStartIter += global_UnusedColIdxOnProc[proc];
    }
    for(GlobalOrdinal k = global_cntUnusedColIdxStartIter; k < global_cntUnusedColIdxStartIter+local_cntUnusedColIdx; k++) {
      local_UnusedColIdxVector[k] = qUnusedGColIdx.front();
      qUnusedGColIdx.pop();
    }
    Teuchos::reduceAll<int,GlobalOrdinal>(*comm,Teuchos::REDUCE_MAX,(int)global_cntUnusedColIdx,&local_UnusedColIdxVector[0],&global_UnusedColIdxVector[0]);
#ifdef DEBUG_OUTPUT
    std::cout << "PROC " << myRank << " global UnusedGColIdx: ";
    for (size_t ljk = 0; ljk < global_UnusedColIdxVector.size(); ++ljk) {
      std::cout << " " << global_UnusedColIdxVector[ljk];
    }
    std::cout << std::endl;
#endif



    // 4) communicate, how many column idx are needed on each processor
    //    to complete Qperm
    std::vector<LocalOrdinal> local_EmptyColIdxOnProc (numProcs);
    std::vector<LocalOrdinal> global_EmptyColIdxOnProc(numProcs);
    local_EmptyColIdxOnProc[myRank] = local_cntFreeColIdx;
    Teuchos::reduceAll<int,LocalOrdinal>(*comm,Teuchos::REDUCE_MAX,numProcs,&local_EmptyColIdxOnProc[0],&global_EmptyColIdxOnProc[0]);

#ifdef DEBUG_OUTPUT
    std::cout << "PROC " << myRank << " global num of needed column indices: ";
    for (size_t ljk = 0; ljk < global_EmptyColIdxOnProc.size(); ++ljk) {
      std::cout << " " << global_EmptyColIdxOnProc[ljk];
    }
    std::cout << std::endl;
#endif

    // 5) determine first index in global_UnusedColIdxVector for unused column indices,
    //    that are marked to be used by this processor
    GlobalOrdinal global_UnusedColStartIdx = 0;
    for(int proc=0; proc<myRank; proc++) {
      global_UnusedColStartIdx += global_EmptyColIdxOnProc[proc];
    }

#ifdef DEBUG_OUTPUT
    GetOStream(Statistics0,0) << "PROC " << myRank << " is allowd to use the following column gids: ";
    for(GlobalOrdinal k = global_UnusedColStartIdx; k < global_UnusedColStartIdx + Teuchos::as<GlobalOrdinal>(cntFreeColIdx); k++) {
      GetOStream(Statistics0,0) << global_UnusedColIdxVector[k] << " ";
    }
    GetOStream(Statistics0,0) << std::endl;
#endif

    // 6.) fix Qperm with global entries
    GlobalOrdinal array_iter = 0;
    for (size_t ljk = 0; ljk < ColIdStatus->getLocalLength(); ++ljk) {

      if(ColIdStatusArray[ljk] == 0.0) {
        ColIdStatusArray[ljk] = 1.0; // use this row id
        Qperm->replaceLocalValue(ljk, global_UnusedColIdxVector[global_UnusedColStartIdx + array_iter]);
        //lQperm->replaceLocalValue(ljk, global_UnusedColIdxVector[global_UnusedColStartIdx + array_iter]);
        ColIdUsed->replaceGlobalValue(global_UnusedColIdxVector[global_UnusedColStartIdx + array_iter],1.0);
        // detect wide range permutations
        if(floor(global_UnusedColIdxVector[global_UnusedColStartIdx + array_iter]/nDofsPerNode) != floor(ColIdStatus->getMap()->getGlobalElement(ljk)/nDofsPerNode)) {
          lWideRangeColPermutations++;
        }
        array_iter++;
        //cntUnusedColIdx--; // check me
      }
    }
    //Qperm->doExport(*lQperm,*QpermExporter,Xpetra::ABSMAX);
    //ColIdStatus->doExport(*lColIdStatus,*QpermExporter,Xpetra::ABSMAX);
  } // end if global_cntFreeColIdx > 0
  /////////////////// Qperm should be fine now...


  // create new empty Matrix
  Teuchos::RCP<CrsMatrixWrap> permPTmatrix = Teuchos::rcp(new CrsMatrixWrap(A->getRowMap(),1,Xpetra::StaticProfile));
  Teuchos::RCP<CrsMatrixWrap> permQTmatrix = Teuchos::rcp(new CrsMatrixWrap(A->getRowMap(),1,Xpetra::StaticProfile));

  for(size_t row=0; row<A->getNodeNumRows(); row++) {
    Teuchos::ArrayRCP<GlobalOrdinal> indoutP(1,Teuchos::as<GO>(PpermData[row])); // column idx for Perm^T
    Teuchos::ArrayRCP<GlobalOrdinal> indoutQ(1,Teuchos::as<GO>(QpermData[row])); // column idx for Qperm
    Teuchos::ArrayRCP<Scalar> valout(1,1.0);
    permPTmatrix->insertGlobalValues(A->getRowMap()->getGlobalElement(row), indoutP.view(0,indoutP.size()), valout.view(0,valout.size()));
    permQTmatrix->insertGlobalValues (A->getRowMap()->getGlobalElement(row), indoutQ.view(0,indoutQ.size()), valout.view(0,valout.size()));
  }

  permPTmatrix->fillComplete();
  permQTmatrix->fillComplete();

  Teuchos::RCP<Matrix> permPmatrix = Utils2::Transpose(permPTmatrix,true);

  for(size_t row=0; row<permPTmatrix->getNodeNumRows(); row++) {
    if(permPTmatrix->getNumEntriesInLocalRow(row) != 1)
      GetOStream(Warnings0,0) <<"#entries in row " << row << " of permPTmatrix is " << permPTmatrix->getNumEntriesInLocalRow(row) << std::endl;
    if(permPmatrix->getNumEntriesInLocalRow(row) != 1)
      GetOStream(Warnings0,0) <<"#entries in row " << row << " of permPmatrix is " << permPmatrix->getNumEntriesInLocalRow(row) << std::endl;
    if(permQTmatrix->getNumEntriesInLocalRow(row) != 1)
      GetOStream(Warnings0,0) <<"#entries in row " << row << " of permQmatrix is " << permQTmatrix->getNumEntriesInLocalRow(row) << std::endl;
  }

  // build permP * A * permQT
  Teuchos::RCP<Matrix> ApermQt = Utils::TwoMatrixMultiply(A, false, permQTmatrix, false);
  Teuchos::RCP<Matrix> permPApermQt = Utils::TwoMatrixMultiply(Teuchos::rcp_dynamic_cast<Matrix>(permPmatrix), false, ApermQt, false);

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
      GetOStream(Statistics0,0) << "MueLu::PermutationFactory: found zero on diagonal in row " << i << std::endl;
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

  //// count row permutations
  // count zeros on diagonal in P -> number of row permutations
  Teuchos::RCP<Vector> diagPVec = VectorFactory::Build(permPmatrix->getRowMap(),true);
  permPmatrix->getLocalDiagCopy(*diagPVec);
  Teuchos::ArrayRCP< const Scalar > diagPVecData = diagPVec->getData(0);
  LocalOrdinal lNumRowPermutations = 0;
  GlobalOrdinal gNumRowPermutations = 0;
  for(size_t i = 0; i<diagPVec->getMap()->getNodeNumElements(); ++i) {
    if(diagPVecData[i] == 0.0) {
      lNumRowPermutations++;
    }
  }

  // sum up all entries in multipleColRequests over all processors
  sumAll(diagPVec->getMap()->getComm(), (LocalOrdinal)lNumRowPermutations, gNumRowPermutations);

  //// count column permutations
  // count zeros on diagonal in Q^T -> number of column permutations
  Teuchos::RCP<Vector> diagQTVec = VectorFactory::Build(permQTmatrix->getRowMap(),true);
  permQTmatrix->getLocalDiagCopy(*diagQTVec);
  Teuchos::ArrayRCP< const Scalar > diagQTVecData = diagQTVec->getData(0);
  LocalOrdinal lNumColPermutations = 0;
  GlobalOrdinal gNumColPermutations = 0;
  for(size_t i = 0; i<diagQTVec->getMap()->getNodeNumElements(); ++i) {
    if(diagQTVecData[i] == 0.0) {
      lNumColPermutations++;
    }
  }

  // sum up all entries in multipleColRequests over all processors
  sumAll(diagQTVec->getMap()->getComm(), (LocalOrdinal)lNumColPermutations, gNumColPermutations);

  currentLevel.Set("#RowPermutations", gNumRowPermutations, this);
  currentLevel.Set("#ColPermutations", gNumColPermutations, this);
  currentLevel.Set("#WideRangeRowPermutations", gWideRangeRowPermutations, this);
  currentLevel.Set("#WideRangeColPermutations", gWideRangeColPermutations, this);

  GetOStream(Statistics0, 0) << "#Row    permutations/max possible permutations: " << gNumRowPermutations << "/" << diagPVec->getMap()->getGlobalNumElements() << std::endl;
  GetOStream(Statistics0, 0) << "#Column permutations/max possible permutations: " << gNumColPermutations << "/" << diagQTVec->getMap()->getGlobalNumElements() << std::endl;
  GetOStream(Runtime1, 0) << "#wide range row permutations: " << gWideRangeRowPermutations << " #wide range column permutations: " << gWideRangeColPermutations << std::endl;

#else
#warning PermutationFactory not compiling/working for Scalar==complex.
#endif // #ifndef HAVE_MUELU_INST_COMPLEX_INT_INT


}

} // namespace MueLu


#endif /* MUELU_PERMUTATIONFACTORY_DEF_HPP_ */
