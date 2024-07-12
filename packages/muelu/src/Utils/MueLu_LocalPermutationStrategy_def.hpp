// @HEADER
// *****************************************************************************
//        MueLu: A package for multigrid based preconditioning
//
// Copyright 2012 NTESS and the MueLu contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/*
 * MueLu_LocalPermutationStrategy_def.hpp
 *
 *  Created on: Feb 19, 2013
 *      Author: tobias
 */

#ifndef MUELU_LOCALPERMUTATIONSTRATEGY_DEF_HPP_
#define MUELU_LOCALPERMUTATIONSTRATEGY_DEF_HPP_

#include <algorithm>

#include <Xpetra_MultiVector.hpp>
#include <Xpetra_Matrix.hpp>
#include <Xpetra_MatrixMatrix.hpp>
#include <Xpetra_CrsGraph.hpp>
#include <Xpetra_Vector.hpp>
#include <Xpetra_VectorFactory.hpp>
#include <Xpetra_CrsMatrixWrap.hpp>

#include "MueLu_Utilities.hpp"
#include "MueLu_LocalPermutationStrategy_decl.hpp"

namespace MueLu {

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void LocalPermutationStrategy<Scalar, LocalOrdinal, GlobalOrdinal, Node>::BuildPermutations(size_t nDofsPerNode) const {
  permWidth_ = nDofsPerNode;

  result_permvecs_.clear();

  // build permutation string
  std::stringstream ss;
  for (size_t t = 0; t < nDofsPerNode; t++)
    ss << t;
  std::string cs = ss.str();
  // std::vector<std::string> result_perms;
  do {
    // result_perms.push_back(cs);

    std::vector<int> newPerm(cs.length(), -1);
    for (size_t len = 0; len < cs.length(); len++) {
      newPerm[len] = Teuchos::as<int>(cs[len] - '0');
    }
    result_permvecs_.push_back(newPerm);

  } while (std::next_permutation(cs.begin(), cs.end()));
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void LocalPermutationStrategy<Scalar, LocalOrdinal, GlobalOrdinal, Node>::BuildPermutation(const Teuchos::RCP<Matrix>& A, const Teuchos::RCP<const Map> /* permRowMap */, Level& currentLevel, const FactoryBase* genFactory) const {
  SC SC_ZERO = Teuchos::ScalarTraits<SC>::zero();

  size_t nDofsPerNode = 1;
  if (A->IsView("stridedMaps")) {
    Teuchos::RCP<const Map> permRowMapStrided = A->getRowMap("stridedMaps");
    nDofsPerNode                              = Teuchos::rcp_dynamic_cast<const StridedMap>(permRowMapStrided)->getFixedBlockSize();
  }

  RCP<const Teuchos::Comm<int> > comm = A->getRowMap()->getComm();

  //////////////////
  std::vector<std::pair<GlobalOrdinal, GlobalOrdinal> > RowColPairs;

  // check whether we have to (re)build the permutation vector
  if (permWidth_ != nDofsPerNode)
    BuildPermutations(nDofsPerNode);

  // declare local variables used inside loop
  LocalOrdinal lonDofsPerNode = Teuchos::as<LocalOrdinal>(nDofsPerNode);
  Teuchos::ArrayView<const LocalOrdinal> indices;
  Teuchos::ArrayView<const Scalar> vals;
  Teuchos::SerialDenseMatrix<LocalOrdinal, Scalar> subBlockMatrix(nDofsPerNode, nDofsPerNode, true);
  std::vector<GlobalOrdinal> growIds(nDofsPerNode);

  // loop over local nodes
  // TODO what about nOffset?
  LocalOrdinal numLocalNodes = A->getRowMap()->getLocalNumElements() / nDofsPerNode;
  for (LocalOrdinal node = 0; node < numLocalNodes; ++node) {
    // zero out block matrix
    subBlockMatrix.putScalar();

    // loop over all DOFs in current node
    // Note: were assuming constant number of Dofs per node here!
    // TODO This is more complicated for variable dofs per node
    for (LocalOrdinal lrdof = 0; lrdof < lonDofsPerNode; ++lrdof) {
      GlobalOrdinal grow = getGlobalDofId(A, node, lrdof);
      growIds[lrdof]     = grow;

      // extract local row information from matrix
      A->getLocalRowView(A->getRowMap()->getLocalElement(grow), indices, vals);

      // find column entry with max absolute value
      typedef typename Teuchos::ScalarTraits<Scalar>::magnitudeType MT;
      MT maxVal = 0.0;
      for (size_t j = 0; j < Teuchos::as<size_t>(indices.size()); j++) {
        if (Teuchos::ScalarTraits<Scalar>::magnitude(vals[j]) > maxVal) {
          maxVal = Teuchos::ScalarTraits<Scalar>::magnitude(vals[j]);
        }
      }

      GlobalOrdinal grnodeid = globalDofId2globalNodeId(A, grow);

      for (size_t j = 0; j < Teuchos::as<size_t>(indices.size()); j++) {
        GlobalOrdinal gcol     = A->getColMap()->getGlobalElement(indices[j]);
        GlobalOrdinal gcnodeid = globalDofId2globalNodeId(A, gcol);  // -> global node id
        if (grnodeid == gcnodeid) {
          if (maxVal != Teuchos::ScalarTraits<MT>::zero()) {
            subBlockMatrix(lrdof, gcol % nDofsPerNode) = vals[j] / maxVal;
          } else {
            subBlockMatrix(lrdof, gcol % nDofsPerNode) = vals[j];  // there is a problem
            std::cout << "maxVal never should be zero!!!!" << std::endl;
          }
        }
      }
    }

    // now we have the sub block matrix

    // build permutation string
    /*std::stringstream ss;
    for(size_t t = 0; t<nDofsPerNode; t++)
      ss << t;
    std::string cs = ss.str();
    std::vector<std::string> result_perms;
    do {
      result_perms.push_back(cs);
      //std::cout << result_perms.back() << std::endl;
    } while (std::next_permutation(cs.begin(),cs.end()));*/

    std::vector<Scalar> performance_vector = std::vector<Scalar>(result_permvecs_.size());
    for (size_t t = 0; t < result_permvecs_.size(); t++) {
      std::vector<int> vv = result_permvecs_[t];
      Scalar value        = 1.0;
      for (size_t j = 0; j < vv.size(); j++) {
        value = value * subBlockMatrix(j, vv[j]);
      }
      performance_vector[t] = value;
    }
    /*for(size_t t = 0; t < result_perms.size(); t++) {
      std::string s = result_perms[t];
      Scalar value = 1.0;
      for(size_t len=0; len<s.length(); len++) {
        int col = Teuchos::as<int>(s[len]-'0');
        value = value * subBlockMatrix(len,col);
      }
      performance_vector[t] = value;
    }*/

    // find permutation with maximum performance value
    typedef typename Teuchos::ScalarTraits<Scalar>::magnitudeType MT;
    MT maxVal                           = Teuchos::ScalarTraits<MT>::zero();
    size_t maxPerformancePermutationIdx = 0;
    for (size_t j = 0; j < Teuchos::as<size_t>(performance_vector.size()); j++) {
      if (Teuchos::ScalarTraits<Scalar>::magnitude(performance_vector[j]) > maxVal) {
        maxVal                       = Teuchos::ScalarTraits<Scalar>::magnitude(performance_vector[j]);
        maxPerformancePermutationIdx = j;
      }
    }

    // build RowColPairs for best permutation
    std::vector<int> bestPerformancePermutation = result_permvecs_[maxPerformancePermutationIdx];
    for (size_t t = 0; t < nDofsPerNode; t++) {
      RowColPairs.push_back(std::make_pair(growIds[t], growIds[bestPerformancePermutation[t]]));
    }

  }  // end loop over local nodes

  // build Pperm and Qperm vectors
  Teuchos::RCP<Vector> Pperm = VectorFactory::Build(A->getRowMap());
  Teuchos::RCP<Vector> Qperm = VectorFactory::Build(A->getDomainMap());

  Pperm->putScalar(SC_ZERO);
  Qperm->putScalar(SC_ZERO);

  Teuchos::ArrayRCP<Scalar> PpermData = Pperm->getDataNonConst(0);
  Teuchos::ArrayRCP<Scalar> QpermData = Qperm->getDataNonConst(0);

  typename std::vector<std::pair<GlobalOrdinal, GlobalOrdinal> >::iterator p = RowColPairs.begin();
  while (p != RowColPairs.end()) {
    GlobalOrdinal ik = (*p).first;
    GlobalOrdinal jk = (*p).second;

    LocalOrdinal lik = A->getRowMap()->getLocalElement(ik);
    LocalOrdinal ljk = A->getDomainMap()->getLocalElement(jk);

    Pperm->replaceLocalValue(lik, ik);
    Qperm->replaceLocalValue(ljk, ik);

    p = RowColPairs.erase(p);
  }

  if (RowColPairs.size() > 0) GetOStream(Warnings0) << "MueLu::LocalPermutationStrategy: There are Row/col pairs left!" << std::endl;

  // Qperm should be fine
  // build matrices

  // create new empty Matrix
  Teuchos::RCP<CrsMatrixWrap> permPTmatrix = Teuchos::rcp(new CrsMatrixWrap(A->getRowMap(), 1));
  Teuchos::RCP<CrsMatrixWrap> permQTmatrix = Teuchos::rcp(new CrsMatrixWrap(A->getRowMap(), 1));

  for (size_t row = 0; row < A->getLocalNumRows(); row++) {
    Teuchos::ArrayRCP<GlobalOrdinal> indoutP(1, Teuchos::as<GO>(Teuchos::ScalarTraits<Scalar>::real(PpermData[row])));  // column idx for Perm^T
    Teuchos::ArrayRCP<GlobalOrdinal> indoutQ(1, Teuchos::as<GO>(Teuchos::ScalarTraits<Scalar>::real(QpermData[row])));  // column idx for Qperm
    Teuchos::ArrayRCP<Scalar> valout(1, Teuchos::ScalarTraits<Scalar>::one());
    permPTmatrix->insertGlobalValues(A->getRowMap()->getGlobalElement(row), indoutP.view(0, indoutP.size()), valout.view(0, valout.size()));
    permQTmatrix->insertGlobalValues(A->getRowMap()->getGlobalElement(row), indoutQ.view(0, indoutQ.size()), valout.view(0, valout.size()));
  }

  permPTmatrix->fillComplete();
  permQTmatrix->fillComplete();

  Teuchos::RCP<Matrix> permPmatrix = Utilities::Transpose(*permPTmatrix, true);

  /*for(size_t row=0; row<permPTmatrix->getLocalNumRows(); row++) {
    if(permPTmatrix->getNumEntriesInLocalRow(row) != 1)
      GetOStream(Warnings0) <<"#entries in row " << row << " of permPTmatrix is " << permPTmatrix->getNumEntriesInLocalRow(row) << std::endl;
    if(permPmatrix->getNumEntriesInLocalRow(row) != 1)
      GetOStream(Warnings0) <<"#entries in row " << row << " of permPmatrix is " << permPmatrix->getNumEntriesInLocalRow(row) << std::endl;
    if(permQTmatrix->getNumEntriesInLocalRow(row) != 1)
      GetOStream(Warnings0) <<"#entries in row " << row << " of permQmatrix is " << permQTmatrix->getNumEntriesInLocalRow(row) << std::endl;
  }*/

  // build permP * A * permQT
  Teuchos::RCP<Matrix> ApermQt      = Xpetra::MatrixMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Multiply(*A, false, *permQTmatrix, false, GetOStream(Statistics2), true, true);
  Teuchos::RCP<Matrix> permPApermQt = Xpetra::MatrixMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Multiply(*permPmatrix, false, *ApermQt, false, GetOStream(Statistics2), true, true);

  /*
  MueLu::Utils<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Write("A.mat", *A);
  MueLu::Utils<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Write("permP.mat", *permPmatrix);
  MueLu::Utils<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Write("permQt.mat", *permQTmatrix);
  MueLu::Utils<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Write("permPApermQt.mat", *permPApermQt);
  */
  // build scaling matrix
  Teuchos::RCP<Vector> diagVec             = VectorFactory::Build(permPApermQt->getRowMap(), true);
  Teuchos::RCP<Vector> invDiagVec          = VectorFactory::Build(permPApermQt->getRowMap(), true);
  Teuchos::ArrayRCP<Scalar> invDiagVecData = invDiagVec->getDataNonConst(0);

  LO lCntZeroDiagonals = 0;
  permPApermQt->getLocalDiagCopy(*diagVec);
  Teuchos::ArrayRCP<const Scalar> diagVecData = diagVec->getData(0);
  for (size_t i = 0; i < diagVec->getMap()->getLocalNumElements(); ++i) {
    if (diagVecData[i] != SC_ZERO)
      invDiagVecData[i] = Teuchos::ScalarTraits<Scalar>::one() / diagVecData[i];
    else {
      invDiagVecData[i] = Teuchos::ScalarTraits<Scalar>::one();
      lCntZeroDiagonals++;
      // GetOStream(Statistics0) << "MueLu::LocalPermutationStrategy: found zero on diagonal in row " << i << std::endl;
    }
  }

#if 0
    GO gCntZeroDiagonals  = 0;
    GO glCntZeroDiagonals = Teuchos::as<GlobalOrdinal>(lCntZeroDiagonals);  /* LO->GO conversion */
    MueLu_sumAll(comm,glCntZeroDiagonals,gCntZeroDiagonals);
    GetOStream(Statistics0) << "MueLu::LocalPermutationStrategy: found " << gCntZeroDiagonals << " zeros on diagonal" << std::endl;
#endif

  Teuchos::RCP<CrsMatrixWrap> diagScalingOp = Teuchos::rcp(new CrsMatrixWrap(permPApermQt->getRowMap(), 1));

  for (size_t row = 0; row < A->getLocalNumRows(); row++) {
    Teuchos::ArrayRCP<GlobalOrdinal> indout(1, permPApermQt->getRowMap()->getGlobalElement(row));  // column idx for Perm^T
    Teuchos::ArrayRCP<Scalar> valout(1, invDiagVecData[row]);
    diagScalingOp->insertGlobalValues(A->getRowMap()->getGlobalElement(row), indout.view(0, indout.size()), valout.view(0, valout.size()));
  }
  diagScalingOp->fillComplete();

  Teuchos::RCP<Matrix> scaledA = Xpetra::MatrixMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Multiply(*diagScalingOp, false, *permPApermQt, false, GetOStream(Statistics2), true, true);
  currentLevel.Set("A", Teuchos::rcp_dynamic_cast<Matrix>(scaledA), genFactory);

  currentLevel.Set("permA", Teuchos::rcp_dynamic_cast<Matrix>(permPApermQt), genFactory);
  currentLevel.Set("permP", Teuchos::rcp_dynamic_cast<Matrix>(permPmatrix), genFactory);
  currentLevel.Set("permQT", Teuchos::rcp_dynamic_cast<Matrix>(permQTmatrix), genFactory);
  currentLevel.Set("permScaling", Teuchos::rcp_dynamic_cast<Matrix>(diagScalingOp), genFactory);

  //// count row permutations
  // count zeros on diagonal in P -> number of row permutations
  Teuchos::RCP<Vector> diagPVec = VectorFactory::Build(permPmatrix->getRowMap(), true);
  permPmatrix->getLocalDiagCopy(*diagPVec);
  Teuchos::ArrayRCP<const Scalar> diagPVecData = diagPVec->getData(0);
  GlobalOrdinal lNumRowPermutations            = 0;
  GlobalOrdinal gNumRowPermutations            = 0;
  for (size_t i = 0; i < diagPVec->getMap()->getLocalNumElements(); ++i) {
    if (diagPVecData[i] == SC_ZERO) {
      lNumRowPermutations++;
    }
  }

  // sum up all entries in multipleColRequests over all processors
  MueLu_sumAll(diagPVec->getMap()->getComm(), lNumRowPermutations, gNumRowPermutations);

  //// count column permutations
  // count zeros on diagonal in Q^T -> number of column permutations
  Teuchos::RCP<Vector> diagQTVec = VectorFactory::Build(permQTmatrix->getRowMap(), true);
  permQTmatrix->getLocalDiagCopy(*diagQTVec);
  Teuchos::ArrayRCP<const Scalar> diagQTVecData = diagQTVec->getData(0);
  GlobalOrdinal lNumColPermutations             = 0;
  GlobalOrdinal gNumColPermutations             = 0;
  for (size_t i = 0; i < diagQTVec->getMap()->getLocalNumElements(); ++i) {
    if (diagQTVecData[i] == SC_ZERO) {
      lNumColPermutations++;
    }
  }

  // sum up all entries in multipleColRequests over all processors
  MueLu_sumAll(diagQTVec->getMap()->getComm(), lNumColPermutations, gNumColPermutations);

  currentLevel.Set("#RowPermutations", gNumRowPermutations, genFactory /*this*/);
  currentLevel.Set("#ColPermutations", gNumColPermutations, genFactory /*this*/);
  currentLevel.Set("#WideRangeRowPermutations", 0, genFactory /*this*/);
  currentLevel.Set("#WideRangeColPermutations", 0, genFactory /*this*/);

  GetOStream(Statistics0) << "#Row    permutations/max possible permutations: " << gNumRowPermutations << "/" << diagPVec->getMap()->getGlobalNumElements() << std::endl;
  GetOStream(Statistics0) << "#Column permutations/max possible permutations: " << gNumColPermutations << "/" << diagQTVec->getMap()->getGlobalNumElements() << std::endl;
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
GlobalOrdinal LocalPermutationStrategy<Scalar, LocalOrdinal, GlobalOrdinal, Node>::getGlobalDofId(const Teuchos::RCP<Matrix>& A, LocalOrdinal localNodeId, LocalOrdinal localDof) const {
  size_t nDofsPerNode = 1;
  if (A->IsView("stridedMaps")) {
    Teuchos::RCP<const Map> permRowMapStrided = A->getRowMap("stridedMaps");
    nDofsPerNode                              = Teuchos::rcp_dynamic_cast<const StridedMap>(permRowMapStrided)->getFixedBlockSize();
  }

  LocalOrdinal localDofId = localNodeId * nDofsPerNode + localDof;

  return A->getRowMap()->getGlobalElement(localDofId);
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
GlobalOrdinal LocalPermutationStrategy<Scalar, LocalOrdinal, GlobalOrdinal, Node>::globalDofId2globalNodeId(const Teuchos::RCP<Matrix>& A, GlobalOrdinal grid) const {
  size_t nDofsPerNode = 1;
  if (A->IsView("stridedMaps")) {
    Teuchos::RCP<const Map> permRowMapStrided = A->getRowMap("stridedMaps");
    nDofsPerNode                              = Teuchos::rcp_dynamic_cast<const StridedMap>(permRowMapStrided)->getFixedBlockSize();
  }

  return (GlobalOrdinal)grid / (GlobalOrdinal)nDofsPerNode;  // TODO what about nOffset???
}

}  // namespace MueLu

#endif /* MUELU_LOCALPERMUTATIONSTRATEGY_DEF_HPP_ */
