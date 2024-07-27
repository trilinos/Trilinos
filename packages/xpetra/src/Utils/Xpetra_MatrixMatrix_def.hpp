// @HEADER
// *****************************************************************************
//             Xpetra: A linear algebra interface package
//
// Copyright 2012 NTESS and the Xpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef PACKAGES_XPETRA_SUP_UTILS_XPETRA_MATRIXMATRIX_DEF_HPP_
#define PACKAGES_XPETRA_SUP_UTILS_XPETRA_MATRIXMATRIX_DEF_HPP_

#include "Xpetra_ConfigDefs.hpp"

#include "Xpetra_BlockedCrsMatrix.hpp"
#include "Xpetra_CrsMatrixWrap.hpp"
#include "Xpetra_MapExtractor.hpp"
#include "Xpetra_Map.hpp"
#include "Xpetra_MatrixFactory.hpp"
#include "Xpetra_Matrix.hpp"
#include "Xpetra_StridedMapFactory.hpp"
#include "Xpetra_StridedMap.hpp"

#ifdef HAVE_XPETRA_EPETRA
#include <Xpetra_EpetraCrsMatrix_fwd.hpp>
#endif

#ifdef HAVE_XPETRA_EPETRAEXT
#include <EpetraExt_MatrixMatrix.h>
#include <EpetraExt_RowMatrixOut.h>
#include <Epetra_RowMatrixTransposer.h>
#endif  // HAVE_XPETRA_EPETRAEXT

#ifdef HAVE_XPETRA_TPETRA
#include <TpetraExt_MatrixMatrix.hpp>
#include <Tpetra_RowMatrixTransposer.hpp>
#include <MatrixMarket_Tpetra.hpp>
#include <Xpetra_TpetraCrsMatrix.hpp>
#include <Xpetra_TpetraBlockCrsMatrix.hpp>
#include <Tpetra_BlockCrsMatrix_Helpers.hpp>
#include <Xpetra_TpetraMultiVector.hpp>
#include <Xpetra_TpetraVector.hpp>
#endif  // HAVE_XPETRA_TPETRA

#include "Xpetra_MatrixMatrix_decl.hpp"

namespace Xpetra {

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void MatrixMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Multiply(const Matrix& A, bool transposeA,
                                                                       const Matrix& B, bool transposeB,
                                                                       Matrix& C,
                                                                       bool call_FillComplete_on_result,
                                                                       bool doOptimizeStorage,
                                                                       const std::string& label,
                                                                       const RCP<ParameterList>& params) {
  TEUCHOS_TEST_FOR_EXCEPTION(transposeA == false && C.getRowMap()->isSameAs(*A.getRowMap()) == false,
                             Exceptions::RuntimeError, "XpetraExt::MatrixMatrix::Multiply: row map of C is not same as row map of A");
  TEUCHOS_TEST_FOR_EXCEPTION(transposeA == true && C.getRowMap()->isSameAs(*A.getDomainMap()) == false,
                             Exceptions::RuntimeError, "XpetraExt::MatrixMatrix::Multiply: row map of C is not same as domain map of A");

  TEUCHOS_TEST_FOR_EXCEPTION(!A.isFillComplete(), Exceptions::RuntimeError, "A is not fill-completed");
  TEUCHOS_TEST_FOR_EXCEPTION(!B.isFillComplete(), Exceptions::RuntimeError, "B is not fill-completed");

  bool haveMultiplyDoFillComplete = call_FillComplete_on_result && doOptimizeStorage;

  if (C.getRowMap()->lib() == Xpetra::UseEpetra) {
#if defined(HAVE_XPETRA_EPETRA) && defined(HAVE_XPETRA_EPETRAEXT)
    throw(Xpetra::Exceptions::RuntimeError("Xpetra::MatrixMatrix::Multiply only available for GO=int or GO=long long with EpetraNode (Serial or OpenMP depending on configuration)"));
#else
    throw(Xpetra::Exceptions::RuntimeError("Xpetra::MatrixMatrix::Multiply requires EpetraExt to be compiled."));

#endif
  } else if (C.getRowMap()->lib() == Xpetra::UseTpetra) {
#ifdef HAVE_XPETRA_TPETRA
    using helpers = Xpetra::Helpers<SC, LO, GO, NO>;
    if (helpers::isTpetraCrs(A) && helpers::isTpetraCrs(B) && helpers::isTpetraCrs(C)) {
      // All matrices are Crs
      const Tpetra::CrsMatrix<SC, LO, GO, NO>& tpA = helpers::Op2TpetraCrs(A);
      const Tpetra::CrsMatrix<SC, LO, GO, NO>& tpB = helpers::Op2TpetraCrs(B);
      Tpetra::CrsMatrix<SC, LO, GO, NO>& tpC       = helpers::Op2NonConstTpetraCrs(C);

      // 18Feb2013 JJH I'm reenabling the code that allows the matrix matrix multiply to do the fillComplete.
      // Previously, Tpetra's matrix matrix multiply did not support fillComplete.
      Tpetra::MatrixMatrix::Multiply(tpA, transposeA, tpB, transposeB, tpC, haveMultiplyDoFillComplete, label, params);
    } else if (helpers::isTpetraBlockCrs(A) && helpers::isTpetraBlockCrs(B)) {
      // All matrices are BlockCrs (except maybe Ac)
      // FIXME: For the moment we're just going to clobber the innards of Ac, so no reuse. Once we have a reuse kernel,
      // we'll need to think about refactoring BlockCrs so we can do something smarter here.
      if (!A.getRowMap()->getComm()->getRank())
        std::cout << "WARNING: Using inefficient BlockCrs Multiply Placeholder" << std::endl;

      const Tpetra::BlockCrsMatrix<SC, LO, GO, NO>& tpA = Xpetra::Helpers<SC, LO, GO, NO>::Op2TpetraBlockCrs(A);
      const Tpetra::BlockCrsMatrix<SC, LO, GO, NO>& tpB = Xpetra::Helpers<SC, LO, GO, NO>::Op2TpetraBlockCrs(B);
      using CRS                                         = Tpetra::CrsMatrix<SC, LO, GO, NO>;
      RCP<const CRS> Acrs                               = Tpetra::convertToCrsMatrix(tpA);
      RCP<const CRS> Bcrs                               = Tpetra::convertToCrsMatrix(tpB);

      // We need the global constants to do the copy back to BlockCrs
      RCP<ParameterList> new_params;
      if (!params.is_null()) {
        new_params = rcp(new Teuchos::ParameterList(*params));
        new_params->set("compute global constants", true);
      }

      // FIXME: The lines below only works because we're assuming Ac is Point
      RCP<CRS> tempAc = Teuchos::rcp(new CRS(Acrs->getRowMap(), 0));
      Tpetra::MatrixMatrix::Multiply(*Acrs, transposeA, *Bcrs, transposeB, *tempAc, haveMultiplyDoFillComplete, label, new_params);

      // Temporary output matrix
      RCP<Tpetra::BlockCrsMatrix<SC, LO, GO, NO>> Ac_t       = Tpetra::convertToBlockCrsMatrix(*tempAc, A.GetStorageBlockSize());
      RCP<Xpetra::TpetraBlockCrsMatrix<SC, LO, GO, NO>> Ac_x = Teuchos::rcp(new Xpetra::TpetraBlockCrsMatrix<SC, LO, GO, NO>(Ac_t));
      RCP<Xpetra::CrsMatrix<SC, LO, GO, NO>> Ac_p            = Ac_x;

      // We can now cheat and replace the innards of Ac
      RCP<Xpetra::CrsMatrixWrap<SC, LO, GO, NO>> Ac_w = Teuchos::rcp_dynamic_cast<Xpetra::CrsMatrixWrap<SC, LO, GO, NO>>(Teuchos::rcpFromRef(C));
      Ac_w->replaceCrsMatrix(Ac_p);
    } else {
      // Mix and match
      TEUCHOS_TEST_FOR_EXCEPTION(1, Exceptions::RuntimeError, "Mix-and-match Crs/BlockCrs Multiply not currently supported");
    }
#else
    throw(Xpetra::Exceptions::RuntimeError("Xpetra must be compiled with Tpetra."));
#endif
  }

  if (call_FillComplete_on_result && !haveMultiplyDoFillComplete) {
    RCP<Teuchos::ParameterList> fillParams = rcp(new Teuchos::ParameterList());
    fillParams->set("Optimize Storage", doOptimizeStorage);
    C.fillComplete((transposeB) ? B.getRangeMap() : B.getDomainMap(),
                   (transposeA) ? A.getDomainMap() : A.getRangeMap(),
                   fillParams);
  }

  // transfer striding information
  RCP<Matrix> rcpA = Teuchos::rcp_const_cast<Matrix>(Teuchos::rcpFromRef(A));
  RCP<Matrix> rcpB = Teuchos::rcp_const_cast<Matrix>(Teuchos::rcpFromRef(B));
  C.CreateView("stridedMaps", rcpA, transposeA, rcpB, transposeB);  // TODO use references instead of RCPs
}  // end Multiply

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
RCP<Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>> MatrixMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Multiply(const Matrix& A, bool transposeA, const Matrix& B, bool transposeB, RCP<Matrix> C_in,
                                                                                                                                 Teuchos::FancyOStream& fos,
                                                                                                                                 bool doFillComplete,
                                                                                                                                 bool doOptimizeStorage,
                                                                                                                                 const std::string& label,
                                                                                                                                 const RCP<ParameterList>& params) {
  TEUCHOS_TEST_FOR_EXCEPTION(!A.isFillComplete(), Exceptions::RuntimeError, "A is not fill-completed");
  TEUCHOS_TEST_FOR_EXCEPTION(!B.isFillComplete(), Exceptions::RuntimeError, "B is not fill-completed");

  // Default case: Xpetra Multiply
  RCP<Matrix> C = C_in;

  if (C == Teuchos::null) {
    double nnzPerRow = Teuchos::as<double>(0);

#if 0
        if (A.getDomainMap()->lib() == Xpetra::UseTpetra) {
          // For now, follow what ML and Epetra do.
          GO numRowsA = A.getGlobalNumRows();
          GO numRowsB = B.getGlobalNumRows();
          nnzPerRow = sqrt(Teuchos::as<double>(A.getGlobalNumEntries())/numRowsA) +
              sqrt(Teuchos::as<double>(B.getGlobalNumEntries())/numRowsB) - 1;
          nnzPerRow *=  nnzPerRow;
          double totalNnz = nnzPerRow * A.getGlobalNumRows() * 0.75 + 100;
          double minNnz = Teuchos::as<double>(1.2 * A.getGlobalNumEntries());
          if (totalNnz < minNnz)
            totalNnz = minNnz;
          nnzPerRow = totalNnz / A.getGlobalNumRows();

          fos << "Matrix product nnz per row estimate = " << Teuchos::as<LO>(nnzPerRow) << std::endl;
        }
#endif
    if (transposeA)
      C = MatrixFactory::Build(A.getDomainMap(), Teuchos::as<LO>(nnzPerRow));
    else
      C = MatrixFactory::Build(A.getRowMap(), Teuchos::as<LO>(nnzPerRow));

  } else {
    C->resumeFill();  // why this is not done inside of Tpetra MxM?
    fos << "Reuse C pattern" << std::endl;
  }

  Multiply(A, transposeA, B, transposeB, *C, doFillComplete, doOptimizeStorage, label, params);  // call Multiply routine from above

  return C;
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
RCP<Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>> MatrixMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Multiply(const Matrix& A, bool transposeA, const Matrix& B, bool transposeB, Teuchos::FancyOStream& fos,
                                                                                                                                 bool callFillCompleteOnResult, bool doOptimizeStorage, const std::string& label,
                                                                                                                                 const RCP<ParameterList>& params) {
  return Multiply(A, transposeA, B, transposeB, Teuchos::null, fos, callFillCompleteOnResult, doOptimizeStorage, label, params);
}

#ifdef HAVE_XPETRA_EPETRAEXT
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
RCP<Epetra_CrsMatrix> MatrixMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::MLTwoMatrixMultiply(const Epetra_CrsMatrix& epA,
                                                                                                   const Epetra_CrsMatrix& epB,
                                                                                                   Teuchos::FancyOStream& fos) {
  throw(Xpetra::Exceptions::RuntimeError("MLTwoMatrixMultiply only available for GO=int or GO=long long with EpetraNode (Serial or OpenMP depending on configuration)"));
  TEUCHOS_UNREACHABLE_RETURN(Teuchos::null);
}
#endif  // ifdef HAVE_XPETRA_EPETRAEXT

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
RCP<Xpetra::BlockedCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>> MatrixMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::TwoMatrixMultiplyBlock(const BlockedCrsMatrix& A, bool transposeA,
                                                                                                                                                         const BlockedCrsMatrix& B, bool transposeB,
                                                                                                                                                         Teuchos::FancyOStream& fos,
                                                                                                                                                         bool doFillComplete,
                                                                                                                                                         bool doOptimizeStorage) {
  TEUCHOS_TEST_FOR_EXCEPTION(transposeA || transposeB, Exceptions::RuntimeError,
                             "TwoMatrixMultiply for BlockedCrsMatrix not implemented for transposeA==true or transposeB==true");

  // Preconditions
  TEUCHOS_TEST_FOR_EXCEPTION(!A.isFillComplete(), Exceptions::RuntimeError, "A is not fill-completed");
  TEUCHOS_TEST_FOR_EXCEPTION(!B.isFillComplete(), Exceptions::RuntimeError, "B is not fill-completed");

  RCP<const MapExtractor> rgmapextractor = A.getRangeMapExtractor();
  RCP<const MapExtractor> domapextractor = B.getDomainMapExtractor();

  RCP<BlockedCrsMatrix> C = rcp(new BlockedCrsMatrix(rgmapextractor, domapextractor, 33 /* TODO fix me */));

  for (size_t i = 0; i < A.Rows(); ++i) {    // loop over all block rows of A
    for (size_t j = 0; j < B.Cols(); ++j) {  // loop over all block columns of B
      RCP<Matrix> Cij;

      for (size_t l = 0; l < B.Rows(); ++l) {  // loop for calculating entry C_{ij}
        RCP<Matrix> crmat1 = A.getMatrix(i, l);
        RCP<Matrix> crmat2 = B.getMatrix(l, j);

        if (crmat1.is_null() || crmat2.is_null())
          continue;

        // try unwrapping 1x1 blocked matrix
        {
          auto unwrap1 = Teuchos::rcp_dynamic_cast<BlockedCrsMatrix>(crmat1);
          auto unwrap2 = Teuchos::rcp_dynamic_cast<BlockedCrsMatrix>(crmat2);

          if (unwrap1.is_null() != unwrap2.is_null()) {
            if (unwrap1 != Teuchos::null && unwrap1->Rows() == 1 && unwrap1->Cols() == 1)
              crmat1 = unwrap1->getCrsMatrix();
            if (unwrap2 != Teuchos::null && unwrap2->Rows() == 1 && unwrap2->Cols() == 1)
              crmat2 = unwrap2->getCrsMatrix();
          }
        }

        RCP<CrsMatrixWrap> crop1 = Teuchos::rcp_dynamic_cast<CrsMatrixWrap>(crmat1);
        RCP<CrsMatrixWrap> crop2 = Teuchos::rcp_dynamic_cast<CrsMatrixWrap>(crmat2);
        TEUCHOS_TEST_FOR_EXCEPTION(crop1.is_null() != crop2.is_null(), Xpetra::Exceptions::RuntimeError,
                                   "A and B must be either both (compatible) BlockedCrsMatrix objects or both CrsMatrixWrap objects.");

        // Forcibly compute the global constants if we don't have them (only works for real CrsMatrices, not nested blocks)
        if (!crop1.is_null())
          Teuchos::rcp_const_cast<CrsGraph>(crmat1->getCrsGraph())->computeGlobalConstants();
        if (!crop2.is_null())
          Teuchos::rcp_const_cast<CrsGraph>(crmat2->getCrsGraph())->computeGlobalConstants();

        TEUCHOS_TEST_FOR_EXCEPTION(!crmat1->haveGlobalConstants(), Exceptions::RuntimeError,
                                   "crmat1 does not have global constants");
        TEUCHOS_TEST_FOR_EXCEPTION(!crmat2->haveGlobalConstants(), Exceptions::RuntimeError,
                                   "crmat2 does not have global constants");

        if (crmat1->getGlobalNumEntries() == 0 || crmat2->getGlobalNumEntries() == 0)
          continue;

        // temporary matrix containing result of local block multiplication
        RCP<Matrix> temp = Teuchos::null;

        if (crop1 != Teuchos::null && crop2 != Teuchos::null)
          temp = Multiply(*crop1, false, *crop2, false, fos);
        else {
          RCP<BlockedCrsMatrix> bop1 = Teuchos::rcp_dynamic_cast<BlockedCrsMatrix>(crmat1);
          RCP<BlockedCrsMatrix> bop2 = Teuchos::rcp_dynamic_cast<BlockedCrsMatrix>(crmat2);
          TEUCHOS_TEST_FOR_EXCEPTION(bop1.is_null() == true, Xpetra::Exceptions::BadCast, "A is not a BlockedCrsMatrix. (TwoMatrixMultiplyBlock)");
          TEUCHOS_TEST_FOR_EXCEPTION(bop2.is_null() == true, Xpetra::Exceptions::BadCast, "B is not a BlockedCrsMatrix. (TwoMatrixMultiplyBlock)");
          TEUCHOS_TEST_FOR_EXCEPTION(bop1->Cols() != bop2->Rows(), Xpetra::Exceptions::RuntimeError, "A has " << bop1->Cols() << " columns and B has " << bop2->Rows() << " rows. Matrices are not compatible! (TwoMatrixMultiplyBlock)");
          TEUCHOS_TEST_FOR_EXCEPTION(bop1->getDomainMap()->isSameAs(*(bop2->getRangeMap())) == false, Xpetra::Exceptions::RuntimeError, "Domain map of A is not the same as range map of B. Matrices are not compatible! (TwoMatrixMultiplyBlock)");

          // recursive multiplication call
          temp = TwoMatrixMultiplyBlock(*bop1, transposeA, *bop2, transposeB, fos, doFillComplete, doOptimizeStorage);

          RCP<BlockedCrsMatrix> btemp = Teuchos::rcp_dynamic_cast<BlockedCrsMatrix>(temp);
          TEUCHOS_TEST_FOR_EXCEPTION(btemp->Rows() != bop1->Rows(), Xpetra::Exceptions::RuntimeError, "Number of block rows of local blocked operator is " << btemp->Rows() << " but should be " << bop1->Rows() << ". (TwoMatrixMultiplyBlock)");
          TEUCHOS_TEST_FOR_EXCEPTION(btemp->Cols() != bop2->Cols(), Xpetra::Exceptions::RuntimeError, "Number of block cols of local blocked operator is " << btemp->Cols() << " but should be " << bop2->Cols() << ". (TwoMatrixMultiplyBlock)");
          TEUCHOS_TEST_FOR_EXCEPTION(btemp->getRangeMapExtractor()->getFullMap()->isSameAs(*(bop1->getRangeMapExtractor()->getFullMap())) == false, Xpetra::Exceptions::RuntimeError, "Range map of local blocked operator should be same as first operator. (TwoMatrixMultiplyBlock)");
          TEUCHOS_TEST_FOR_EXCEPTION(btemp->getDomainMapExtractor()->getFullMap()->isSameAs(*(bop2->getDomainMapExtractor()->getFullMap())) == false, Xpetra::Exceptions::RuntimeError, "Domain map of local blocked operator should be same as second operator. (TwoMatrixMultiplyBlock)");
          TEUCHOS_TEST_FOR_EXCEPTION(btemp->getRangeMapExtractor()->getThyraMode() != bop1->getRangeMapExtractor()->getThyraMode(), Xpetra::Exceptions::RuntimeError, "Thyra mode of local range map extractor incompatible with range map extractor of A (TwoMatrixMultiplyBlock)");
          TEUCHOS_TEST_FOR_EXCEPTION(btemp->getDomainMapExtractor()->getThyraMode() != bop2->getDomainMapExtractor()->getThyraMode(), Xpetra::Exceptions::RuntimeError, "Thyra mode of local domain map extractor incompatible with domain map extractor of B (TwoMatrixMultiplyBlock)");
        }

        TEUCHOS_TEST_FOR_EXCEPTION(temp->isFillComplete() == false, Xpetra::Exceptions::RuntimeError, "Local block is not filled. (TwoMatrixMultiplyBlock)");

        RCP<Matrix> addRes = null;
        if (Cij.is_null())
          Cij = temp;
        else {
          MatrixMatrix::TwoMatrixAdd(*temp, false, 1.0, *Cij, false, 1.0, addRes, fos);
          Cij = addRes;
        }
      }

      if (!Cij.is_null()) {
        if (Cij->isFillComplete())
          Cij->resumeFill();
        Cij->fillComplete(B.getDomainMap(j), A.getRangeMap(i));
        C->setMatrix(i, j, Cij);
      } else {
        C->setMatrix(i, j, Teuchos::null);
      }
    }
  }

  if (doFillComplete)
    C->fillComplete();  // call default fillComplete for BlockCrsMatrixWrap objects

  return C;
}  // TwoMatrixMultiplyBlock

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void MatrixMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::TwoMatrixAdd(const Matrix& A, bool transposeA, SC alpha, Matrix& B, SC beta) {
  if (!(A.getRowMap()->isSameAs(*(B.getRowMap()))))
    throw Exceptions::Incompatible("TwoMatrixAdd: matrix row maps are not the same.");

  if (A.getRowMap()->lib() == Xpetra::UseEpetra) {
    throw Exceptions::RuntimeError("TwoMatrixAdd for Epetra matrices needs <double,int,int> for Scalar, LocalOrdinal and GlobalOrdinal.");
  } else if (A.getRowMap()->lib() == Xpetra::UseTpetra) {
#ifdef HAVE_XPETRA_TPETRA
    const Tpetra::CrsMatrix<SC, LO, GO, NO>& tpA = Xpetra::Helpers<SC, LO, GO, NO>::Op2TpetraCrs(A);
    Tpetra::CrsMatrix<SC, LO, GO, NO>& tpB       = Xpetra::Helpers<SC, LO, GO, NO>::Op2NonConstTpetraCrs(B);

    Tpetra::MatrixMatrix::Add(tpA, transposeA, alpha, tpB, beta);
#else
    throw Exceptions::RuntimeError("Xpetra must be compiled with Tpetra.");
#endif
  }
}  // MatrixMatrix::TwoMatrixAdd()

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void MatrixMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::TwoMatrixAdd(const Matrix& A, bool transposeA, const SC& alpha,
                                                                           const Matrix& B, bool transposeB, const SC& beta,
                                                                           RCP<Matrix>& C, Teuchos::FancyOStream& fos, bool AHasFixedNnzPerRow) {
  RCP<const Matrix> rcpA              = Teuchos::rcpFromRef(A);
  RCP<const Matrix> rcpB              = Teuchos::rcpFromRef(B);
  RCP<const BlockedCrsMatrix> rcpBopA = Teuchos::rcp_dynamic_cast<const BlockedCrsMatrix>(rcpA);
  RCP<const BlockedCrsMatrix> rcpBopB = Teuchos::rcp_dynamic_cast<const BlockedCrsMatrix>(rcpB);
  // We have to distinguish 4 cases:
  // both matrices are CrsMatrixWrap based, only one of them is CrsMatrixWrap based, or none.

  // C can be null, so just use A to get the lib
  auto lib = A.getRowMap()->lib();

  // Both matrices are CrsMatrixWrap
  if (rcpBopA == Teuchos::null && rcpBopB == Teuchos::null) {
    if (!(A.getRowMap()->isSameAs(*(B.getRowMap()))))
      throw Exceptions::Incompatible("TwoMatrixAdd: matrix row maps are not the same.");
    if (lib == Xpetra::UseEpetra) {
      throw Exceptions::RuntimeError("MatrixMatrix::Add for Epetra only available with Scalar = double, LO = GO = int.");
    } else if (lib == Xpetra::UseTpetra) {
#ifdef HAVE_XPETRA_TPETRA
      using tcrs_matrix_type      = Tpetra::CrsMatrix<SC, LO, GO, NO>;
      using helpers               = Xpetra::Helpers<SC, LO, GO, NO>;
      const tcrs_matrix_type& tpA = helpers::Op2TpetraCrs(A);
      const tcrs_matrix_type& tpB = helpers::Op2TpetraCrs(B);
      C                           = helpers::tpetraAdd(tpA, transposeA, alpha, tpB, transposeB, beta);
#else
      throw Exceptions::RuntimeError("Xpetra must be compiled with Tpetra.");
#endif
    }
    ///////////////////////// EXPERIMENTAL
    if (A.IsView("stridedMaps")) C->CreateView("stridedMaps", rcpFromRef(A));
    if (B.IsView("stridedMaps")) C->CreateView("stridedMaps", rcpFromRef(B));
    ///////////////////////// EXPERIMENTAL
  }
  // the first matrix is of type CrsMatrixWrap, the second is a blocked operator
  else if (rcpBopA == Teuchos::null && rcpBopB != Teuchos::null) {
    RCP<const MapExtractor> rgmapextractor = rcpBopB->getRangeMapExtractor();
    RCP<const MapExtractor> domapextractor = rcpBopB->getDomainMapExtractor();

    C                        = rcp(new BlockedCrsMatrix(rgmapextractor, domapextractor, 33 /* TODO fix me */));
    RCP<BlockedCrsMatrix> bC = Teuchos::rcp_dynamic_cast<BlockedCrsMatrix>(C);

    size_t i = 0;
    for (size_t j = 0; j < rcpBopB->Cols(); ++j) {  // loop over all block columns of B
      RCP<Matrix> Cij = Teuchos::null;
      if (rcpA != Teuchos::null &&
          rcpBopB->getMatrix(i, j) != Teuchos::null) {
        // recursive call
        TwoMatrixAdd(*rcpA, transposeA, alpha,
                     *(rcpBopB->getMatrix(i, j)), transposeB, beta,
                     Cij, fos, AHasFixedNnzPerRow);
      } else if (rcpA == Teuchos::null &&
                 rcpBopB->getMatrix(i, j) != Teuchos::null) {
        Cij = rcpBopB->getMatrix(i, j);
      } else if (rcpA != Teuchos::null &&
                 rcpBopB->getMatrix(i, j) == Teuchos::null) {
        Cij = Teuchos::rcp_const_cast<Matrix>(rcpA);
      } else {
        Cij = Teuchos::null;
      }

      if (!Cij.is_null()) {
        if (Cij->isFillComplete())
          Cij->resumeFill();
        Cij->fillComplete();
        bC->setMatrix(i, j, Cij);
      } else {
        bC->setMatrix(i, j, Teuchos::null);
      }
    }  // loop over columns j
  }
  // the second matrix is of type CrsMatrixWrap, the first is a blocked operator
  else if (rcpBopA != Teuchos::null && rcpBopB == Teuchos::null) {
    RCP<const MapExtractor> rgmapextractor = rcpBopA->getRangeMapExtractor();
    RCP<const MapExtractor> domapextractor = rcpBopA->getDomainMapExtractor();

    C                        = rcp(new BlockedCrsMatrix(rgmapextractor, domapextractor, 33 /* TODO fix me */));
    RCP<BlockedCrsMatrix> bC = Teuchos::rcp_dynamic_cast<BlockedCrsMatrix>(C);
    size_t j                 = 0;
    for (size_t i = 0; i < rcpBopA->Rows(); ++i) {  // loop over all block rows of A
      RCP<Matrix> Cij = Teuchos::null;
      if (rcpBopA->getMatrix(i, j) != Teuchos::null &&
          rcpB != Teuchos::null) {
        // recursive call
        TwoMatrixAdd(*(rcpBopA->getMatrix(i, j)), transposeA, alpha,
                     *rcpB, transposeB, beta,
                     Cij, fos, AHasFixedNnzPerRow);
      } else if (rcpBopA->getMatrix(i, j) == Teuchos::null &&
                 rcpB != Teuchos::null) {
        Cij = Teuchos::rcp_const_cast<Matrix>(rcpB);
      } else if (rcpBopA->getMatrix(i, j) != Teuchos::null &&
                 rcpB == Teuchos::null) {
        Cij = rcpBopA->getMatrix(i, j);
      } else {
        Cij = Teuchos::null;
      }

      if (!Cij.is_null()) {
        if (Cij->isFillComplete())
          Cij->resumeFill();
        Cij->fillComplete();
        bC->setMatrix(i, j, Cij);
      } else {
        bC->setMatrix(i, j, Teuchos::null);
      }
    }  // loop over rows i
  } else {
    // This is the version for blocked matrices

    // check the compatibility of the blocked operators
    TEUCHOS_TEST_FOR_EXCEPTION(rcpBopA.is_null() == true, Xpetra::Exceptions::BadCast, "A is not a BlockedCrsMatrix. (TwoMatrixAdd)");
    TEUCHOS_TEST_FOR_EXCEPTION(rcpBopB.is_null() == true, Xpetra::Exceptions::BadCast, "B is not a BlockedCrsMatrix. (TwoMatrixAdd)");
    TEUCHOS_TEST_FOR_EXCEPTION(rcpBopA->Rows() != rcpBopB->Rows(), Xpetra::Exceptions::RuntimeError, "A has " << rcpBopA->Rows() << " rows and B has " << rcpBopA->Rows() << " rows. Matrices are not compatible! (TwoMatrixAdd)");
    TEUCHOS_TEST_FOR_EXCEPTION(rcpBopA->Rows() != rcpBopB->Rows(), Xpetra::Exceptions::RuntimeError, "A has " << rcpBopA->Cols() << " cols and B has " << rcpBopA->Cols() << " cols. Matrices are not compatible! (TwoMatrixAdd)");
    TEUCHOS_TEST_FOR_EXCEPTION(rcpBopA->getRangeMap()->isSameAs(*(rcpBopB->getRangeMap())) == false, Xpetra::Exceptions::RuntimeError, "Range map of A is not the same as range map of B. Matrices are not compatible! (TwoMatrixAdd)");
    TEUCHOS_TEST_FOR_EXCEPTION(rcpBopA->getDomainMap()->isSameAs(*(rcpBopB->getDomainMap())) == false, Xpetra::Exceptions::RuntimeError, "Domain map of A is not the same as domain map of B. Matrices are not compatible! (TwoMatrixAdd)");
    TEUCHOS_TEST_FOR_EXCEPTION(rcpBopA->getRangeMapExtractor()->getThyraMode() != rcpBopB->getRangeMapExtractor()->getThyraMode(), Xpetra::Exceptions::RuntimeError, "Different Thyra/Xpetra style gids in RangeMapExtractor of A and B. Matrices are not compatible! (TwoMatrixAdd)");
    TEUCHOS_TEST_FOR_EXCEPTION(rcpBopA->getDomainMapExtractor()->getThyraMode() != rcpBopB->getDomainMapExtractor()->getThyraMode(), Xpetra::Exceptions::RuntimeError, "Different Thyra/Xpetra style gids in DomainMapExtractor of A and B. Matrices are not compatible! (TwoMatrixAdd)");

    RCP<const MapExtractor> rgmapextractor = rcpBopA->getRangeMapExtractor();
    RCP<const MapExtractor> domapextractor = rcpBopB->getDomainMapExtractor();

    C                        = rcp(new BlockedCrsMatrix(rgmapextractor, domapextractor, 33 /* TODO fix me */));
    RCP<BlockedCrsMatrix> bC = Teuchos::rcp_dynamic_cast<BlockedCrsMatrix>(C);
    for (size_t i = 0; i < rcpBopA->Rows(); ++i) {    // loop over all block rows of A
      for (size_t j = 0; j < rcpBopB->Cols(); ++j) {  // loop over all block columns of B

        RCP<Matrix> Cij = Teuchos::null;
        if (rcpBopA->getMatrix(i, j) != Teuchos::null &&
            rcpBopB->getMatrix(i, j) != Teuchos::null) {
          // recursive call
          TwoMatrixAdd(*(rcpBopA->getMatrix(i, j)), transposeA, alpha,
                       *(rcpBopB->getMatrix(i, j)), transposeB, beta,
                       Cij, fos, AHasFixedNnzPerRow);
        } else if (rcpBopA->getMatrix(i, j) == Teuchos::null &&
                   rcpBopB->getMatrix(i, j) != Teuchos::null) {
          Cij = rcpBopB->getMatrix(i, j);
        } else if (rcpBopA->getMatrix(i, j) != Teuchos::null &&
                   rcpBopB->getMatrix(i, j) == Teuchos::null) {
          Cij = rcpBopA->getMatrix(i, j);
        } else {
          Cij = Teuchos::null;
        }

        if (!Cij.is_null()) {
          if (Cij->isFillComplete())
            Cij->resumeFill();
          Cij->fillComplete();
          bC->setMatrix(i, j, Cij);
        } else {
          bC->setMatrix(i, j, Teuchos::null);
        }
      }  // loop over columns j
    }    // loop over rows i

  }  // end blocked recursive algorithm
}  // MatrixMatrix::TwoMatrixAdd()

}  // end namespace Xpetra

#endif /* PACKAGES_XPETRA_SUP_UTILS_XPETRA_MATRIXMATRIX_DEF_HPP_ */
