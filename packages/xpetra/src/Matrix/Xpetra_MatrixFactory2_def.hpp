// @HEADER
// *****************************************************************************
//             Xpetra: A linear algebra interface package
//
// Copyright 2012 NTESS and the Xpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

// WARNING: This code is experimental. Backwards compatibility should not be expected.

#ifndef XPETRA_MATRIXFACTORY2_DEF_HPP
#define XPETRA_MATRIXFACTORY2_DEF_HPP

#include "Xpetra_MatrixFactory2_decl.hpp"
#include "Xpetra_BlockedCrsMatrix.hpp"

namespace Xpetra {

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
RCP<Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>> MatrixFactory2<Scalar, LocalOrdinal, GlobalOrdinal, Node>::BuildCopy(const RCP<const Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>> A, bool setFixedBlockSize) {
  RCP<const CrsMatrixWrap> oldOp = Teuchos::rcp_dynamic_cast<const CrsMatrixWrap>(A);
  if (oldOp == Teuchos::null)
    throw Exceptions::BadCast("Cast from Xpetra::Matrix to Xpetra::CrsMatrixWrap failed");

  RCP<const CrsMatrix> oldCrsOp = oldOp->getCrsMatrix();

  UnderlyingLib lib = A->getRowMap()->lib();

  TEUCHOS_TEST_FOR_EXCEPTION(lib != UseEpetra && lib != UseTpetra, Exceptions::RuntimeError,
                             "Not Epetra or Tpetra matrix");

#ifdef HAVE_XPETRA_EPETRA
  if (lib == UseEpetra) {
    // NOTE: The proper Epetra conversion in Xpetra_MatrixFactory.cpp
    throw Exceptions::RuntimeError("Xpetra::BuildCopy(): matrix templates are incompatible with Epetra");
  }
#endif

#ifdef HAVE_XPETRA_TPETRA
  if (lib == UseTpetra) {
    // Underlying matrix is Tpetra
    RCP<const TpetraCrsMatrix> oldTCrsOp = Teuchos::rcp_dynamic_cast<const TpetraCrsMatrix>(oldCrsOp);

    if (oldTCrsOp != Teuchos::null) {
      RCP<TpetraCrsMatrix> newTCrsOp(new TpetraCrsMatrix(*oldTCrsOp));
      RCP<CrsMatrixWrap> newOp(new CrsMatrixWrap(Teuchos::as<RCP<CrsMatrix>>(newTCrsOp)));
      if (setFixedBlockSize)
        newOp->SetFixedBlockSize(A->GetFixedBlockSize());

      return newOp;
    } else {
      throw Exceptions::BadCast("Cast from Xpetra::Matrix to Xpetra::TpetraCrsMatrix failed");
    }
  }
#endif

  return Teuchos::null;
}

}  // namespace Xpetra

#endif
