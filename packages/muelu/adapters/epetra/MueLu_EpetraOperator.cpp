// @HEADER
// *****************************************************************************
//        MueLu: A package for multigrid based preconditioning
//
// Copyright 2012 NTESS and the MueLu contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include <Xpetra_Matrix.hpp>
#include <Xpetra_CrsMatrixWrap.hpp>
#include <Xpetra_BlockedCrsMatrix.hpp>
#include <Xpetra_EpetraMultiVector.hpp>

#include "MueLu_EpetraOperator.hpp"
#include "MueLu_Level.hpp"

#if defined(HAVE_MUELU_SERIAL) and defined(HAVE_MUELU_EPETRA)

namespace MueLu {

int EpetraOperator::ApplyInverse(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const {
  try {
    // There is no rcpFromRef(const T&), so we need to do const_cast
    const Xpetra::EpetraMultiVectorT<GO, NO> eX(rcpFromRef(const_cast<Epetra_MultiVector&>(X)));
    Xpetra::EpetraMultiVectorT<GO, NO> eY(rcpFromRef(Y));
    // Generally, we assume two different vectors, but AztecOO uses a single vector
    if (X.Values() == Y.Values()) {
      // X and Y point to the same memory, use an additional vector
      RCP<Xpetra::EpetraMultiVectorT<GO, NO>> tmpY = Teuchos::rcp(new Xpetra::EpetraMultiVectorT<GO, NO>(eY.getMap(), eY.getNumVectors()));
      // InitialGuessIsZero in MueLu::Hierarchy.Iterate() does not zero out components, it
      // only assumes that user provided an already zeroed out vector
      bool initialGuessZero = true;
      tmpY->putScalar(0.0);
      // apply one V-cycle as preconditioner
      Hierarchy_->Iterate(eX, *tmpY, 1, initialGuessZero);
      // deep copy solution from MueLu
      eY.update(1.0, *tmpY, 0.0);
    } else {
      // X and Y point to different memory, pass the vectors through

      // InitialGuessIsZero in MueLu::Hierarchy.Iterate() does not zero out components, it
      // only assumes that user provided an already zeroed out vector
      bool initialGuessZero = true;
      eY.putScalar(0.0);
      Hierarchy_->Iterate(eX, eY, 1, initialGuessZero);
    }

  } catch (std::exception& e) {
    // TODO: error msg directly on std::cerr?
    std::cerr << "Caught an exception in MueLu::EpetraOperator::ApplyInverse():" << std::endl
              << e.what() << std::endl;
    return -1;
  }
  return 0;
}

const Epetra_Comm& EpetraOperator::Comm() const {
  RCP<Matrix> A = Hierarchy_->GetLevel(0)->Get<RCP<Matrix>>("A");

  // TODO: This code is not pretty
  RCP<Xpetra::BlockedCrsMatrix<SC, LO, GO, NO>> epbA = Teuchos::rcp_dynamic_cast<Xpetra::BlockedCrsMatrix<SC, LO, GO, NO>>(A);
  if (epbA != Teuchos::null) {
    RCP<const Xpetra::Matrix<SC, LO, GO, NO>> blockMat            = epbA->getMatrix(0, 0);
    RCP<const Xpetra::CrsMatrixWrap<SC, LO, GO, NO>> blockCrsWrap = Teuchos::rcp_dynamic_cast<const Xpetra::CrsMatrixWrap<SC, LO, GO, NO>>(blockMat);
    if (blockCrsWrap == Teuchos::null)
      throw Exceptions::BadCast("MueLu::EpetraOperator::Comm(): Cast from block (0,0) to CrsMatrixWrap failed. Could be a block matrix. TODO implement recursive support for block matrices.");
    RCP<const Xpetra::EpetraCrsMatrixT<GO, NO>> tmp_ECrsMtx = rcp_dynamic_cast<const Xpetra::EpetraCrsMatrixT<GO, NO>>(blockCrsWrap->getCrsMatrix());
    if (tmp_ECrsMtx == Teuchos::null)
      throw Exceptions::BadCast("MueLu::EpetraOperator::Comm(): Cast from Xpetra::CrsMatrix to Xpetra::EpetraCrsMatrix failed");
    RCP<Epetra_CrsMatrix> epA = tmp_ECrsMtx->getEpetra_CrsMatrixNonConst();
    return epA->Comm();
  }

  RCP<const Xpetra::CrsMatrixWrap<SC, LO, GO, NO>> crsOp = rcp_dynamic_cast<const Xpetra::CrsMatrixWrap<SC, LO, GO, NO>>(A);
  if (crsOp == Teuchos::null)
    throw Exceptions::BadCast("Cast from Xpetra::Matrix to Xpetra::CrsMatrixWrap failed");
  const RCP<const Xpetra::EpetraCrsMatrixT<GO, NO>>& tmp_ECrsMtx = rcp_dynamic_cast<const Xpetra::EpetraCrsMatrixT<GO, NO>>(crsOp->getCrsMatrix());
  if (tmp_ECrsMtx == Teuchos::null)
    throw Exceptions::BadCast("Cast from Xpetra::CrsMatrix to Xpetra::EpetraCrsMatrix failed");
  return tmp_ECrsMtx->getEpetra_CrsMatrixNonConst()->Comm();
}

const Epetra_Map& EpetraOperator::OperatorDomainMap() const {
  RCP<Xpetra::Matrix<SC, LO, GO, NO>> A = Hierarchy_->GetLevel(0)->Get<RCP<Matrix>>("A");

  RCP<Xpetra::BlockedCrsMatrix<SC, LO, GO, NO>> epbA = Teuchos::rcp_dynamic_cast<Xpetra::BlockedCrsMatrix<SC, LO, GO, NO>>(A);
  if (epbA != Teuchos::null)
    return Xpetra::toEpetra(epbA->getFullDomainMap());  // TODO check me

  RCP<const Xpetra::CrsMatrixWrap<SC, LO, GO, NO>> crsOp = rcp_dynamic_cast<const Xpetra::CrsMatrixWrap<SC, LO, GO, NO>>(A);
  if (crsOp == Teuchos::null)
    throw Exceptions::BadCast("Cast from Xpetra::Matrix to Xpetra::CrsMatrixWrap failed");
  const RCP<const Xpetra::EpetraCrsMatrixT<GO, NO>>& tmp_ECrsMtx = rcp_dynamic_cast<const Xpetra::EpetraCrsMatrixT<GO, NO>>(crsOp->getCrsMatrix());
  if (tmp_ECrsMtx == Teuchos::null)
    throw Exceptions::BadCast("Cast from Xpetra::CrsMatrix to Xpetra::EpetraCrsMatrix failed");
  return tmp_ECrsMtx->getEpetra_CrsMatrixNonConst()->DomainMap();
}

const Epetra_Map& EpetraOperator::OperatorRangeMap() const {
  RCP<Xpetra::Matrix<SC, LO, GO, NO>> A = Hierarchy_->GetLevel(0)->Get<RCP<Matrix>>("A");

  RCP<Xpetra::BlockedCrsMatrix<SC, LO, GO, NO>> epbA = Teuchos::rcp_dynamic_cast<Xpetra::BlockedCrsMatrix<SC, LO, GO, NO>>(A);
  if (epbA != Teuchos::null)
    return Xpetra::toEpetra(epbA->getFullRangeMap());

  RCP<const Xpetra::CrsMatrixWrap<SC, LO, GO, NO>> crsOp = rcp_dynamic_cast<const Xpetra::CrsMatrixWrap<SC, LO, GO, NO>>(A);
  if (crsOp == Teuchos::null)
    throw Exceptions::BadCast("Cast from Xpetra::Matrix to Xpetra::CrsMatrixWrap failed");
  const RCP<const Xpetra::EpetraCrsMatrixT<GO, NO>>& tmp_ECrsMtx = rcp_dynamic_cast<const Xpetra::EpetraCrsMatrixT<GO, NO>>(crsOp->getCrsMatrix());
  if (tmp_ECrsMtx == Teuchos::null)
    throw Exceptions::BadCast("Cast from Xpetra::CrsMatrix to Xpetra::EpetraCrsMatrix failed");
  return tmp_ECrsMtx->getEpetra_CrsMatrixNonConst()->RangeMap();
}

}  // namespace MueLu

#endif  // #if defined(HAVE_MUELU_SERIAL) and defined(HAVE_MUELU_EPETRA)
