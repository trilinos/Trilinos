// @HEADER
// *****************************************************************************
//        MueLu: A package for multigrid based preconditioning
//
// Copyright 2012 NTESS and the MueLu contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef MUELU_CGSOLVER_DEF_HPP
#define MUELU_CGSOLVER_DEF_HPP

#include <Xpetra_MatrixFactory.hpp>
#include <Xpetra_MatrixFactory2.hpp>
#include <Xpetra_MatrixMatrix.hpp>

#include "MueLu_Utilities.hpp"
#include "MueLu_Constraint.hpp"
#include "MueLu_Monitor.hpp"

#include "MueLu_CGSolver_decl.hpp"

namespace MueLu {

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
CGSolver<Scalar, LocalOrdinal, GlobalOrdinal, Node>::CGSolver(size_t Its)
  : nIts_(Its) {}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void CGSolver<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Iterate(const Matrix& Aref, const Constraint& C, const Matrix& P0, RCP<Matrix>& finalP) const {
  // Note: this function matrix notations follow Saad's "Iterative methods", ed. 2, pg. 246
  // So, X is the unknown prolongator, P's are conjugate directions, Z's are preconditioned P's
  PrintMonitor m(*this, "CG iterations");

  if (nIts_ == 0) {
    finalP = MatrixFactory2::BuildCopy(rcpFromRef(P0));
    return;
  }

  RCP<const Matrix> A = rcpFromRef(Aref);
  const auto rowMap   = A->getRowMap();
  auto invDiag        = Xpetra::VectorFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Build(rowMap, true);
  A->getLocalDiagCopy(*invDiag);
  invDiag->reciprocal(*invDiag);

  bool useTpetra = (A->getRowMap()->lib() == Xpetra::UseTpetra);

  Teuchos::FancyOStream& mmfancy = this->GetOStream(Statistics2);

  SC one = Teuchos::ScalarTraits<SC>::one();

  RCP<Matrix> X, P, R, Z, AP;
  RCP<Matrix> newX, tmpAP;
#ifndef TWO_ARG_MATRIX_ADD
  RCP<Matrix> newR, newP;
#endif

  SC oldRZ, newRZ, alpha, beta, app;

  // T is used only for projecting onto
  RCP<CrsMatrix> T_ = CrsMatrixFactory::Build(C.GetPattern());
  T_->fillComplete(P0.getDomainMap(), P0.getRangeMap());
  RCP<Matrix> T = rcp(new CrsMatrixWrap(T_));

  // Initial P0 would only be used for multiplication
  X = rcp_const_cast<Matrix>(rcpFromRef(P0));

  tmpAP = MatrixMatrix::Multiply(*A, false, *X, false, mmfancy, true /*doFillComplete*/, true /*optimizeStorage*/);
  C.Apply(*tmpAP, *T);

  // R_0 = -A*X_0
  R = Xpetra::MatrixFactory2<Scalar, LocalOrdinal, GlobalOrdinal, Node>::BuildCopy(T);

  R->scale(-one);

  // Z_0 = M^{-1}R_0
  Z = Xpetra::MatrixFactory2<Scalar, LocalOrdinal, GlobalOrdinal, Node>::BuildCopy(R);
  Z->leftScale(*invDiag);

  // P_0 = Z_0
  P = Xpetra::MatrixFactory2<Scalar, LocalOrdinal, GlobalOrdinal, Node>::BuildCopy(Z);

  oldRZ = Utilities::Frobenius(*R, *Z);

  for (size_t i = 0; i < nIts_; i++) {
    // AP = constrain(A*P)
    if (i == 0 || useTpetra) {
      // Construct the MxM pattern from scratch
      // This is done by default for Tpetra as the three argument version requires tmpAP
      // to *not* be locally indexed which defeats the purpose
      // TODO: need a three argument Tpetra version which allows reuse of already fill-completed matrix
      tmpAP = MatrixMatrix::Multiply(*A, false, *P, false, mmfancy, true /*doFillComplete*/, true /*optimizeStorage*/);
    } else {
      // Reuse the MxM pattern
      tmpAP = MatrixMatrix::Multiply(*A, false, *P, false, tmpAP, mmfancy, true /*doFillComplete*/, true /*optimizeStorage*/);
    }
    C.Apply(*tmpAP, *T);
    AP = T;

    app = Utilities::Frobenius(*AP, *P);
    if (Teuchos::ScalarTraits<SC>::magnitude(app) < Teuchos::ScalarTraits<SC>::sfmin()) {
      // It happens, for instance, if P = 0
      // For example, if we use TentativePFactory for both nonzero pattern and initial guess
      // I think it might also happen because of numerical breakdown, but we don't test for that yet
      if (i == 0)
        X = MatrixFactory2::BuildCopy(rcpFromRef(P0));
      break;
    }

    // alpha = (R_i, Z_i)/(A*P_i, P_i)
    alpha = oldRZ / app;
    this->GetOStream(Runtime1, 1) << "alpha = " << alpha << std::endl;

    // X_{i+1} = X_i + alpha*P_i
#ifndef TWO_ARG_MATRIX_ADD
    newX = Teuchos::null;
    MatrixMatrix::TwoMatrixAdd(*P, false, alpha, *X, false, one, newX, mmfancy);
    newX->fillComplete(P0.getDomainMap(), P0.getRangeMap());
    X.swap(newX);
#else
    MatrixMatrix::TwoMatrixAdd(*P, false, alpha, *X, one);
#endif

    if (i == nIts_ - 1)
      break;

      // R_{i+1} = R_i - alpha*A*P_i
#ifndef TWO_ARG_MATRIX_ADD
    newR = Teuchos::null;
    MatrixMatrix::TwoMatrixAdd(*AP, false, -alpha, *R, false, one, newR, mmfancy);
    newR->fillComplete(P0.getDomainMap(), P0.getRangeMap());
    R.swap(newR);
#else
    MatrixMatrix::TwoMatrixAdd(*AP, false, -alpha, *R, one);
#endif

    // Z_{i+1} = M^{-1} R_{i+1}
    Z = MatrixFactory2::BuildCopy(R);
    Z->leftScale(*invDiag);

    // beta = (R_{i+1}, Z_{i+1})/(R_i, Z_i)
    newRZ = Utilities::Frobenius(*R, *Z);
    beta  = newRZ / oldRZ;

    // P_{i+1} = Z_{i+1} + beta*P_i
#ifndef TWO_ARG_MATRIX_ADD
    newP = Teuchos::null;
    MatrixMatrix::TwoMatrixAdd(*P, false, beta, *Z, false, one, newP, mmfancy);
    newP->fillComplete(P0.getDomainMap(), P0.getRangeMap());
    P.swap(newP);
#else
    MatrixMatrix::TwoMatrixAdd(*Z, false, one, *P, beta);
#endif

    oldRZ = newRZ;
  }

  finalP = X;
}

}  // namespace MueLu

#endif  // ifndef MUELU_CGSOLVER_DECL_HPP
