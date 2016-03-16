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
#ifndef MUELU_CGSOLVER_DEF_HPP
#define MUELU_CGSOLVER_DEF_HPP

#include <Xpetra_MatrixFactory.hpp>
#include <Xpetra_MatrixMatrix.hpp>

#include "MueLu_Utilities.hpp"
#include "MueLu_Constraint.hpp"
#include "MueLu_Monitor.hpp"


#include "MueLu_CGSolver.hpp"



namespace MueLu {

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  CGSolver<Scalar, LocalOrdinal, GlobalOrdinal, Node>::CGSolver(size_t Its)
  : nIts_(Its)
  { }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void CGSolver<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Iterate(const Matrix& Aref, const Constraint& C, const Matrix& P0, RCP<Matrix>& finalP) const {
    PrintMonitor m(*this, "CG iterations");

    if (nIts_ == 0) {
      finalP = Xpetra::MatrixFactory2<Scalar, LocalOrdinal, GlobalOrdinal, Node>::BuildCopy(rcpFromRef(P0));
      return;
    }

    // Note: this function matrix notations follow Saad's "Iterative methods", ed. 2, pg. 246
    // So, X is the unknown prolongator, P's are conjugate directions, Z's are preconditioned P's
    RCP<const Matrix> A = rcpFromRef(Aref);

    RCP<Matrix> X, P, R, Z, AP;
    RCP<Matrix> newX, tmpAP;
#ifndef TWO_ARG_MATRIX_ADD
    RCP<Matrix> newR, newP;
#endif

    SC oldRZ, newRZ, alpha, beta, app;

    bool useTpetra = (A->getRowMap()->lib() == Xpetra::UseTpetra);

    Teuchos::FancyOStream& mmfancy = this->GetOStream(Statistics2);

    // T is used only for projecting onto
    RCP<CrsMatrix> T_ = Xpetra::CrsMatrixFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Build(C.GetPattern());
    T_->fillComplete(P0.getDomainMap(), P0.getRangeMap());
    RCP<Matrix>    T = rcp(new CrsMatrixWrap(T_));

    SC one = Teuchos::ScalarTraits<SC>::one();

    Teuchos::ArrayRCP<const SC> D = Utilities::GetMatrixDiagonal(*A);

    // Initial P0 would only be used for multiplication
    X = rcp_const_cast<Matrix>(rcpFromRef(P0));

    bool doFillComplete  = true;
    // bool optimizeStorage = false;
    bool optimizeStorage = true;

    tmpAP = MatrixMatrix::Multiply(*A, false, *X, false, mmfancy, doFillComplete, optimizeStorage);
    C.Apply(*tmpAP, *T);

    // R_0 = -A*X_0
    R = Xpetra::MatrixFactory2<Scalar, LocalOrdinal, GlobalOrdinal, Node>::BuildCopy(T);
#ifdef HAVE_MUELU_TPETRA
#if 0 //def HAVE_MUELU_TPETRA_INST_INT_INT
    // TAW: Oct 16 2015: MueLu::Utilities returns the Tpetra::CrsMatrix object which would not be instantiated!
    //                   Catching this in Op2NonConstTpetraCrs is not possible as this does not affect the return type
    //                   Tpetra::CrsMatrix!
    if (useTpetra)
      Utilities::Op2NonConstTpetraCrs(R)->resumeFill();
#else
    this->GetOStream(Warnings0) << "WARNING: MueLu_CGSolver: calling Xpetra::CrsMatrix::resumeFill instead of Tpetra::CrsMatrix::resumeFill. The results should be verified in this case." << std::endl;
    R->resumeFill();
#endif
#endif
    R->scale(-one);
    if (!R->isFillComplete())
      R->fillComplete(R->getDomainMap(), R->getRangeMap());

    // Z_0 = M^{-1}R_0
    Z = Xpetra::MatrixFactory2<Scalar, LocalOrdinal, GlobalOrdinal, Node>::BuildCopy(R);
    Utilities::MyOldScaleMatrix(*Z, D, true, true, false);

    // P_0 = Z_0
    P = Xpetra::MatrixFactory2<Scalar, LocalOrdinal, GlobalOrdinal, Node>::BuildCopy(Z);

    oldRZ = Utilities::Frobenius(*R, *Z);

    for (size_t k = 0; k < nIts_; k++) {
      // AP = constrain(A*P)
      if (k == 0 || useTpetra) {
        // Construct the MxM pattern from scratch
        // This is done by default for Tpetra as the three argument version requires tmpAP
        // to *not* be locally indexed which defeats the purpose
        // TODO: need a three argument Tpetra version which allows reuse of already fill-completed matrix
        tmpAP = MatrixMatrix::Multiply(*A, false, *P, false,        mmfancy, doFillComplete, optimizeStorage);
      } else {
        // Reuse the MxM pattern
        tmpAP = MatrixMatrix::Multiply(*A, false, *P, false, tmpAP, mmfancy, doFillComplete, optimizeStorage);
      }
      C.Apply(*tmpAP, *T);
      AP = T;

      app = Utilities::Frobenius(*AP, *P);
      if (Teuchos::ScalarTraits<SC>::magnitude(app) < Teuchos::ScalarTraits<SC>::sfmin()) {
        // It happens, for instance, if P = 0
        // For example, if we use TentativePFactory for both nonzero pattern and initial guess
        // I think it might also happen because of numerical breakdown, but we don't test for that yet
        if (k == 0)
          X = Xpetra::MatrixFactory2<Scalar, LocalOrdinal, GlobalOrdinal, Node>::BuildCopy(rcpFromRef(P0));
        break;
      }

      // alpha = (R_k, Z_k)/(A*P_k, P_k)
      alpha = oldRZ / app;
      this->GetOStream(Runtime1,1) << "alpha = " << alpha << std::endl;

      // X_{k+1} = X_k + alpha*P_k
#ifndef TWO_ARG_MATRIX_ADD
      newX = Teuchos::null;
      Xpetra::MatrixMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::TwoMatrixAdd(*P, false, alpha, *X, false, Teuchos::ScalarTraits<Scalar>::one(), newX, mmfancy);
      newX->fillComplete(P0.getDomainMap(), P0.getRangeMap());
      X.swap(newX);
#else
      MatrixMatrix::TwoMatrixAdd(*P, false, alpha, *X, one);
#endif

      if (k == nIts_ - 1)
        break;

      // R_{k+1} = R_k - alpha*A*P_k
#ifndef TWO_ARG_MATRIX_ADD
      newR = Teuchos::null;
      Xpetra::MatrixMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::TwoMatrixAdd(*AP, false, -alpha, *R, false, Teuchos::ScalarTraits<Scalar>::one(), newR, mmfancy);
      newR->fillComplete(P0.getDomainMap(), P0.getRangeMap());
      R.swap(newR);
#else
      MatrixMatrix::TwoMatrixAdd(*AP, false, -alpha, *R, one);
#endif

      // Z_{k+1} = M^{-1} R_{k+1}
      Z = Xpetra::MatrixFactory2<Scalar, LocalOrdinal, GlobalOrdinal, Node>::BuildCopy(R);
      Utilities::MyOldScaleMatrix(*Z, D, true, true, false);

      // beta = (R_{k+1}, Z_{k+1})/(R_k, Z_k)
      newRZ = Utilities::Frobenius(*R, *Z);
      beta = newRZ / oldRZ;

      // P_{k+1} = Z_{k+1} + beta*P_k
#ifndef TWO_ARG_MATRIX_ADD
      newP = Teuchos::null;
      MatrixMatrix::TwoMatrixAdd(*P, false, beta, *Z, false, Teuchos::ScalarTraits<Scalar>::one(), newP, mmfancy);
      newP->fillComplete(P0.getDomainMap(), P0.getRangeMap());
      P.swap(newP);
#else
      MatrixMatrix::TwoMatrixAdd(*Z, false, one, *P, beta);
#endif

      oldRZ = newRZ;
    }

    finalP = X;
  }

} // namespace MueLu

#endif //ifndef MUELU_CGSOLVER_DECL_HPP
