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
#ifndef MUELU_CGSOLVER_DEF_HPP
#define MUELU_CGSOLVER_DEF_HPP

#include <Xpetra_MatrixFactory.hpp>

#include "MueLu_CGSolver_decl.hpp"

#include "MueLu_Constraint.hpp"
#include "MueLu_Utilities.hpp"

namespace MueLu {

  using Teuchos::rcp_const_cast;

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  CGSolver<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::CGSolver(size_t Its)
  : nIts_(Its)
  { }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  void CGSolver<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::Iterate(const Matrix& Aref, const Constraint& C, const Matrix& P0, const MultiVector& B, RCP<Matrix>& finalP) const {
    // Note: this function matrix notations follow Saad's "Iterative methods", ed. 2, pg. 246
    // So, X is the unknown prolongator, P's are conjugate directions, Z's are preconditioned P's
    RCP<const Matrix> A = rcpFromRef(Aref);

    RCP<Matrix> X, P, R, Z, AP;
    RCP<Matrix> newX, tmpAP;
    SC oldRZ, newRZ, alpha, beta, app;

    bool useTpetra = (A->getRowMap()->lib() == Xpetra::UseTpetra);

    // T is used only for projecting onto
    RCP<CrsMatrix> T_ = CrsMatrixFactory::Build(C.GetPattern());
    T_->fillComplete(P0.getDomainMap(), P0.getRangeMap());
    RCP<Matrix>    T = rcp(new CrsMatrixWrap(T_));

    SC one = Teuchos::ScalarTraits<SC>::one();

    Teuchos::ArrayRCP<const SC> D = Utils::GetMatrixDiagonal(*A);

    // Initial P0 would only be used for multiplication
    X = rcp_const_cast<Matrix>(rcpFromRef(P0));

    tmpAP = Utils::Multiply(*A, false, *X, false, true, false);
    C.Apply(*tmpAP, *T);

    // R_0 = -A*X_0
    R = MatrixFactory::BuildCopy(T);
#ifdef HAVE_MUELU_TPETRA
    if (useTpetra)
      Utils::Op2NonConstTpetraCrs(R)->resumeFill();
#endif
    R->scale(-one);
    if (useTpetra)
      R->fillComplete(R->getDomainMap(), R->getRangeMap());

    // Z_0 = M^{-1}R_0
    Z = MatrixFactory::BuildCopy(R);
    Utils::MyOldScaleMatrix(Z, D, true, true, false);

    // P_0 = Z_0
    P = MatrixFactory::BuildCopy(Z);

    oldRZ = Frobenius(*R, *Z);

    for (size_t k = 0; k < nIts_; k++) {
      // AP = constrain(A*P)
      tmpAP = Utils::Multiply(*A, false, *P, false, true, false);
      C.Apply(*tmpAP, *T);
      AP = T;

      app = Frobenius(*AP, *P);
      if (Teuchos::ScalarTraits<SC>::magnitude(app) < Teuchos::ScalarTraits<SC>::sfmin()) {
        // It happens, for instance, if P = 0
        // For example, if we use TentativePFactory for both nonzero pattern and initial guess
        // I think it might also happen because of numerical breakdown, but we don't test for that yet
        if (k == 0)
          X = MatrixFactory::BuildCopy(rcpFromRef(P0));
        break;
      }

      // alpha = (R_k, Z_k)/(A*P_k, P_k)
      alpha = oldRZ / app;
      std::cout << "emin: alpha = " << alpha << std::endl;

      // X_{k+1} = X_k + alpha*P_k
#if 0
      Utils2::TwoMatrixAdd(P, false, alpha, X, one);
#else
      newX = Teuchos::null;
      Utils2::TwoMatrixAdd(P, false, alpha, X, false, Teuchos::ScalarTraits<Scalar>::one(), newX);
      newX->fillComplete(P0.getDomainMap(), P0.getRangeMap());
      X.swap(newX);
#endif

      if (k == nIts_ - 1)
        break;

      // R_{k+1} = R_k - alpha*A*P_k
      Utils2::TwoMatrixAdd(AP, false, -alpha, R, one);

      // Z_{k+1} = M^{-1} R_{k+1}
      Z = MatrixFactory::BuildCopy(R);
      Utils::MyOldScaleMatrix(Z, D, true, true, false);

      // beta = (R_{k+1}, Z_{k+1})/(R_k, Z_k)
      newRZ = Frobenius(*R, *Z);
      beta = newRZ / oldRZ;

      // P_{k+1} = Z_{k+1} + beta*P_k
      Utils2::TwoMatrixAdd(Z, false, one, P, beta);

      oldRZ = newRZ;
    }

    finalP = X;
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  Scalar CGSolver<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::Frobenius(const Matrix& Aref, const Matrix& Bref) const {
    RCP<const Matrix> A = rcpFromRef(Aref), B = rcpFromRef(Bref);

    size_t numRows = A->getNodeNumRows();

    // In the future, one might need to restrict this test for only row maps
    // For instance, if matrix B = M*A then they would have different colmaps
    // See comments in the loop how to update the algorithm
    TEUCHOS_TEST_FOR_EXCEPTION(!A->getRowMap()->isSameAs(*(B->getRowMap())), Exceptions::Incompatible, "MueLu::CGSolver::Frobenius: row maps are incompatible");
    TEUCHOS_TEST_FOR_EXCEPTION(!A->getColMap()->isSameAs(*(B->getColMap())), Exceptions::Incompatible, "MueLu::CGSolver::Frobenius: col maps are incompatible");

    SC f = Teuchos::ScalarTraits<SC>::zero();
    for (size_t i = 0; i < numRows; i++) {
      Teuchos::ArrayView<const LO> indA, indB;
      Teuchos::ArrayView<const SC> valA, valB;

      // If A and B have different colmaps, we need to replace
      // indA and indB by their corresponding global indices
      // It can be done, for instance, using getGlobalElement() function.
      // We would also probably need to sort those GIDs.
      A->getLocalRowView(i, indA, valA);
      B->getLocalRowView(i, indB, valB);

      size_t nnzA = indA.size()/*, nnzB = indB.size()*/;

      // We assume that indA and indB are sorted in increasing order
      for (size_t j0 = 0, j1 = 0; j0 < nnzA;) {
        if (indA[j0] < indB[j1])
          j0++;
        else if (indA[j0] > indB[j1])
          j1++;
        else {
          f += valA[j0]*valB[j1];
          j0++;
          j1++;
        }
      }
    }

    SC gf;
    sumAll(A->getRowMap()->getComm(), f, gf);

    return gf;
  }

} // namespace MueLu

#endif //ifndef MUELU_CGSOLVER_DECL_HPP
