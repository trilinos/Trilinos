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
    // Note: this function matrix notations follow Saad's "Iterative methods", ed. 1, pg. 246
    // So, X is the unknown prolongator, P's are conjugate directions, Z's are preconditioned P's
    RCP<const Matrix> A = rcpFromRef(Aref);

    RCP<Matrix> X, P, R, Z, AP;
    RCP<Matrix> newX, newR, newP, tmpAP;
    SC oldRZ, newRZ, alpha, beta, app;

    // T is used only for projecting onto
    RCP<CrsMatrix> T_ = CrsMatrixFactory::Build(C.GetPattern());
    T_->fillComplete(P0.getDomainMap(), P0.getRangeMap());
    RCP<Matrix>    T = rcp(new CrsMatrixWrap(T_));

    Teuchos::ArrayRCP<const SC> D = Utils::GetMatrixDiagonal(*A);

    // Initial P0 would only be used for multiplication
    X = rcp_const_cast<Matrix>(rcpFromRef(P0));

    tmpAP = Utils::Multiply(*A, false, *X, false, true, false);
    C.Apply(*tmpAP, *T);

    R = MatrixFactory::BuildCopy(T);

    bool useTpetra = (A->getRowMap()->lib() == Xpetra::UseTpetra);

#ifdef HAVE_MUELU_TPETRA
    if (useTpetra)
      Utils::Op2NonConstTpetraCrs(R)->resumeFill();
#endif
    R->scale(-Teuchos::ScalarTraits<SC>::one());
    if (useTpetra)
      R->fillComplete(R->getDomainMap(), R->getRangeMap());

    Z = MatrixFactory::BuildCopy(R);
    Utils::MyOldScaleMatrix(Z, D, true, true, false);

    P = MatrixFactory::BuildCopy(Z);

    oldRZ = Frobenius(*R, *Z);

    for (size_t k = 0; k < nIts_; k++) {
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

      alpha = oldRZ / app;

      newX = Teuchos::null;
      Utils2::TwoMatrixAdd(P, false, alpha, X, false, Teuchos::ScalarTraits<Scalar>::one(), newX);
      newX->fillComplete(P0.getDomainMap(), P0.getRangeMap());
      X.swap(newX);

      if (k == nIts_ - 1)
        break;

      newR = Teuchos::null;
      Utils2::TwoMatrixAdd(AP, false, -alpha, R, false, Teuchos::ScalarTraits<Scalar>::one(), newR);
      newR->fillComplete(P0.getDomainMap(), P0.getRangeMap());
      R.swap(newR);

      Z = MatrixFactory::BuildCopy(R);
      Utils::MyOldScaleMatrix(Z, D, true, true, false);

      newRZ = Frobenius(*R, *Z);
      beta = newRZ / oldRZ;

      newP = Teuchos::null;
      Utils2::TwoMatrixAdd(P, false, beta, Z, false, Teuchos::ScalarTraits<Scalar>::one(), newP);
      newP->fillComplete(P0.getDomainMap(), P0.getRangeMap());
      P.swap(newP);

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
    TEUCHOS_TEST_FOR_EXCEPTION(!A->getRowMap()->isSameAs(*(B->getRowMap())) ||
                               !A->getColMap()->isSameAs(*(B->getColMap())),
                               Exceptions::Incompatible, "MueLu::CGSolver: matrices' maps are incompatible");

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
      for (size_t j0 = 0, j1 = 0; j0 < nnzA; j0++) {
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
