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
#ifndef MUELU_GMRESSOLVER_DEF_HPP
#define MUELU_GMRESSOLVER_DEF_HPP

#include <Teuchos_LAPACK.hpp>

#include <Xpetra_MatrixFactory.hpp>
#include <Xpetra_MatrixMatrix.hpp>
#include <Xpetra_IO.hpp>

#include "MueLu_GMRESSolver_decl.hpp"

#include "MueLu_Constraint.hpp"
#include "MueLu_Monitor.hpp"
#include "MueLu_Utilities.hpp"

namespace MueLu {

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
GMRESSolver<Scalar, LocalOrdinal, GlobalOrdinal, Node>::GMRESSolver(size_t Its)
  : nIts_(Its) {}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void GMRESSolver<Scalar, LocalOrdinal, GlobalOrdinal, Node>::givapp(Scalar* c, Scalar* s, Scalar* v, int k) const {
  for (int i = 0; i < k; i++) {
    SC w1    = c[i] * v[i] - s[i] * v[i + 1];
    SC w2    = s[i] * v[i] + c[i] * v[i + 1];
    v[i]     = w1;
    v[i + 1] = w2;
  }
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void GMRESSolver<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Iterate(const Matrix& Aref, const Constraint& C, const Matrix& P0, RCP<Matrix>& finalP) const {
  PrintMonitor m(*this, "GMRES iterations");

  finalP = MatrixFactory2::BuildCopy(rcpFromRef(P0));
  if (nIts_ == 0)
    return;

  TEUCHOS_TEST_FOR_EXCEPTION(nIts_ > 1, Exceptions::RuntimeError,
                             "For now, solving Hessenberg system works only for a single iteration");

  SC one = Teuchos::ScalarTraits<SC>::one(), zero = Teuchos::ScalarTraits<SC>::zero();

  RCP<const Matrix> A = rcpFromRef(Aref);
  // bool               useTpetra = (A->getRowMap()->lib() == Xpetra::UseTpetra);

  // FIXME: Don't know why, but in the MATLAB code we have D = I. Follow that for now.
#if 0
    ArrayRCP<const SC> D         = Utilities::GetMatrixDiagonal_arcp(*A);
#else
  ArrayRCP<const SC> D(A->getLocalNumRows(), one);
#endif

  Teuchos::FancyOStream& mmfancy = this->GetOStream(Statistics2);

  // Initial P0 would only be used for multiplication
  RCP<Matrix> X = rcp_const_cast<Matrix>(rcpFromRef(P0)), tmpAP, newV;
  std::vector<RCP<Matrix> > V(nIts_ + 1);

  // T is used only for projecting onto
  RCP<CrsMatrix> T_ = CrsMatrixFactory::Build(C.GetPattern());
  T_->fillComplete(P0.getDomainMap(), P0.getRangeMap());
  RCP<Matrix> T = rcp(new CrsMatrixWrap(T_));

  SC rho;
  {
    // R_0 = -D^{-1}*A*X_0
    // V_0 = R_0 / ||R_0||_F
    tmpAP = MatrixMatrix::Multiply(*A, false, *X, false, mmfancy, true /*doFillComplete*/, true /*optimizeStorage*/);
    C.Apply(*tmpAP, *T);

    V[0] = MatrixFactory2::BuildCopy(T);
    Utilities::MyOldScaleMatrix(*V[0], D, true /*doInverse*/, true /*doFillComplete*/, false /*doOptimizeStorage*/);

    rho = sqrt(Utilities::Frobenius(*V[0], *V[0]));

    V[0]->scale(-one / rho);
  }

  std::vector<SC> h((nIts_ + 1) * (nIts_ + 1));
  std::vector<SC> c(nIts_ + 1, 0.0);
  std::vector<SC> s(nIts_ + 1, 0.0);
  std::vector<SC> g(nIts_ + 1, 0.0);
  g[0] = rho;

#define I(i, j) ((i) + (j) * (nIts_ + 1))  // column ordering
  for (size_t i = 0; i < nIts_; i++) {
    // V_{i+1} = D^{-1}*A*V_i
    tmpAP = MatrixMatrix::Multiply(*A, false, *V[i], false, mmfancy, true /*doFillComplete*/, true /*optimizeStorage*/);
    C.Apply(*tmpAP, *T);

    V[i + 1] = MatrixFactory2::BuildCopy(T);
    Utilities::MyOldScaleMatrix(*V[i + 1], D, true /*doInverse*/, true /*doFillComplete*/, false /*doOptimizeStorage*/);

    // Update Hessenberg matrix
    for (size_t j = 0; j <= i; j++) {
      h[I(j, i)] = Utilities::Frobenius(*V[i + 1], *V[j]);

      // V_{i+1} = V_{i+1} - h(j,i+1)*V_j
#ifndef TWO_ARG_MATRIX_ADD
      newV = Teuchos::null;
      MatrixMatrix::TwoMatrixAdd(*V[j], false, -h[I(j, i)], *V[i + 1], false, one, newV, mmfancy);
      newV->fillComplete(V[i + 1]->getDomainMap(), V[i + 1]->getRangeMap());
      V[i + 1].swap(newV);
#else
      // FIXME: this does not work now. Fails with the following exception:
      //   what():  ../../packages/tpetra/core/ext/TpetraExt_MatrixMatrix_def.hpp:408:
      //
      //   Throw number = 1
      //
      //   Throw test that evaluated to true: B.isLocallyIndexed()
      //
      //   TpetraExt::MatrixMatrix::Add(): ERROR, input matrix B must not be locally indexed
      MatrixMatrix::TwoMatrixAdd(*V[j], false, -h[I(j, i)], *V[i + 1], one);
#endif
    }
    h[I(i + 1, i)] = sqrt(Utilities::Frobenius(*V[i + 1], *V[i + 1]));

    // NOTE: potentially we'll need some reorthogonalization code here
    // The matching MATLAB code is
    //    normav  = norm(v.num(k+1).matrix, 'fro');
    //    normav2 = h(k+1,k);
    //    if  (reorth == -1 && normav + .001*normav2 == normav)
    //        for j = 1:k
    //            hr       = v(:,j)'*v(:,k+1);    % hr=v(:,k+1)'*v(:,j);
    //            h(j,k)   = h(j,k)+hr;
    //            v(:,k+1) = v(:,k+1)-hr*v(:,j);
    //        end
    //        h(k+1,k) = norm(v(:,k+1));
    //    end

    // Check for nonsymmetric case
    if (h[I(i + 1, i)] != zero) {
      // Normalize V_i
      V[i + 1]->scale(one / h[I(i + 1, i)]);
    }

    if (i > 0)
      givapp(&c[0], &s[0], &h[I(0, i)], i);  // Due to column ordering &h[...] is a column

    SC nu = sqrt(h[I(i, i)] * h[I(i, i)] + h[I(i + 1, i)] * h[I(i + 1, i)]);
    if (nu != zero) {
      c[i]           = h[I(i, i)] / nu;
      s[i]           = -h[I(i + 1, i)] / nu;
      h[I(i, i)]     = c[i] * h[I(i, i)] - s[i] * h[I(i + 1, i)];
      h[I(i + 1, i)] = zero;

      givapp(&c[i], &s[i], &g[i], 1);
    }
  }

  // Solve Hessenberg system
  //   y = solve(H, \rho e_1)
  std::vector<SC> y(nIts_);
  if (nIts_ == 1) {
    y[0] = g[0] / h[I(0, 0)];
  }
#undef I

  // Compute final
  for (size_t i = 0; i < nIts_; i++) {
#ifndef TWO_ARG_MATRIX_ADD
    newV = Teuchos::null;
    MatrixMatrix::TwoMatrixAdd(*V[i], false, y[i], *finalP, false, one, newV, mmfancy);
    newV->fillComplete(finalP->getDomainMap(), finalP->getRangeMap());
    finalP.swap(newV);
#else
    MatrixMatrix::TwoMatrixAdd(*V[i], false, y[i], *finalP, one);
#endif
  }
}

}  // namespace MueLu

#endif  // ifndef MUELU_GMRESSOLVER_DECL_HPP
