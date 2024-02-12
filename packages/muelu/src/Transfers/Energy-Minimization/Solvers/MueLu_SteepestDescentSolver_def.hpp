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
#ifndef MUELU_STEEPESTDESCENTSOLVER_DEF_HPP
#define MUELU_STEEPESTDESCENTSOLVER_DEF_HPP

#include <Xpetra_CrsMatrixFactory.hpp>
#include <Xpetra_CrsMatrixWrap.hpp>
#include <Xpetra_MatrixMatrix.hpp>

#include "MueLu_Constraint.hpp"
#include "MueLu_Monitor.hpp"
#include "MueLu_Utilities.hpp"

#include "MueLu_SteepestDescentSolver_decl.hpp"

namespace MueLu {

using Teuchos::rcp_const_cast;

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
SteepestDescentSolver<Scalar, LocalOrdinal, GlobalOrdinal, Node>::SteepestDescentSolver(size_t Its, SC StepLength)
  : nIts_(Its)
  , stepLength_(StepLength) {}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void SteepestDescentSolver<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Iterate(const Matrix& Aref, const Constraint& C, const Matrix& P0, RCP<Matrix>& P) const {
  PrintMonitor m(*this, "SD iterations");

  RCP<const Matrix> A = rcpFromRef(Aref);
  RCP<Matrix> AP, G;

  Teuchos::FancyOStream& mmfancy = this->GetOStream(Statistics2);

  Teuchos::ArrayRCP<const SC> D = Utilities::GetMatrixDiagonal_arcp(*A);

  RCP<CrsMatrix> Ptmp_ = CrsMatrixFactory::Build(C.GetPattern());
  Ptmp_->fillComplete(P0.getDomainMap(), P0.getRangeMap());
  RCP<Matrix> Ptmp = rcp(new CrsMatrixWrap(Ptmp_));

  // Initial P0 would only be used for multiplication
  P = rcp_const_cast<Matrix>(rcpFromRef(P0));

  for (size_t k = 0; k < nIts_; k++) {
    AP = Xpetra::MatrixMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Multiply(*A, false, *P, false, mmfancy, true, true);
#if 0
      // gradient = -2 A^T * A * P
      SC stepLength = 2*stepLength_;
      G = Xpetra::MatrixMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Multiply(*A, true, *AP, false, true, true);
      C.Apply(*G, *Ptmp);
#else
    // gradient = - A * P
    SC stepLength = stepLength_;
    Utilities::MyOldScaleMatrix(*AP, D, true, true, false);
    C.Apply(*AP, *Ptmp);
#endif

    RCP<Matrix> newP;
    Xpetra::MatrixMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::TwoMatrixAdd(*Ptmp, false, -stepLength, *P, false, Teuchos::ScalarTraits<Scalar>::one(), newP, mmfancy);
    newP->fillComplete(P->getDomainMap(), P->getRangeMap());
    P = newP;
  }
}
}  // namespace MueLu

#endif  // ifndef MUELU_STEEPESTDESCENTSOLVER_DECL_HPP
