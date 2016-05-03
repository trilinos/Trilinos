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
#ifndef MUELU_SCHURCOMPLEMENTFACTORY_DEF_HPP_
#define MUELU_SCHURCOMPLEMENTFACTORY_DEF_HPP_

#ifdef HAVE_MUELU_EXPERIMENTAL

#include <Xpetra_BlockedCrsMatrix.hpp>
#include <Xpetra_MultiVectorFactory.hpp>
#include <Xpetra_VectorFactory.hpp>
#include <Xpetra_MatrixFactory.hpp>
#include <Xpetra_Matrix.hpp>
#include <Xpetra_MatrixMatrix.hpp>
#include <Xpetra_CrsMatrixWrap.hpp>
#include <Xpetra_BlockedCrsMatrix.hpp>
#include <Xpetra_CrsMatrix.hpp>
#include "MueLu_Level.hpp"
#include "MueLu_Monitor.hpp"
#include "MueLu_Utilities.hpp"
//#include "MueLu_HierarchyHelpers.hpp"

#include "MueLu_SchurComplementFactory.hpp"

namespace MueLu {

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  RCP<const ParameterList> SchurComplementFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::GetValidParameterList() const {
    RCP<ParameterList> validParamList = rcp(new ParameterList());

    SC one = Teuchos::ScalarTraits<SC>::one();

    validParamList->set<RCP<const FactoryBase> >("A", NoFactory::getRCP()/*null*/, "Generating factory of the matrix A used for building Schur complement\n"
                                                                                   "(must be a 2x2 blocked operator)");
    validParamList->set<SC>                     ("omega",                     one, "Scaling parameter in S = A(1,1) - 1/omega A(1,0) diag{A(0,0)}^{-1} A(0,1)");
    validParamList->set<bool>                   ("lumping",                 false, "Use lumping to construct diag(A(0,0)), i.e. use row sum of the abs values on the diagonal "
                                                                                   "as approximation of A00 (and A00^{-1})");
    validParamList->set<bool>                   ("fixing",                  false, "Fix diagonal by replacing small entries with 1.0");

    return validParamList;
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void SchurComplementFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::DeclareInput(Level& currentLevel) const {
    Input(currentLevel, "A");
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void SchurComplementFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Build(Level& currentLevel) const {
    FactoryMonitor m(*this, "Build", currentLevel);

    typedef Teuchos::ScalarTraits<SC> STS;
    SC zero = STS::zero(), one = STS::one();

    RCP<Matrix>            A = Get<RCP<Matrix> >(currentLevel, "A");
    RCP<BlockedCrsMatrix> bA = rcp_dynamic_cast<BlockedCrsMatrix>(A);
    TEUCHOS_TEST_FOR_EXCEPTION(bA.is_null(), Exceptions::BadCast,
                               "MueLu::SchurComplementFactory::Build: input matrix A is not of type BlockedCrsMatrix!");

    RCP<Matrix> A00 = bA->getMatrix(0,0);
    RCP<Matrix> A01 = bA->getMatrix(0,1);
    RCP<Matrix> A10 = bA->getMatrix(1,0);
    RCP<Matrix> A11 = bA->getMatrix(1,1);

    const ParameterList& pL = GetParameterList();
    SC omega = pL.get<Scalar>("omega");

    // Copy the value of A01 so we can do the left scale.
    RCP<Matrix> T = MatrixFactory2::BuildCopy(A01);

    bool lumping = pL.get<bool>("lumping");
    bool fixing  = pL.get<bool>("fixing");
    ArrayRCP<SC> D;
    if (!lumping) {
      D = Utilities::GetMatrixDiagonal(*A00);

      if (fixing) {
        for (size_t k = 0; k < as<size_t>(D.size()); k++)
          if (STS::magnitude(D[k]) < 1e-4)
            D[k] = STS::one();
      }

    } else {
      D = Utilities::GetLumpedMatrixDiagonal(*A00);
    }
    // Update D to use omega
    // As D is going to be used as D^{-1}, we must use -omega instead of -one/omega
    for (size_t k = 0; k < as<size_t>(D.size()); k++)
      D[k] *= -omega;
    Utilities::MyOldScaleMatrix(*T, D, true/*doInverse*/, true/*doFillComplete*/, false/*doOptimizeStorage*/);

    RCP<Matrix> S = Xpetra::MatrixMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Multiply(*A10, false, *T, false, GetOStream(Statistics2));

    if (!A11.is_null()) {
      T = Teuchos::null;
      Xpetra::MatrixMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::TwoMatrixAdd(*A11, false, one, *S, false, one, T, GetOStream(Statistics2));
      T->fillComplete();
      S.swap(T);
    }

    // NOTE: "A" generated by this factory is actually the Schur complement
    // matrix, but it is required as all smoothers expect "A"
    Set(currentLevel, "A", S);
  }

} // namespace MueLu

#endif /* HAVE_MUELU_EXPERIMENTAL */
#endif /* MUELU_SCHURCOMPLEMENTFACTORY_DEF_HPP_ */
