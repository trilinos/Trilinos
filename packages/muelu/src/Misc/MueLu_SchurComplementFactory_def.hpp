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

#include <Xpetra_BlockedCrsMatrix.hpp>
#include <Xpetra_MultiVectorFactory.hpp>
#include <Xpetra_VectorFactory.hpp>
#include <Xpetra_MatrixFactory.hpp>
#include <Xpetra_Matrix.hpp>
#include <Xpetra_MatrixMatrix.hpp>
#include <Xpetra_TripleMatrixMultiply.hpp>
#include <Xpetra_CrsMatrixWrap.hpp>
#include <Xpetra_CrsMatrix.hpp>

#include "MueLu_Level.hpp"
#include "MueLu_Monitor.hpp"
#include "MueLu_SchurComplementFactory_decl.hpp"

namespace MueLu {

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  RCP<const ParameterList> SchurComplementFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::GetValidParameterList() const {
    RCP<ParameterList> validParamList = rcp(new ParameterList());

    const SC one = Teuchos::ScalarTraits<SC>::one();

    validParamList->set<RCP<const FactoryBase> >("A"    , NoFactory::getRCP(), "Generating factory of the matrix A used for building Schur complement (must be a 2x2 blocked operator)");
    validParamList->set<RCP<const FactoryBase> >("Ainv" , Teuchos::null,       "Generating factory of the inverse matrix used in the Schur complement");

    validParamList->set<SC>                     ("omega", one, "Scaling parameter in S = A(1,1) - 1/omega A(1,0) Ainv A(0,1)");

    return validParamList;
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void SchurComplementFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::DeclareInput(Level& currentLevel) const {
    Input(currentLevel, "A");
    
    // Get default or user-given inverse approximation factory
    RCP<const FactoryBase> AinvFact = GetFactory("Ainv");
    currentLevel.DeclareInput("Ainv", AinvFact.get(), this);
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void SchurComplementFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Build(Level& currentLevel) const {
    FactoryMonitor m(*this, "Build", currentLevel);

    RCP<Matrix> A = Get<RCP<Matrix> >(currentLevel, "A");
    RCP<BlockedCrsMatrix> bA = rcp_dynamic_cast<BlockedCrsMatrix>(A);

    TEUCHOS_TEST_FOR_EXCEPTION(bA.is_null(), Exceptions::BadCast,
                               "MueLu::SchurComplementFactory::Build: input matrix A is not of type BlockedCrsMatrix!");
    TEUCHOS_TEST_FOR_EXCEPTION(bA->Rows() != 2 || bA->Cols() != 2, Exceptions::RuntimeError,
                               "MueLu::SchurComplementFactory::Build: input matrix A is a " << bA->Rows() << "x" << bA->Cols() << " block matrix. We expect a 2x2 blocked operator.");

    // Calculate Schur Complement
    RCP<Matrix> Ainv = currentLevel.Get<RCP<Matrix> >("Ainv", this->GetFactory("Ainv").get());
    RCP<Matrix> S = ComputeSchurComplement(bA, Ainv);

    GetOStream(Statistics1) << "S has " << S->getGlobalNumRows() << "x" << S->getGlobalNumCols() << " rows and columns." << std::endl;

    // NOTE: "A" generated by this factory is actually the Schur complement
    // matrix, but it is required as all smoothers expect "A"
    Set(currentLevel, "A", S);
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  RCP<Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>>
  SchurComplementFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::ComputeSchurComplement(RCP<BlockedCrsMatrix>& bA, RCP<Matrix>& Ainv) const {

    using STS = Teuchos::ScalarTraits<SC>;
    const SC zero = STS::zero(), one = STS::one();

    RCP<Matrix> A01 = bA->getMatrix(0,1);
    RCP<Matrix> A10 = bA->getMatrix(1,0);
    RCP<Matrix> A11 = bA->getMatrix(1,1);

    RCP<BlockedCrsMatrix> bA01 = Teuchos::rcp_dynamic_cast<BlockedCrsMatrix>(A01);
    const bool isBlocked = (bA01 == Teuchos::null ? false : true);

    const ParameterList& pL = GetParameterList();
    const SC omega = pL.get<Scalar>("omega");

    TEUCHOS_TEST_FOR_EXCEPTION(omega == zero, Exceptions::RuntimeError,
                               "MueLu::SchurComplementFactory::Build: Scaling parameter omega must not be zero to avoid division by zero.");

    RCP<Matrix> S = Teuchos::null; // Schur complement
    RCP<Matrix> D = Teuchos::null; // temporary result for A10*Ainv*A01

    // only if the off-diagonal blocks A10 and A01 are non-zero we have to do the MM multiplication
    if(A01.is_null() == false && A10.is_null() == false) {
      // scale with -1/omega
      Ainv->scale(Teuchos::as<Scalar>(-one/omega));

      // build Schur complement operator
      if (!isBlocked) {
        RCP<ParameterList> myparams = rcp(new ParameterList);
        myparams->set("compute global constants", true);

        // -1/omega*Ainv*A01
        TEUCHOS_TEST_FOR_EXCEPTION(A01->getRangeMap()->isSameAs(*(Ainv->getDomainMap())) == false, Exceptions::RuntimeError,
                                   "MueLu::SchurComplementFactory::Build: RangeMap of A01 and domain map of Ainv are not the same.");
        RCP<Matrix> C = MatrixMatrix::Multiply(*Ainv, false, *A01, false, GetOStream(Statistics2), true, true, std::string("SchurComplementFactory"), myparams);

        // -1/omega*A10*Ainv*A01
        TEUCHOS_TEST_FOR_EXCEPTION(A01->getRangeMap()->isSameAs(*(A10->getDomainMap())) == false, Exceptions::RuntimeError,
                                   "MueLu::SchurComplementFactory::Build: RangeMap of A10 and domain map A01 are not the same.");
        D = MatrixMatrix::Multiply(*A10, false, *C, false, GetOStream(Statistics2), true, true, std::string("SchurComplementFactory"), myparams);
      }
      else {
        // nested blocking
        auto bA10  = Teuchos::rcp_dynamic_cast<BlockedCrsMatrix>(A10);
        auto bAinv = Teuchos::rcp_dynamic_cast<BlockedCrsMatrix>(Ainv);
        TEUCHOS_TEST_FOR_EXCEPTION(bAinv == Teuchos::null, Exceptions::RuntimeError,
                                   "MueLu::SchurComplementFactory::Build: Casting Ainv to BlockedCrsMatrix not possible.");

        // -1/omega*bAinv*bA01
        TEUCHOS_TEST_FOR_EXCEPTION(bA01->Rows() != bAinv->Cols(), Exceptions::RuntimeError,
                                   "MueLu::SchurComplementFactory::Build: Block rows and cols of bA01 and bAinv are not compatible.");
        RCP<BlockedCrsMatrix> C = MatrixMatrix::TwoMatrixMultiplyBlock(*bAinv, false, *bA01, false, GetOStream(Statistics2));

        // -1/omega*A10*Ainv*A01
        TEUCHOS_TEST_FOR_EXCEPTION(bA10->Rows() != bA01->Cols(), Exceptions::RuntimeError,
                                   "MueLu::SchurComplementFactory::Build: Block rows and cols of bA10 and bA01 are not compatible.");
        D = MatrixMatrix::TwoMatrixMultiplyBlock(*bA10, false, *C, false, GetOStream(Statistics2));
      }
      if (!A11.is_null()) {
        MatrixMatrix::TwoMatrixAdd(*A11, false, one, *D, false, one, S, GetOStream(Statistics2));
        S->fillComplete();

        TEUCHOS_TEST_FOR_EXCEPTION(A11->getRangeMap()->isSameAs(*(S->getRangeMap())) == false, Exceptions::RuntimeError,
                                   "MueLu::SchurComplementFactory::Build: RangeMap of A11 and S are not the same.");
        TEUCHOS_TEST_FOR_EXCEPTION(A11->getDomainMap()->isSameAs(*(S->getDomainMap())) == false, Exceptions::RuntimeError,
                                   "MueLu::SchurComplementFactory::Build: DomainMap of A11 and S are not the same.");
      }
      else {
        S = MatrixFactory::BuildCopy(D);
      }
    }
    else {
      if (!A11.is_null()) {
        S = MatrixFactory::BuildCopy(A11);
      } else {
        S = MatrixFactory::Build(A11->getRowMap(), 10 /*A11->getLocalMaxNumRowEntries()*/);
        S->fillComplete(A11->getDomainMap(),A11->getRangeMap());
      }
    }

    // Check whether Schur complement operator is a 1x1 block matrix.
    // If so, unwrap it and return the CrsMatrix based Matrix object
    // We need this, as single-block smoothers expect it this way.
    // In case of Thyra GIDs we obtain a Schur complement operator in Thyra GIDs
    // This may make some special handling in feeding the SchurComplement solver Apply routine
    // necessary!
    if (isBlocked) {
      RCP<BlockedCrsMatrix> bS = Teuchos::rcp_dynamic_cast<BlockedCrsMatrix>(S);

      if (bS != Teuchos::null && bS->Rows() == 1 && bS->Cols() == 1) {
        RCP<Matrix> temp = bS->getCrsMatrix();
        S.swap(temp);
      }
    }

    return S;
  }

} // namespace MueLu

#endif /* MUELU_SCHURCOMPLEMENTFACTORY_DEF_HPP_ */
