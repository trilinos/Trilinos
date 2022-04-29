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
#ifndef MUELU_INVERSEAPPROXIMATIONFACTORY_DEF_HPP_
#define MUELU_INVERSEAPPROXIMATIONFACTORY_DEF_HPP_

#include <Xpetra_BlockedCrsMatrix.hpp>
#include <Xpetra_MultiVectorFactory.hpp>
#include <Xpetra_VectorFactory.hpp>
#include <Xpetra_MatrixFactory.hpp>
#include <Xpetra_Matrix.hpp>
#include <Xpetra_MatrixMatrix.hpp>
#include <Xpetra_TripleMatrixMultiply.hpp>
#include <Xpetra_CrsMatrixWrap.hpp>
#include <Xpetra_BlockedCrsMatrix.hpp>
#include <Xpetra_CrsMatrix.hpp>

#include "MueLu_Level.hpp"
#include "MueLu_Monitor.hpp"
#include "MueLu_Utilities.hpp"
#include "MueLu_InverseApproximationFactory_decl.hpp"

namespace MueLu {

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  RCP<const ParameterList> InverseApproximationFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::GetValidParameterList() const {
    RCP<ParameterList> validParamList = rcp(new ParameterList());

    validParamList->set<RCP<const FactoryBase> >("A", NoFactory::getRCP(), "Matrix to build the approximate inverse on.\n");

    validParamList->set<std::string>            ("inverse: approximation type",  "diagonal", "Method used to approximate the inverse.");
    validParamList->set<bool>                   ("inverse: fixing",              false     , "Fix diagonal by replacing small entries with 1.0");

    return validParamList;
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void InverseApproximationFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::DeclareInput(Level& currentLevel) const {
    Input(currentLevel, "A");
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void InverseApproximationFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Build(Level& currentLevel) const {
    FactoryMonitor m(*this, "Build", currentLevel);

    typedef Teuchos::ScalarTraits<SC> STS;
    SC one = STS::one();

    const ParameterList& pL = GetParameterList();

    // check which approximation type to use
    bool fixing = pL.get<bool>("inverse: fixing");
    std::string method = pL.get<std::string>("inverse: approximation type");
    TEUCHOS_TEST_FOR_EXCEPTION(method != "diagonal" && method != "lumping", Exceptions::RuntimeError,
                               "MueLu::InverseApproximationFactory::Build: Approximation type can be 'diagonal' or 'lumping'.");

    RCP<Matrix> A = Get<RCP<Matrix> >(currentLevel, "A");
    RCP<BlockedCrsMatrix> bA = Teuchos::rcp_dynamic_cast<BlockedCrsMatrix>(A);
    bool isBlocked = (bA == Teuchos::null ? false : true);

    // if blocked operator is used, defaults to A(0,0)
    if(isBlocked) A = bA->getMatrix(0,0);

    RCP<Vector> Ainv = Teuchos::null;
    if(method=="diagonal") {
      auto diag = VectorFactory::Build(A->getRangeMap(), true);
      A->getLocalDiagCopy(*diag);
      Ainv = (!fixing ? Utilities::GetInverse(diag) : Utilities::GetInverse(diag, 1e-4, one));
    }
    else if(method=="lumping") {
      auto diag = Utilities::GetLumpedMatrixDiagonal(*A);
      Ainv = (!fixing ? Utilities::GetInverse(diag) : Utilities::GetInverse(diag, 1e-4, one));
    }
    
    GetOStream(Statistics1) << "Approximate inverse calculated by: " << method << "." << std::endl;

    Set(currentLevel, "Ainv", Ainv);
  }

} // namespace MueLu

#endif /* MUELU_INVERSEAPPROXIMATIONFACTORY_DEF_HPP_ */
