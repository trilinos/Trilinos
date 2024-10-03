// @HEADER
// *****************************************************************************
//        MueLu: A package for multigrid based preconditioning
//
// Copyright 2012 NTESS and the MueLu contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef MUELU_INVERSEAPPROXIMATIONFACTORY_DEF_HPP_
#define MUELU_INVERSEAPPROXIMATIONFACTORY_DEF_HPP_

#include <Xpetra_BlockedCrsMatrix.hpp>
#include <Xpetra_CrsGraph.hpp>
#include <Xpetra_CrsGraphFactory.hpp>
#include <Xpetra_CrsMatrixWrap.hpp>
#include <Xpetra_CrsMatrix.hpp>
#include <Xpetra_MultiVectorFactory.hpp>
#include <Xpetra_VectorFactory.hpp>
#include <Xpetra_MatrixFactory.hpp>
#include <Xpetra_Matrix.hpp>
#include <Xpetra_MatrixMatrix.hpp>
#include <Xpetra_TripleMatrixMultiply.hpp>

#include <Teuchos_SerialDenseVector.hpp>
#include <Teuchos_SerialDenseMatrix.hpp>
#include <Teuchos_SerialQRDenseSolver.hpp>

#include "MueLu_Level.hpp"
#include "MueLu_Monitor.hpp"
#include "MueLu_Utilities.hpp"
#include "MueLu_InverseApproximationFactory_decl.hpp"

namespace MueLu {

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
RCP<const ParameterList> InverseApproximationFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::GetValidParameterList() const {
  RCP<ParameterList> validParamList = rcp(new ParameterList());
  using Magnitude                   = typename Teuchos::ScalarTraits<Scalar>::magnitudeType;

  validParamList->set<RCP<const FactoryBase>>("A", NoFactory::getRCP(), "Matrix to build the approximate inverse on.\n");

  validParamList->set<std::string>("inverse: approximation type", "diagonal", "Method used to approximate the inverse.");
  validParamList->set<Magnitude>("inverse: drop tolerance", 0.0, "Values below this threshold  are dropped from the matrix (or fixed if diagonal fixing is active).");
  validParamList->set<bool>("inverse: fixing", false, "Keep diagonal and fix small entries with 1.0");

  return validParamList;
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void InverseApproximationFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::DeclareInput(Level& currentLevel) const {
  Input(currentLevel, "A");
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void InverseApproximationFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Build(Level& currentLevel) const {
  FactoryMonitor m(*this, "Build", currentLevel);

  using STS       = Teuchos::ScalarTraits<SC>;
  const SC one    = STS::one();
  using Magnitude = typename Teuchos::ScalarTraits<Scalar>::magnitudeType;

  const ParameterList& pL = GetParameterList();
  const bool fixing       = pL.get<bool>("inverse: fixing");

  // check which approximation type to use
  const std::string method = pL.get<std::string>("inverse: approximation type");
  TEUCHOS_TEST_FOR_EXCEPTION(method != "diagonal" && method != "lumping" && method != "sparseapproxinverse", Exceptions::RuntimeError,
                             "MueLu::InverseApproximationFactory::Build: Approximation type can be 'diagonal' or 'lumping' or "
                             "'sparseapproxinverse'.");

  RCP<Matrix> A            = Get<RCP<Matrix>>(currentLevel, "A");
  RCP<BlockedCrsMatrix> bA = Teuchos::rcp_dynamic_cast<BlockedCrsMatrix>(A);
  const bool isBlocked     = (bA == Teuchos::null ? false : true);

  // if blocked operator is used, defaults to A(0,0)
  if (isBlocked) A = bA->getMatrix(0, 0);

  const Magnitude tol = pL.get<Magnitude>("inverse: drop tolerance");
  RCP<Matrix> Ainv    = Teuchos::null;

  if (method == "diagonal") {
    const auto diag = VectorFactory::Build(A->getRangeMap(), true);
    A->getLocalDiagCopy(*diag);
    const RCP<const Vector> D = (!fixing ? Utilities::GetInverse(diag) : Utilities::GetInverse(diag, tol, one));
    Ainv                      = MatrixFactory::Build(D);
  } else if (method == "lumping") {
    const auto diag           = Utilities::GetLumpedMatrixDiagonal(*A);
    const RCP<const Vector> D = (!fixing ? Utilities::GetInverse(diag) : Utilities::GetInverse(diag, tol, one));
    Ainv                      = MatrixFactory::Build(D);
  } else if (method == "sparseapproxinverse") {
    RCP<CrsGraph> sparsityPattern = Utilities::GetThresholdedGraph(A, tol, A->getGlobalMaxNumRowEntries());
    GetOStream(Statistics1) << "NNZ Graph(A): " << A->getCrsGraph()->getGlobalNumEntries() << " , NNZ Tresholded Graph(A): " << sparsityPattern->getGlobalNumEntries() << std::endl;
    RCP<Matrix> pAinv = GetSparseInverse(A, sparsityPattern);
    Ainv              = Utilities::GetThresholdedMatrix(pAinv, tol, fixing, pAinv->getGlobalMaxNumRowEntries());
    GetOStream(Statistics1) << "NNZ Ainv: " << pAinv->getGlobalNumEntries() << ", NNZ Tresholded Ainv (parameter: " << tol << "): " << Ainv->getGlobalNumEntries() << std::endl;
  }

  GetOStream(Statistics1) << "Approximate inverse calculated by: " << method << "." << std::endl;
  GetOStream(Statistics1) << "Ainv has " << Ainv->getGlobalNumRows() << "x" << Ainv->getGlobalNumCols() << " rows and columns." << std::endl;

  Set(currentLevel, "Ainv", Ainv);
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
RCP<Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>>
InverseApproximationFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::GetSparseInverse(const RCP<Matrix>& Aorg, const RCP<const CrsGraph>& sparsityPattern) const {
  // construct the inverse matrix with the given sparsity pattern
  RCP<Matrix> Ainv = MatrixFactory::Build(sparsityPattern);
  Ainv->resumeFill();

  // gather missing rows from other procs to generate an overlapping map
  RCP<Import> rowImport = ImportFactory::Build(sparsityPattern->getRowMap(), sparsityPattern->getColMap());
  RCP<Matrix> A         = MatrixFactory::Build(Aorg, *rowImport);

  // loop over all rows of the inverse sparsity pattern (this can be done in parallel)
  for (size_t k = 0; k < sparsityPattern->getLocalNumRows(); k++) {
    // 1. get column indices Ik of local row k
    ArrayView<const LO> Ik;
    sparsityPattern->getLocalRowView(k, Ik);

    // 2. get all local A(Ik,:) rows
    Array<ArrayView<const LO>> J(Ik.size());
    Array<ArrayView<const SC>> Ak(Ik.size());
    Array<LO> Jk;
    for (LO i = 0; i < Ik.size(); i++) {
      A->getLocalRowView(Ik[i], J[i], Ak[i]);
      for (LO j = 0; j < J[i].size(); j++)
        Jk.append(J[i][j]);
    }
    // set of unique column indices Jk
    std::sort(Jk.begin(), Jk.end());
    Jk.erase(std::unique(Jk.begin(), Jk.end()), Jk.end());
    // create map
    std::map<LO, LO> G;
    for (LO i = 0; i < Jk.size(); i++) G.insert(std::pair<LO, LO>(Jk[i], i));

    // 3. merge rows together
    Teuchos::SerialDenseMatrix<LO, SC> localA(Jk.size(), Ik.size(), true);
    for (LO i = 0; i < Ik.size(); i++) {
      for (LO j = 0; j < J[i].size(); j++) {
        localA(G.at(J[i][j]), i) = Ak[i][j];
      }
    }

    // 4. get direction-vector
    // diagonal needs an entry!
    Teuchos::SerialDenseVector<LO, SC> ek(Jk.size(), true);
    ek[std::find(Jk.begin(), Jk.end(), k) - Jk.begin()] = Teuchos::ScalarTraits<Scalar>::one();
    ;

    // 5. solve linear system for x
    Teuchos::SerialDenseVector<LO, SC> localX(Ik.size());
    Teuchos::SerialQRDenseSolver<LO, SC> qrSolver;
    qrSolver.setMatrix(Teuchos::rcp(&localA, false));
    qrSolver.setVectors(Teuchos::rcp(&localX, false), Teuchos::rcp(&ek, false));
    const int err = qrSolver.solve();
    TEUCHOS_TEST_FOR_EXCEPTION(err != 0, Exceptions::RuntimeError,
                               "MueLu::InverseApproximationFactory::GetSparseInverse: Error in serial QR solve.");

    // 6. set calculated row into Ainv
    ArrayView<const SC> Mk(localX.values(), localX.length());
    Ainv->replaceLocalValues(k, Ik, Mk);
  }
  Ainv->fillComplete();

  return Ainv;
}

}  // namespace MueLu

#endif /* MUELU_INVERSEAPPROXIMATIONFACTORY_DEF_HPP_ */
