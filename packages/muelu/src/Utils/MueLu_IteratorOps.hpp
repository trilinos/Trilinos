// @HEADER
// *****************************************************************************
//             MueLu: A linear algebra interface package
//
// Copyright 2012 NTESS and the MueLu contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef MUELU_ITERATOROPS_HPP_
#define MUELU_ITERATOROPS_HPP_

#include "MueLu_ConfigDefs.hpp"

#include "Xpetra_Matrix.hpp"
#include "Xpetra_MatrixMatrix.hpp"
#include "Xpetra_CrsMatrixWrap.hpp"

#include "Xpetra_Map.hpp"
#include "Xpetra_StridedMap.hpp"
#include "Xpetra_StridedMapFactory.hpp"
#include "Xpetra_MapExtractor.hpp"
#include "Xpetra_MatrixFactory.hpp"
#include "Xpetra_BlockedCrsMatrix.hpp"

#include "MueLu_Exceptions.hpp"

namespace MueLu {

// General implementation
// Epetra variant throws
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void Jacobi(
    Scalar omega,
    const Xpetra::Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node>& Dinv,
    const Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>& A,
    const Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>& B,
    Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>& C,
    bool call_FillComplete_on_result                   = true,
    bool doOptimizeStorage                             = true,
    const std::string& label                           = std::string(),
    const Teuchos::RCP<Teuchos::ParameterList>& params = Teuchos::null) {
  typedef Scalar SC;
  typedef LocalOrdinal LO;
  typedef GlobalOrdinal GO;
  typedef Node NO;

  TEUCHOS_TEST_FOR_EXCEPTION(C.getRowMap()->isSameAs(*A.getRowMap()) == false, MueLu::Exceptions::RuntimeError,
                             "MueLu::Jacobi: row map of C is not same as row map of A")
  TEUCHOS_TEST_FOR_EXCEPTION(C.getRowMap()->isSameAs(*B.getRowMap()) == false, MueLu::Exceptions::RuntimeError,
                             "MueLu::Jacobi: row map of C is not same as row map of B");
  TEUCHOS_TEST_FOR_EXCEPTION(!A.isFillComplete(), MueLu::Exceptions::RuntimeError, "A is not fill-completed");
  TEUCHOS_TEST_FOR_EXCEPTION(!B.isFillComplete(), MueLu::Exceptions::RuntimeError, "B is not fill-completed");

  bool haveMultiplyDoFillComplete = call_FillComplete_on_result && doOptimizeStorage;

  if (C.getRowMap()->lib() == Xpetra::UseTpetra) {
    const Tpetra::CrsMatrix<SC, LO, GO, NO>& tpA    = Xpetra::Helpers<SC, LO, GO, NO>::Op2TpetraCrs(A);
    const Tpetra::CrsMatrix<SC, LO, GO, NO>& tpB    = Xpetra::Helpers<SC, LO, GO, NO>::Op2TpetraCrs(B);
    Tpetra::CrsMatrix<SC, LO, GO, NO>& tpC          = Xpetra::Helpers<SC, LO, GO, NO>::Op2NonConstTpetraCrs(C);
    const RCP<Tpetra::Vector<SC, LO, GO, NO> >& tpD = toTpetra(Dinv);
    Tpetra::MatrixMatrix::Jacobi(omega, *tpD, tpA, tpB, tpC, haveMultiplyDoFillComplete, label, params);
  }

  if (call_FillComplete_on_result && !haveMultiplyDoFillComplete) {
    RCP<Teuchos::ParameterList> fillParams = rcp(new Teuchos::ParameterList());
    fillParams->set("Optimize Storage", doOptimizeStorage);
    C.fillComplete(B.getDomainMap(), B.getRangeMap(), fillParams);
  }

  // transfer striding information
  RCP<Xpetra::Matrix<SC, LO, GO, NO> > rcpA = Teuchos::rcp_const_cast<Xpetra::Matrix<SC, LO, GO, NO> >(Teuchos::rcpFromRef(A));
  RCP<Xpetra::Matrix<SC, LO, GO, NO> > rcpB = Teuchos::rcp_const_cast<Xpetra::Matrix<SC, LO, GO, NO> >(Teuchos::rcpFromRef(B));
  C.CreateView("stridedMaps", rcpA, false, rcpB, false);  // TODO use references instead of RCPs
}  // end Jacobi

/*!
  @class IteratorOps
  @brief MueLu utility class containing iteration operators.

  Currently it only contains routines to generate the Jacobi iteration operator

*/
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
class IteratorOps {
 public:
  static RCP<Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> >
  Jacobi(Scalar omega, const Xpetra::Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node>& Dinv, const Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>& A, const Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>& B, RCP<Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> > C_in, Teuchos::FancyOStream& fos, const std::string& label, RCP<ParameterList>& params) {
    TEUCHOS_TEST_FOR_EXCEPTION(!A.isFillComplete(), MueLu::Exceptions::RuntimeError, "A is not fill-completed");
    TEUCHOS_TEST_FOR_EXCEPTION(!B.isFillComplete(), MueLu::Exceptions::RuntimeError, "B is not fill-completed");

    RCP<Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> > C = C_in;
    if (C == Teuchos::null) {
      C = Xpetra::MatrixFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Build(B.getRowMap(), Teuchos::OrdinalTraits<LocalOrdinal>::zero());
    } else {
      C->resumeFill();  // why this is not done inside of Tpetra Jacobi?
      fos << "Reuse C pattern" << std::endl;
    }

    MueLu::Jacobi<Scalar, LocalOrdinal, GlobalOrdinal, Node>(omega, Dinv, A, B, *C, true, true, label, params);
    C->CreateView("stridedMaps", rcpFromRef(A), false, rcpFromRef(B), false);

    return C;
  }
};

}  // end namespace MueLu

#define MUELU_ITERATOROPS_SHORT

#endif /* PACKAGES_MUELU_SUP_UTILS_XPETRA_ITERATOROPS_HPP_ */
