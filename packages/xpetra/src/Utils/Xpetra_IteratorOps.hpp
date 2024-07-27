// @HEADER
// *****************************************************************************
//             Xpetra: A linear algebra interface package
//
// Copyright 2012 NTESS and the Xpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef PACKAGES_XPETRA_SUP_UTILS_XPETRA_ITERATOROPS_HPP_
#define PACKAGES_XPETRA_SUP_UTILS_XPETRA_ITERATOROPS_HPP_

#include "Xpetra_ConfigDefs.hpp"

#include "Xpetra_Matrix.hpp"
#include "Xpetra_MatrixMatrix.hpp"
#include "Xpetra_CrsMatrixWrap.hpp"

#include "Xpetra_Map.hpp"
#include "Xpetra_StridedMap.hpp"
#include "Xpetra_StridedMapFactory.hpp"
#include "Xpetra_MapExtractor.hpp"
#include "Xpetra_MatrixFactory.hpp"
#include "Xpetra_BlockedCrsMatrix.hpp"

namespace Xpetra {

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

  TEUCHOS_TEST_FOR_EXCEPTION(C.getRowMap()->isSameAs(*A.getRowMap()) == false, Exceptions::RuntimeError,
                             "XpetraExt::MatrixMatrix::Jacobi: row map of C is not same as row map of A")
  TEUCHOS_TEST_FOR_EXCEPTION(C.getRowMap()->isSameAs(*B.getRowMap()) == false, Exceptions::RuntimeError,
                             "XpetraExt::MatrixMatrix::Jacobi: row map of C is not same as row map of B");
  TEUCHOS_TEST_FOR_EXCEPTION(!A.isFillComplete(), Exceptions::RuntimeError, "A is not fill-completed");
  TEUCHOS_TEST_FOR_EXCEPTION(!B.isFillComplete(), Exceptions::RuntimeError, "B is not fill-completed");

  bool haveMultiplyDoFillComplete = call_FillComplete_on_result && doOptimizeStorage;

  if (C.getRowMap()->lib() == Xpetra::UseEpetra) {
#ifndef HAVE_XPETRA_EPETRAEXT
    throw(Xpetra::Exceptions::RuntimeError("Xpetra::MatrixMatrix::Jacobi requires EpetraExt to be compiled."));
#else
    throw(Xpetra::Exceptions::RuntimeError("Xpetra::MatrixMatrix::Jacobi requires you to use an Epetra-compatible data type."));
#endif
  } else if (C.getRowMap()->lib() == Xpetra::UseTpetra) {
#ifdef HAVE_XPETRA_TPETRA
    const Tpetra::CrsMatrix<SC, LO, GO, NO>& tpA    = Xpetra::Helpers<SC, LO, GO, NO>::Op2TpetraCrs(A);
    const Tpetra::CrsMatrix<SC, LO, GO, NO>& tpB    = Xpetra::Helpers<SC, LO, GO, NO>::Op2TpetraCrs(B);
    Tpetra::CrsMatrix<SC, LO, GO, NO>& tpC          = Xpetra::Helpers<SC, LO, GO, NO>::Op2NonConstTpetraCrs(C);
    const RCP<Tpetra::Vector<SC, LO, GO, NO> >& tpD = toTpetra(Dinv);
    Tpetra::MatrixMatrix::Jacobi(omega, *tpD, tpA, tpB, tpC, haveMultiplyDoFillComplete, label, params);
#else
    throw(Xpetra::Exceptions::RuntimeError("Xpetra must be compiled with Tpetra."));
#endif
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

#if defined(HAVE_XPETRA_EPETRA) && !defined(XPETRA_EPETRA_NO_32BIT_GLOBAL_INDICES)
template <>
void Jacobi<double, int, int, EpetraNode>(double omega,
                                          const Xpetra::Vector<double, int, int, EpetraNode>& Dinv,
                                          const Xpetra::Matrix<double, int, int, EpetraNode>& A,
                                          const Xpetra::Matrix<double, int, int, EpetraNode>& B,
                                          Xpetra::Matrix<double, int, int, EpetraNode>& C,
                                          bool call_FillComplete_on_result,
                                          bool doOptimizeStorage,
                                          const std::string& label,
                                          const Teuchos::RCP<Teuchos::ParameterList>& params);
#endif

#if defined(HAVE_XPETRA_EPETRA) && !defined(XPETRA_EPETRA_NO_64BIT_GLOBAL_INDICES)
template <>
void Jacobi<double, int, long long, EpetraNode>(double omega,
                                                const Xpetra::Vector<double, int, long long, EpetraNode>& Dinv,
                                                const Xpetra::Matrix<double, int, long long, EpetraNode>& A,
                                                const Xpetra::Matrix<double, int, long long, EpetraNode>& B,
                                                Xpetra::Matrix<double, int, long long, EpetraNode>& C,
                                                bool call_FillComplete_on_result,
                                                bool doOptimizeStorage,
                                                const std::string& label,
                                                const Teuchos::RCP<Teuchos::ParameterList>& params);
#endif

/*!
  @class IteratorOps
  @brief Xpetra utility class containing iteration operators.

  Currently it only contains routines to generate the Jacobi iteration operator

*/
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
class IteratorOps {
#undef XPETRA_ITERATOROPS_SHORT
#include "Xpetra_UseShortNames.hpp"
 public:
  static RCP<Matrix>
  Jacobi(SC omega, const Vector& Dinv, const Matrix& A, const Matrix& B, RCP<Matrix> C_in, Teuchos::FancyOStream& fos, const std::string& label, RCP<ParameterList>& params) {
    TEUCHOS_TEST_FOR_EXCEPTION(!A.isFillComplete(), Exceptions::RuntimeError, "A is not fill-completed");
    TEUCHOS_TEST_FOR_EXCEPTION(!B.isFillComplete(), Exceptions::RuntimeError, "B is not fill-completed");

    RCP<Matrix> C = C_in;
    if (C == Teuchos::null) {
      C = MatrixFactory::Build(B.getRowMap(), Teuchos::OrdinalTraits<LO>::zero());
    } else {
      C->resumeFill();  // why this is not done inside of Tpetra Jacobi?
      fos << "Reuse C pattern" << std::endl;
    }

    Xpetra::Jacobi<Scalar, LocalOrdinal, GlobalOrdinal, Node>(omega, Dinv, A, B, *C, true, true, label, params);
    C->CreateView("stridedMaps", rcpFromRef(A), false, rcpFromRef(B), false);

    return C;
  }
};

}  // end namespace Xpetra

#define XPETRA_ITERATOROPS_SHORT

#endif /* PACKAGES_XPETRA_SUP_UTILS_XPETRA_ITERATOROPS_HPP_ */
