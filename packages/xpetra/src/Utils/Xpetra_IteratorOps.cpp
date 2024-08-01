// @HEADER
// *****************************************************************************
//             Xpetra: A linear algebra interface package
//
// Copyright 2012 NTESS and the Xpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Xpetra_IteratorOps.hpp"

namespace Xpetra {

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
                                          const Teuchos::RCP<Teuchos::ParameterList>& params) {
  typedef double SC;
  typedef int LO;
  typedef int GO;
  typedef EpetraNode NO;

  TEUCHOS_TEST_FOR_EXCEPTION(C.getRowMap()->isSameAs(*A.getRowMap()) == false, Exceptions::RuntimeError,
                             "XpetraExt::MatrixMatrix::Jacobi: row map of C is not same as row map of A")
  TEUCHOS_TEST_FOR_EXCEPTION(C.getRowMap()->isSameAs(*B.getRowMap()) == false, Exceptions::RuntimeError,
                             "XpetraExt::MatrixMatrix::Jacobi: row map of C is not same as row map of B");
  TEUCHOS_TEST_FOR_EXCEPTION(!A.isFillComplete(), Exceptions::RuntimeError, "A is not fill-completed");
  TEUCHOS_TEST_FOR_EXCEPTION(!B.isFillComplete(), Exceptions::RuntimeError, "B is not fill-completed");

  bool haveMultiplyDoFillComplete = call_FillComplete_on_result && doOptimizeStorage;

  if (C.getRowMap()->lib() == Xpetra::UseEpetra) {
#ifndef HAVE_XPETRA_EPETRAEXT
    throw(Xpetra::Exceptions::RuntimeError("Xpetra::IteratorOps::Jacobi requires EpetraExt to be compiled."));
#else
    Epetra_CrsMatrix& epA = Xpetra::Helpers<SC, LO, GO, NO>::Op2NonConstEpetraCrs(A);
    Epetra_CrsMatrix& epB = Xpetra::Helpers<SC, LO, GO, NO>::Op2NonConstEpetraCrs(B);
    Epetra_CrsMatrix& epC = Xpetra::Helpers<SC, LO, GO, NO>::Op2NonConstEpetraCrs(C);
    // FIXME
    XPETRA_DYNAMIC_CAST(const EpetraVectorT<GO XPETRA_COMMA NO>, Dinv, epD, "Xpetra::IteratorOps::Jacobi() only accepts Xpetra::EpetraVector as input argument.");

    int i = EpetraExt::MatrixMatrix::Jacobi(omega, *epD.getEpetra_Vector(), epA, epB, epC, haveMultiplyDoFillComplete);
    if (haveMultiplyDoFillComplete) {
      // Due to Epetra wrapper intricacies, we need to explicitly call
      // fillComplete on Xpetra matrix here. Specifically, EpetraCrsMatrix
      // only keeps an internal variable to check whether we are in resumed
      // state or not, but never touches the underlying Epetra object. As
      // such, we need to explicitly update the state of Xpetra matrix to
      // that of Epetra one afterwords
      C.fillComplete();
    }

    if (i != 0) {
      std::ostringstream buf;
      buf << i;
      std::string msg = "EpetraExt::MatrixMatrix::Jacobi return value of " + buf.str();
      throw(Exceptions::RuntimeError(msg));
    }
#endif
  } else if (C.getRowMap()->lib() == Xpetra::UseTpetra) {
#ifdef HAVE_XPETRA_TPETRA
#if ((defined(EPETRA_HAVE_OMP) && (!defined(HAVE_TPETRA_INST_OPENMP) || !defined(HAVE_TPETRA_INST_INT_INT))) || \
     (!defined(EPETRA_HAVE_OMP) && (!defined(HAVE_TPETRA_INST_SERIAL) || !defined(HAVE_TPETRA_INST_INT_INT))))
    throw(Xpetra::Exceptions::RuntimeError("Xpetra must be compiled with Tpetra GO=<double,int,int> enabled."));
#else
    const Tpetra::CrsMatrix<SC, LO, GO, NO>& tpA    = Xpetra::Helpers<SC, LO, GO, NO>::Op2TpetraCrs(A);
    const Tpetra::CrsMatrix<SC, LO, GO, NO>& tpB    = Xpetra::Helpers<SC, LO, GO, NO>::Op2TpetraCrs(B);
    Tpetra::CrsMatrix<SC, LO, GO, NO>& tpC          = Xpetra::Helpers<SC, LO, GO, NO>::Op2NonConstTpetraCrs(C);
    const RCP<Tpetra::Vector<SC, LO, GO, NO> >& tpD = toTpetra(Dinv);
    Tpetra::MatrixMatrix::Jacobi(omega, *tpD, tpA, tpB, tpC, haveMultiplyDoFillComplete, label, params);
#endif
#else
    throw(Xpetra::Exceptions::RuntimeError("Xpetra must be compiled with Tpetra."));
#endif
  }

  if (call_FillComplete_on_result && !haveMultiplyDoFillComplete) {
    RCP<Teuchos::ParameterList> ppp = rcp(new Teuchos::ParameterList());
    ppp->set("Optimize Storage", doOptimizeStorage);
    C.fillComplete(B.getDomainMap(), B.getRangeMap(), ppp);
  }

  // transfer striding information
  Teuchos::RCP<Xpetra::Matrix<SC, LO, GO, NO> > rcpA = Teuchos::rcp_const_cast<Xpetra::Matrix<SC, LO, GO, NO> >(Teuchos::rcpFromRef(A));
  Teuchos::RCP<Xpetra::Matrix<SC, LO, GO, NO> > rcpB = Teuchos::rcp_const_cast<Xpetra::Matrix<SC, LO, GO, NO> >(Teuchos::rcpFromRef(B));
  C.CreateView("stridedMaps", rcpA, false, rcpB, false);  // TODO use references instead of RCPs
}
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
                                                const Teuchos::RCP<Teuchos::ParameterList>& params) {
  typedef double SC;
  typedef int LO;
  typedef long long GO;
  typedef EpetraNode NO;

  TEUCHOS_TEST_FOR_EXCEPTION(C.getRowMap()->isSameAs(*A.getRowMap()) == false, Exceptions::RuntimeError,
                             "XpetraExt::MatrixMatrix::Jacobi: row map of C is not same as row map of A")
  TEUCHOS_TEST_FOR_EXCEPTION(C.getRowMap()->isSameAs(*B.getRowMap()) == false, Exceptions::RuntimeError,
                             "XpetraExt::MatrixMatrix::Jacobi: row map of C is not same as row map of B");
  TEUCHOS_TEST_FOR_EXCEPTION(!A.isFillComplete(), Exceptions::RuntimeError, "A is not fill-completed");
  TEUCHOS_TEST_FOR_EXCEPTION(!B.isFillComplete(), Exceptions::RuntimeError, "B is not fill-completed");

  bool haveMultiplyDoFillComplete = call_FillComplete_on_result && doOptimizeStorage;

  if (C.getRowMap()->lib() == Xpetra::UseEpetra) {
#ifndef HAVE_XPETRA_EPETRAEXT
    throw(Xpetra::Exceptions::RuntimeError("Xpetra::IteratorOps::Jacobi requires EpetraExt to be compiled."));
#else
    Epetra_CrsMatrix& epA = Xpetra::Helpers<SC, LO, GO, NO>::Op2NonConstEpetraCrs(A);
    Epetra_CrsMatrix& epB = Xpetra::Helpers<SC, LO, GO, NO>::Op2NonConstEpetraCrs(B);
    Epetra_CrsMatrix& epC = Xpetra::Helpers<SC, LO, GO, NO>::Op2NonConstEpetraCrs(C);
    // FIXME
    XPETRA_DYNAMIC_CAST(const EpetraVectorT<GO XPETRA_COMMA NO>, Dinv, epD, "Xpetra::IteratorOps::Jacobi() only accepts Xpetra::EpetraVector as input argument.");

    int i = EpetraExt::MatrixMatrix::Jacobi(omega, *epD.getEpetra_Vector(), epA, epB, epC, haveMultiplyDoFillComplete);
    if (haveMultiplyDoFillComplete) {
      // Due to Epetra wrapper intricacies, we need to explicitly call
      // fillComplete on Xpetra matrix here. Specifically, EpetraCrsMatrix
      // only keeps an internal variable to check whether we are in resumed
      // state or not, but never touches the underlying Epetra object. As
      // such, we need to explicitly update the state of Xpetra matrix to
      // that of Epetra one afterwords
      C.fillComplete();
    }

    if (i != 0) {
      std::ostringstream buf;
      buf << i;
      std::string msg = "EpetraExt::MatrixMatrix::Jacobi return value of " + buf.str();
      throw(Exceptions::RuntimeError(msg));
    }
#endif
  } else if (C.getRowMap()->lib() == Xpetra::UseTpetra) {
#ifdef HAVE_XPETRA_TPETRA
#if ((defined(EPETRA_HAVE_OMP) && (!defined(HAVE_TPETRA_INST_OPENMP) || !defined(HAVE_TPETRA_INST_INT_LONG_LONG))) || \
     (!defined(EPETRA_HAVE_OMP) && (!defined(HAVE_TPETRA_INST_SERIAL) || !defined(HAVE_TPETRA_INST_INT_LONG_LONG))))
    throw(Xpetra::Exceptions::RuntimeError("Xpetra must be compiled with Tpetra GO=<double,int,long long> enabled."));
#else
    const Tpetra::CrsMatrix<SC, LO, GO, NO>& tpA    = Xpetra::Helpers<SC, LO, GO, NO>::Op2TpetraCrs(A);
    const Tpetra::CrsMatrix<SC, LO, GO, NO>& tpB    = Xpetra::Helpers<SC, LO, GO, NO>::Op2TpetraCrs(B);
    Tpetra::CrsMatrix<SC, LO, GO, NO>& tpC          = Xpetra::Helpers<SC, LO, GO, NO>::Op2NonConstTpetraCrs(C);
    const RCP<Tpetra::Vector<SC, LO, GO, NO> >& tpD = toTpetra(Dinv);
    Tpetra::MatrixMatrix::Jacobi(omega, *tpD, tpA, tpB, tpC, haveMultiplyDoFillComplete, label, params);
#endif
#else
    throw(Xpetra::Exceptions::RuntimeError("Xpetra must be compiled with Tpetra."));
#endif
  }

  if (call_FillComplete_on_result && !haveMultiplyDoFillComplete) {
    RCP<Teuchos::ParameterList> ppp = rcp(new Teuchos::ParameterList());
    ppp->set("Optimize Storage", doOptimizeStorage);
    C.fillComplete(B.getDomainMap(), B.getRangeMap(), ppp);
  }

  // transfer striding information
  Teuchos::RCP<Xpetra::Matrix<SC, LO, GO, NO> > rcpA = Teuchos::rcp_const_cast<Xpetra::Matrix<SC, LO, GO, NO> >(Teuchos::rcpFromRef(A));
  Teuchos::RCP<Xpetra::Matrix<SC, LO, GO, NO> > rcpB = Teuchos::rcp_const_cast<Xpetra::Matrix<SC, LO, GO, NO> >(Teuchos::rcpFromRef(B));
  C.CreateView("stridedMaps", rcpA, false, rcpB, false);  // TODO use references instead of RCPs
}
#endif

}  // namespace Xpetra
