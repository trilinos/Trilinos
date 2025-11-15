// @HEADER
// *****************************************************************************
//             MueLu: A linear algebra interface package
//
// Copyright 2012 NTESS and the MueLu contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "MueLu_IteratorOps.hpp"

namespace MueLu {

#if defined(HAVE_MUELU_EPETRA) && !defined(XPETRA_EPETRA_NO_32BIT_GLOBAL_INDICES)
template <>
void Jacobi<double, int, int, Xpetra::EpetraNode>(double omega,
                                                  const Xpetra::Vector<double, int, int, Xpetra::EpetraNode>& Dinv,
                                                  const Xpetra::Matrix<double, int, int, Xpetra::EpetraNode>& A,
                                                  const Xpetra::Matrix<double, int, int, Xpetra::EpetraNode>& B,
                                                  Xpetra::Matrix<double, int, int, Xpetra::EpetraNode>& C,
                                                  bool call_FillComplete_on_result,
                                                  bool doOptimizeStorage,
                                                  const std::string& label,
                                                  const Teuchos::RCP<Teuchos::ParameterList>& params) {
  typedef double SC;
  typedef int LO;
  typedef int GO;
  typedef Xpetra::EpetraNode NO;

  TEUCHOS_TEST_FOR_EXCEPTION(C.getRowMap()->isSameAs(*A.getRowMap()) == false, MueLu::Exceptions::RuntimeError,
                             "MueLu::Jacobi: row map of C is not same as row map of A")
  TEUCHOS_TEST_FOR_EXCEPTION(C.getRowMap()->isSameAs(*B.getRowMap()) == false, MueLu::Exceptions::RuntimeError,
                             "MueLu::Jacobi: row map of C is not same as row map of B");
  TEUCHOS_TEST_FOR_EXCEPTION(!A.isFillComplete(), MueLu::Exceptions::RuntimeError, "A is not fill-completed");
  TEUCHOS_TEST_FOR_EXCEPTION(!B.isFillComplete(), MueLu::Exceptions::RuntimeError, "B is not fill-completed");

  bool haveMultiplyDoFillComplete = call_FillComplete_on_result && doOptimizeStorage;

  if (C.getRowMap()->lib() == Xpetra::UseEpetra) {
#ifndef HAVE_MUELU_EPETRAEXT
    throw(MueLu::Exceptions::RuntimeError("MueLu::Jacobi requires EpetraExt to be compiled."));
#else
    Epetra_CrsMatrix& epA = Xpetra::Helpers<SC, LO, GO, NO>::Op2NonConstEpetraCrs(A);
    Epetra_CrsMatrix& epB = Xpetra::Helpers<SC, LO, GO, NO>::Op2NonConstEpetraCrs(B);
    Epetra_CrsMatrix& epC = Xpetra::Helpers<SC, LO, GO, NO>::Op2NonConstEpetraCrs(C);
    // FIXME
    XPETRA_DYNAMIC_CAST(const Xpetra::EpetraVectorT<GO XPETRA_COMMA NO>, Dinv, epD, "MueLu::Jacobi() only accepts Xpetra::EpetraVector as input argument.");

    int i = EpetraExt::MatrixMatrix::Jacobi(omega, *epD.getEpetra_Vector(), epA, epB, epC, haveMultiplyDoFillComplete);
    if (haveMultiplyDoFillComplete) {
      // Due to Epetra wrapper intricacies, we need to explicitly call
      // fillComplete on MueLu matrix here. Specifically, EpetraCrsMatrix
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
#if ((defined(EPETRA_HAVE_OMP) && (!defined(HAVE_TPETRA_INST_OPENMP) || !defined(HAVE_TPETRA_INST_INT_INT))) || \
     (!defined(EPETRA_HAVE_OMP) && (!defined(HAVE_TPETRA_INST_SERIAL) || !defined(HAVE_TPETRA_INST_INT_INT))))
    throw(MueLu::Exceptions::RuntimeError("Xpetra must be compiled with Tpetra GO=<double,int,int> enabled."));
#else
    const Tpetra::CrsMatrix<SC, LO, GO, NO>& tpA    = Xpetra::Helpers<SC, LO, GO, NO>::Op2TpetraCrs(A);
    const Tpetra::CrsMatrix<SC, LO, GO, NO>& tpB    = Xpetra::Helpers<SC, LO, GO, NO>::Op2TpetraCrs(B);
    Tpetra::CrsMatrix<SC, LO, GO, NO>& tpC          = Xpetra::Helpers<SC, LO, GO, NO>::Op2NonConstTpetraCrs(C);
    const RCP<Tpetra::Vector<SC, LO, GO, NO> >& tpD = toTpetra(Dinv);
    Tpetra::MatrixMatrix::Jacobi(omega, *tpD, tpA, tpB, tpC, haveMultiplyDoFillComplete, label, params);
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

#if defined(HAVE_MUELU_EPETRA) && !defined(XPETRA_EPETRA_NO_64BIT_GLOBAL_INDICES)
template <>
void Jacobi<double, int, long long, Xpetra::EpetraNode>(double omega,
                                                        const Xpetra::Vector<double, int, long long, Xpetra::EpetraNode>& Dinv,
                                                        const Xpetra::Matrix<double, int, long long, Xpetra::EpetraNode>& A,
                                                        const Xpetra::Matrix<double, int, long long, Xpetra::EpetraNode>& B,
                                                        Xpetra::Matrix<double, int, long long, Xpetra::EpetraNode>& C,
                                                        bool call_FillComplete_on_result,
                                                        bool doOptimizeStorage,
                                                        const std::string& label,
                                                        const Teuchos::RCP<Teuchos::ParameterList>& params) {
  typedef double SC;
  typedef int LO;
  typedef long long GO;
  typedef Xpetra::EpetraNode NO;

  TEUCHOS_TEST_FOR_EXCEPTION(C.getRowMap()->isSameAs(*A.getRowMap()) == false, MueLu::Exceptions::RuntimeError,
                             "MueLu::Jacobi: row map of C is not same as row map of A")
  TEUCHOS_TEST_FOR_EXCEPTION(C.getRowMap()->isSameAs(*B.getRowMap()) == false, MueLu::Exceptions::RuntimeError,
                             "MueLu::Jacobi: row map of C is not same as row map of B");
  TEUCHOS_TEST_FOR_EXCEPTION(!A.isFillComplete(), MueLu::Exceptions::RuntimeError, "A is not fill-completed");
  TEUCHOS_TEST_FOR_EXCEPTION(!B.isFillComplete(), MueLu::Exceptions::RuntimeError, "B is not fill-completed");

  bool haveMultiplyDoFillComplete = call_FillComplete_on_result && doOptimizeStorage;

  if (C.getRowMap()->lib() == Xpetra::UseEpetra) {
#ifndef HAVE_MUELU_EPETRAEXT
    throw(MueLu::Exceptions::RuntimeError("MueLu::Jacobi requires EpetraExt to be compiled."));
#else
    Epetra_CrsMatrix& epA = Xpetra::Helpers<SC, LO, GO, NO>::Op2NonConstEpetraCrs(A);
    Epetra_CrsMatrix& epB = Xpetra::Helpers<SC, LO, GO, NO>::Op2NonConstEpetraCrs(B);
    Epetra_CrsMatrix& epC = Xpetra::Helpers<SC, LO, GO, NO>::Op2NonConstEpetraCrs(C);
    // FIXME
    XPETRA_DYNAMIC_CAST(const Xpetra::EpetraVectorT<GO XPETRA_COMMA NO>, Dinv, epD, "MueLu::Jacobi() only accepts Xpetra::EpetraVector as input argument.");

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
#if ((defined(EPETRA_HAVE_OMP) && (!defined(HAVE_TPETRA_INST_OPENMP) || !defined(HAVE_TPETRA_INST_INT_LONG_LONG))) || \
     (!defined(EPETRA_HAVE_OMP) && (!defined(HAVE_TPETRA_INST_SERIAL) || !defined(HAVE_TPETRA_INST_INT_LONG_LONG))))
    throw(MueLu::Exceptions::RuntimeError("Xpetra must be compiled with Tpetra GO=<double,int,long long> enabled."));
#else
    const Tpetra::CrsMatrix<SC, LO, GO, NO>& tpA    = Xpetra::Helpers<SC, LO, GO, NO>::Op2TpetraCrs(A);
    const Tpetra::CrsMatrix<SC, LO, GO, NO>& tpB    = Xpetra::Helpers<SC, LO, GO, NO>::Op2TpetraCrs(B);
    Tpetra::CrsMatrix<SC, LO, GO, NO>& tpC          = Xpetra::Helpers<SC, LO, GO, NO>::Op2NonConstTpetraCrs(C);
    const RCP<Tpetra::Vector<SC, LO, GO, NO> >& tpD = toTpetra(Dinv);
    Tpetra::MatrixMatrix::Jacobi(omega, *tpD, tpA, tpB, tpC, haveMultiplyDoFillComplete, label, params);
#endif
  }

  if (call_FillComplete_on_result && !haveMultiplyDoFillComplete) {
    RCP<Teuchos::ParameterList> ppp = rcp(new Teuchos::ParameterList());
    ppp->set("Optimize Storage", doOptimizeStorage);
    C.fillComplete(B.getDomainMap(), B.getRangeMap(), ppp);
  }

  // transfer striding information
  Teuchos::RCP<Xpetra::Matrix<SC, LO, GO, NO> > rcpA  = Teuchos::rcp_const_cast<Xpetra::Matrix<SC, LO, GO, NO> >(Teuchos::rcpFromRef(A));
  Teuchos::RCP<Xpetra : Matrix<SC, LO, GO, NO> > rcpB = Teuchos::rcp_const_cast<Xpetra::Matrix<SC, LO, GO, NO> >(Teuchos::rcpFromRef(B));
  C.CreateView("stridedMaps", rcpA, false, rcpB, false);  // TODO use references instead of RCPs
}
#endif

}  // namespace MueLu
