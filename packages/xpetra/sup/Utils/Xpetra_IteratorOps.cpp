// @HEADER
//
// ***********************************************************************
//
//             Xpetra: A linear algebra interface package
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
//                    Tobias Wiesner    (tawiesn@sandia.gov)
//
// ***********************************************************************
//
// @HEADER
#include "Xpetra_IteratorOps.hpp"

namespace Xpetra {

#if defined(HAVE_XPETRA_EPETRA) && !defined(XPETRA_EPETRA_NO_32BIT_GLOBAL_INDICES)
  template<>
  void Jacobi<double,int,int,EpetraNode>(double omega,
                                         const Xpetra::Vector<double,int,int,EpetraNode> & Dinv,
                                         const Xpetra::Matrix<double,int,int,EpetraNode> & A,
                                         const Xpetra::Matrix<double,int,int,EpetraNode> & B,
                                         Xpetra::Matrix<double,int,int,EpetraNode> &C,
                                         bool call_FillComplete_on_result,
                                         bool doOptimizeStorage,
                                         const std::string & label,
                                         const Teuchos::RCP<Teuchos::ParameterList>& params) {
    typedef double        SC;
    typedef int           LO;
    typedef int           GO;
    typedef EpetraNode    NO;

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
      Epetra_CrsMatrix& epA = Xpetra::Helpers<SC,LO,GO,NO>::Op2NonConstEpetraCrs(A);
      Epetra_CrsMatrix& epB = Xpetra::Helpers<SC,LO,GO,NO>::Op2NonConstEpetraCrs(B);
      Epetra_CrsMatrix& epC = Xpetra::Helpers<SC,LO,GO,NO>::Op2NonConstEpetraCrs(C);
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
# if ((defined(EPETRA_HAVE_OMP) && (!defined(HAVE_TPETRA_INST_OPENMP) || !defined(HAVE_TPETRA_INST_INT_INT))) || \
     (!defined(EPETRA_HAVE_OMP) && (!defined(HAVE_TPETRA_INST_SERIAL) || !defined(HAVE_TPETRA_INST_INT_INT))))
      throw(Xpetra::Exceptions::RuntimeError("Xpetra must be compiled with Tpetra GO=<double,int,int> enabled."));
# else
      const Tpetra::CrsMatrix<SC,LO,GO,NO>    & tpA = Xpetra::Helpers<SC,LO,GO,NO>::Op2TpetraCrs(A);
      const Tpetra::CrsMatrix<SC,LO,GO,NO>    & tpB = Xpetra::Helpers<SC,LO,GO,NO>::Op2TpetraCrs(B);
      Tpetra::CrsMatrix<SC,LO,GO,NO>          & tpC = Xpetra::Helpers<SC,LO,GO,NO>::Op2NonConstTpetraCrs(C);
      const RCP<Tpetra::Vector<SC,LO,GO,NO> > & tpD = toTpetra(Dinv);
      Tpetra::MatrixMatrix::Jacobi(omega, *tpD, tpA, tpB, tpC, haveMultiplyDoFillComplete, label, params);
# endif
    }

    if (call_FillComplete_on_result && !haveMultiplyDoFillComplete) {
      RCP<Teuchos::ParameterList> ppp = rcp(new Teuchos::ParameterList());
      ppp->set("Optimize Storage", doOptimizeStorage);
      C.fillComplete(B.getDomainMap(), B.getRangeMap(), ppp);
    }

    // transfer striding information
    Teuchos::RCP<Xpetra::Matrix<SC,LO,GO,NO> > rcpA = Teuchos::rcp_const_cast<Xpetra::Matrix<SC,LO,GO,NO> >(Teuchos::rcpFromRef(A));
    Teuchos::RCP<Xpetra::Matrix<SC,LO,GO,NO> > rcpB = Teuchos::rcp_const_cast<Xpetra::Matrix<SC,LO,GO,NO> >(Teuchos::rcpFromRef(B));
    C.CreateView("stridedMaps", rcpA, false, rcpB, false); // TODO use references instead of RCPs
  }
#endif

#if defined(HAVE_XPETRA_EPETRA) && !defined(XPETRA_EPETRA_NO_64BIT_GLOBAL_INDICES)
  template<>
  void Jacobi<double,int,long long,EpetraNode>(double omega,
                                               const Xpetra::Vector<double,int,long long,EpetraNode> & Dinv,
                                               const Xpetra::Matrix<double,int,long long,EpetraNode> & A,
                                               const Xpetra::Matrix<double,int,long long,EpetraNode> & B,
                                               Xpetra::Matrix<double,int,long long,EpetraNode> &C,
                                               bool call_FillComplete_on_result,
                                               bool doOptimizeStorage,
                                               const std::string & label,
                                               const Teuchos::RCP<Teuchos::ParameterList>& params) {
    typedef double        SC;
    typedef int           LO;
    typedef long long     GO;
    typedef EpetraNode    NO;

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
      Epetra_CrsMatrix& epA = Xpetra::Helpers<SC,LO,GO,NO>::Op2NonConstEpetraCrs(A);
      Epetra_CrsMatrix& epB = Xpetra::Helpers<SC,LO,GO,NO>::Op2NonConstEpetraCrs(B);
      Epetra_CrsMatrix& epC = Xpetra::Helpers<SC,LO,GO,NO>::Op2NonConstEpetraCrs(C);
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
# if ((defined(EPETRA_HAVE_OMP) && (!defined(HAVE_TPETRA_INST_OPENMP) || !defined(HAVE_TPETRA_INST_INT_LONG_LONG))) || \
     (!defined(EPETRA_HAVE_OMP) && (!defined(HAVE_TPETRA_INST_SERIAL) || !defined(HAVE_TPETRA_INST_INT_LONG_LONG))))
      throw(Xpetra::Exceptions::RuntimeError("Xpetra must be compiled with Tpetra GO=<double,int,long long> enabled."));
# else
      const Tpetra::CrsMatrix<SC,LO,GO,NO>    & tpA = Xpetra::Helpers<SC,LO,GO,NO>::Op2TpetraCrs(A);
      const Tpetra::CrsMatrix<SC,LO,GO,NO>    & tpB = Xpetra::Helpers<SC,LO,GO,NO>::Op2TpetraCrs(B);
      Tpetra::CrsMatrix<SC,LO,GO,NO>          & tpC = Xpetra::Helpers<SC,LO,GO,NO>::Op2NonConstTpetraCrs(C);
      const RCP<Tpetra::Vector<SC,LO,GO,NO> > & tpD = toTpetra(Dinv);
      Tpetra::MatrixMatrix::Jacobi(omega, *tpD, tpA, tpB, tpC, haveMultiplyDoFillComplete, label, params);
# endif
    }

    if (call_FillComplete_on_result && !haveMultiplyDoFillComplete) {
      RCP<Teuchos::ParameterList> ppp = rcp(new Teuchos::ParameterList());
      ppp->set("Optimize Storage", doOptimizeStorage);
      C.fillComplete(B.getDomainMap(), B.getRangeMap(), ppp);
    }

    // transfer striding information
    Teuchos::RCP<Xpetra::Matrix<SC,LO,GO,NO> > rcpA = Teuchos::rcp_const_cast<Xpetra::Matrix<SC,LO,GO,NO> >(Teuchos::rcpFromRef(A));
    Teuchos::RCP<Xpetra::Matrix<SC,LO,GO,NO> > rcpB = Teuchos::rcp_const_cast<Xpetra::Matrix<SC,LO,GO,NO> >(Teuchos::rcpFromRef(B));
    C.CreateView("stridedMaps", rcpA, false, rcpB, false); // TODO use references instead of RCPs
  }
#endif

} // namespace Xpetra
