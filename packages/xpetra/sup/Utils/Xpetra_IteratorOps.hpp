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


template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void Jacobi(
    Scalar omega,
    const Xpetra::Vector<Scalar,LocalOrdinal,GlobalOrdinal,Node> & Dinv,
    const Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>& A,
    const Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>& B,
    Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>& C,
    bool call_FillComplete_on_result = true,
    bool doOptimizeStorage = true,
    const std::string & label = std::string()) {

  if(C.getRowMap()->isSameAs(*A.getRowMap()) == false) {
    std::string msg = "XpetraExt::MatrixMatrix::Jacobi: row map of C is not same as row map of A";
    throw(Xpetra::Exceptions::RuntimeError(msg));
  }
  else if(C.getRowMap()->isSameAs(*B.getRowMap()) == false) {
    std::string msg = "XpetraExt::MatrixMatrix::Jacobi: row map of C is not same as row map of B";
    throw(Xpetra::Exceptions::RuntimeError(msg));
  }

  if (!A.isFillComplete())
    throw(Xpetra::Exceptions::RuntimeError("A is not fill-completed"));
  if (!B.isFillComplete())
    throw(Xpetra::Exceptions::RuntimeError("B is not fill-completed"));

  bool haveMultiplyDoFillComplete = call_FillComplete_on_result && doOptimizeStorage;

  if (C.getRowMap()->lib() == Xpetra::UseEpetra) {
#ifndef HAVE_XPETRA_EPETRAEXT
    throw(Xpetra::Exceptions::RuntimeError("Xpetra::MatrixMatrix::Jacobi requires EpetraExt to be compiled."));
#else
    throw(Xpetra::Exceptions::RuntimeError("Xpetra::MatrixMatrix::Jacobi requires you to use an Epetra-compatible data type."));
#endif
  } else if (C.getRowMap()->lib() == Xpetra::UseTpetra) {
#ifdef HAVE_XPETRA_TPETRA
    const Tpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> & tpA = Xpetra::Helpers<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Op2TpetraCrs(A);
    const Tpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> & tpB = Xpetra::Helpers<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Op2TpetraCrs(B);
    Tpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>       & tpC = Xpetra::Helpers<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Op2NonConstTpetraCrs(C);
    const RCP<Tpetra::Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node>  >          & tpD = toTpetra(Dinv);
    Tpetra::MatrixMatrix::Jacobi(omega,*tpD,tpA,tpB,tpC,haveMultiplyDoFillComplete,label);
#else
    throw(Xpetra::Exceptions::RuntimeError("Xpetra must be compiled with Tpetra."));
#endif
  }

  if(call_FillComplete_on_result && !haveMultiplyDoFillComplete) {
    RCP<Teuchos::ParameterList> params = rcp(new Teuchos::ParameterList());
    params->set("Optimize Storage",doOptimizeStorage);
    C.fillComplete(B.getDomainMap(),B.getRangeMap(),params);
  }

  // transfer striding information
  RCP<Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> > rcpA = Teuchos::rcp_const_cast<Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> >(Teuchos::rcpFromRef(A));
  RCP<Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> > rcpB = Teuchos::rcp_const_cast<Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> >(Teuchos::rcpFromRef(B));
  C.CreateView("stridedMaps", rcpA, false, rcpB, false); // TODO use references instead of RCPs
} // end Jacobi


template <class GlobalOrdinal, class Node>
inline void JacobiT(
    double omega,
    const Xpetra::Vector<double,int,GlobalOrdinal,Node> & Dinv,
    const Xpetra::Matrix<double,int,GlobalOrdinal,Node> & A,
    const Xpetra::Matrix<double,int,GlobalOrdinal,Node> & B,
    Xpetra::Matrix<double,int,GlobalOrdinal,Node> &C,
    bool call_FillComplete_on_result,
    bool doOptimizeStorage,
    const std::string & label) {

  typedef double        SC;
  typedef int           LO;
  typedef GlobalOrdinal GO;
  typedef Node          NO;

  if(C.getRowMap()->isSameAs(*A.getRowMap()) == false) {
    std::string msg = "XpetraExt::MatrixMatrix::Jacobi: row map of C is not same as row map of A";
    throw(Xpetra::Exceptions::RuntimeError(msg));
  }
  else if(C.getRowMap()->isSameAs(*B.getRowMap()) == false) {
    std::string msg = "XpetraExt::MatrixMatrix::Jacobi: row map of C is not same as row map of B";
    throw(Xpetra::Exceptions::RuntimeError(msg));
  }

  if (!A.isFillComplete())
    throw(Xpetra::Exceptions::RuntimeError("A is not fill-completed"));
  if (!B.isFillComplete())
    throw(Xpetra::Exceptions::RuntimeError("B is not fill-completed"));

  bool haveMultiplyDoFillComplete = call_FillComplete_on_result && doOptimizeStorage;

  if (C.getRowMap()->lib() == Xpetra::UseEpetra) {
#       ifndef HAVE_XPETRA_EPETRAEXT
    throw(Xpetra::Exceptions::RuntimeError("Xpetra::IteratorOps::Jacobi requires EpetraExt to be compiled."));
#else
    Epetra_CrsMatrix & epA = Xpetra::Helpers<SC,LO,GO,NO>::Op2NonConstEpetraCrs(A);
    Epetra_CrsMatrix & epB = Xpetra::Helpers<SC,LO,GO,NO>::Op2NonConstEpetraCrs(B);
    Epetra_CrsMatrix & epC = Xpetra::Helpers<SC,LO,GO,NO>::Op2NonConstEpetraCrs(C);
    //    const Epetra_Vector & epD = toEpetra(Dinv);
    XPETRA_DYNAMIC_CAST(const EpetraVectorT<GO COMMA NO>, Dinv, epD, "Xpetra::IteratorOps::Jacobi() only accepts Xpetra::EpetraVector as input argument.");

    int i = EpetraExt::MatrixMatrix::Jacobi(omega,*epD.getEpetra_Vector(),epA,epB,epC,haveMultiplyDoFillComplete);
    if (i != 0) {
      std::ostringstream buf;
      buf << i;
      std::string msg = "EpetraExt::MatrixMatrix::Jacobi return value of " + buf.str();
      throw(Exceptions::RuntimeError(msg));
    }
#endif
  } else if (C.getRowMap()->lib() == Xpetra::UseTpetra) {
#ifdef HAVE_XPETRA_TPETRA
    const Tpetra::CrsMatrix<SC, LO, GO, NO> & tpA = Xpetra::Helpers<SC,LO,GO,NO>::Op2TpetraCrs(A);
    const Tpetra::CrsMatrix<SC, LO, GO, NO> & tpB = Xpetra::Helpers<SC,LO,GO,NO>::Op2TpetraCrs(B);
    Tpetra::CrsMatrix<SC, LO, GO, NO>       & tpC = Xpetra::Helpers<SC,LO,GO,NO>::Op2NonConstTpetraCrs(C);
    const RCP<Tpetra::Vector<SC, LO, GO, NO>  >          & tpD = toTpetra(Dinv);
    Tpetra::MatrixMatrix::Jacobi(omega,*tpD,tpA,tpB,tpC,haveMultiplyDoFillComplete,label);
#else
    throw(Xpetra::Exceptions::RuntimeError("Xpetra must be compiled with Tpetra."));
#endif
  }

  if(call_FillComplete_on_result && !haveMultiplyDoFillComplete) {
    RCP<Teuchos::ParameterList> params = rcp(new Teuchos::ParameterList());
    params->set("Optimize Storage",doOptimizeStorage);
    C.fillComplete(B.getDomainMap(),B.getRangeMap(),params);
  }

  // transfer striding information
  RCP<Xpetra::Matrix<SC, LO, GO, NO> > rcpA = Teuchos::rcp_const_cast<Xpetra::Matrix<SC, LO, GO, NO> >(Teuchos::rcpFromRef(A));
  RCP<Xpetra::Matrix<SC, LO, GO, NO> > rcpB = Teuchos::rcp_const_cast<Xpetra::Matrix<SC, LO, GO, NO> >(Teuchos::rcpFromRef(B));
  C.CreateView("stridedMaps", rcpA, false, rcpB, false); // TODO use references instead of RCPs
} // end Jacobi

inline void JacobiInt(
    double omega,
    const Xpetra::Vector<double,int,int,Kokkos::Compat::KokkosSerialWrapperNode> & Dinv,
    const Xpetra::Matrix<double,int,int,Kokkos::Compat::KokkosSerialWrapperNode> & A,
    const Xpetra::Matrix<double,int,int,Kokkos::Compat::KokkosSerialWrapperNode> & B,
    Xpetra::Matrix<double,int,int,Kokkos::Compat::KokkosSerialWrapperNode> &C,
    bool call_FillComplete_on_result,
    bool doOptimizeStorage,
    const std::string & label) {

  typedef Xpetra::Vector<double, int, int, Kokkos::Compat::KokkosSerialWrapperNode>::node_type NO;

  if(C.getRowMap()->isSameAs(*A.getRowMap()) == false) {
    std::string msg = "XpetraExt::MatrixMatrix::Jacobi: row map of C is not same as row map of A";
    throw(Xpetra::Exceptions::RuntimeError(msg));
  }
  else if(C.getRowMap()->isSameAs(*B.getRowMap()) == false) {
    std::string msg = "XpetraExt::MatrixMatrix::Jacobi: row map of C is not same as row map of B";
    throw(Xpetra::Exceptions::RuntimeError(msg));
  }

  if (!A.isFillComplete())
    throw(Xpetra::Exceptions::RuntimeError("A is not fill-completed"));
  if (!B.isFillComplete())
    throw(Xpetra::Exceptions::RuntimeError("B is not fill-completed"));

  bool haveMultiplyDoFillComplete = call_FillComplete_on_result && doOptimizeStorage;

  if (C.getRowMap()->lib() == Xpetra::UseEpetra) {
#       ifndef HAVE_XPETRA_EPETRAEXT
    throw(Xpetra::Exceptions::RuntimeError("Xpetra::IteratorOps::Jacobi requires EpetraExt to be compiled."));
#else
    Epetra_CrsMatrix & epA = Xpetra::Helpers<double,int,int,Kokkos::Compat::KokkosSerialWrapperNode>::Op2NonConstEpetraCrs(A);
    Epetra_CrsMatrix & epB = Xpetra::Helpers<double,int,int,Kokkos::Compat::KokkosSerialWrapperNode>::Op2NonConstEpetraCrs(B);
    Epetra_CrsMatrix & epC = Xpetra::Helpers<double,int,int,Kokkos::Compat::KokkosSerialWrapperNode>::Op2NonConstEpetraCrs(C);
    //    const Epetra_Vector & epD = toEpetra(Dinv);
    XPETRA_DYNAMIC_CAST(const EpetraVectorT<int COMMA typename Xpetra::Map<int COMMA Kokkos::Compat::KokkosSerialWrapperNode>::node_type>, Dinv, epD, "Xpetra::IteratorOps::Jacobi() only accepts Xpetra::EpetraVector as input argument.");

    int i = EpetraExt::MatrixMatrix::Jacobi(omega,*epD.getEpetra_Vector(),epA,epB,epC,haveMultiplyDoFillComplete);
    if (i != 0) {
      std::ostringstream buf;
      buf << i;
      std::string msg = "EpetraExt::MatrixMatrix::Jacobi return value of " + buf.str();
      throw(Exceptions::RuntimeError(msg));
    }
#endif
  } else if (C.getRowMap()->lib() == Xpetra::UseTpetra) {
#ifdef HAVE_XPETRA_TPETRA
#ifdef HAVE_XPETRA_TPETRA_INST_INT_INT
    const Tpetra::CrsMatrix<double, int, int, Kokkos::Compat::KokkosSerialWrapperNode> & tpA = Xpetra::Helpers<double,int,int,Kokkos::Compat::KokkosSerialWrapperNode>::Op2TpetraCrs(A);
    const Tpetra::CrsMatrix<double, int, int, Kokkos::Compat::KokkosSerialWrapperNode> & tpB = Xpetra::Helpers<double,int,int,Kokkos::Compat::KokkosSerialWrapperNode>::Op2TpetraCrs(B);
    Tpetra::CrsMatrix<double, int, int, Kokkos::Compat::KokkosSerialWrapperNode>       & tpC = Xpetra::Helpers<double,int,int,Kokkos::Compat::KokkosSerialWrapperNode>::Op2NonConstTpetraCrs(C);
    const RCP<Tpetra::Vector<double, int, int, Kokkos::Compat::KokkosSerialWrapperNode>  >          & tpD = toTpetra(Dinv);
    Tpetra::MatrixMatrix::Jacobi(omega,*tpD,tpA,tpB,tpC,haveMultiplyDoFillComplete,label);
#else
    throw(Xpetra::Exceptions::RuntimeError("Xpetra must be compiled with Teptra GO=int enabled."));
#endif
#else
    throw(Xpetra::Exceptions::RuntimeError("Xpetra must be compiled with Tpetra."));
#endif
  }

  if(call_FillComplete_on_result && !haveMultiplyDoFillComplete) {
    RCP<Teuchos::ParameterList> params = rcp(new Teuchos::ParameterList());
    params->set("Optimize Storage",doOptimizeStorage);
    C.fillComplete(B.getDomainMap(),B.getRangeMap(),params);
  }

  // transfer striding information
  Teuchos::RCP<Xpetra::Matrix<double, int, int, Kokkos::Compat::KokkosSerialWrapperNode> > rcpA = Teuchos::rcp_const_cast<Xpetra::Matrix<double, int, int, Kokkos::Compat::KokkosSerialWrapperNode> >(Teuchos::rcpFromRef(A));
  Teuchos::RCP<Xpetra::Matrix<double, int, int, Kokkos::Compat::KokkosSerialWrapperNode> > rcpB = Teuchos::rcp_const_cast<Xpetra::Matrix<double, int, int, Kokkos::Compat::KokkosSerialWrapperNode> >(Teuchos::rcpFromRef(B));
  C.CreateView("stridedMaps", rcpA, false, rcpB, false); // TODO use references instead of RCPs
} // end Jacobi


#ifndef XPETRA_EPETRA_NO_32BIT_GLOBAL_INDICES
inline void Jacobi(
    double omega,
    const Xpetra::Vector<double,int,int,Kokkos::Compat::KokkosSerialWrapperNode> & Dinv,
    const Xpetra::Matrix<double,int,int,Kokkos::Compat::KokkosSerialWrapperNode> & A,
    const Xpetra::Matrix<double,int,int,Kokkos::Compat::KokkosSerialWrapperNode> & B,
    Xpetra::Matrix<double,int,int,Kokkos::Compat::KokkosSerialWrapperNode> &C,
    bool call_FillComplete_on_result,
    bool doOptimizeStorage,
    const std::string & label) {
  JacobiInt(omega, Dinv, A, B, C, call_FillComplete_on_result, doOptimizeStorage, label);
}
#endif

#ifdef HAVE_XPETRA_INT_LONG_LONG
inline void Jacobi(
    double omega,
    const Xpetra::Vector<double,int,long long,Kokkos::Compat::KokkosSerialWrapperNode> & Dinv,
    const Xpetra::Matrix<double,int,long long,Kokkos::Compat::KokkosSerialWrapperNode> & A,
    const Xpetra::Matrix<double,int,long long,Kokkos::Compat::KokkosSerialWrapperNode> & B,
    Xpetra::Matrix<double,int,long long,Kokkos::Compat::KokkosSerialWrapperNode> &C,
    bool call_FillComplete_on_result,
    bool doOptimizeStorage,
    const std::string & label) {
  JacobiT(omega, Dinv, A, B, C, call_FillComplete_on_result, doOptimizeStorage, label);
}
#endif // HAVE_XPETRA_INT_LONG_LONG

/*!
    @class IteratorOps
    @brief Xpetra utility class containing iteration operators.

    Currently it only contains routines to generate the Jacobi iteration operator

 */
template <class Scalar,
class LocalOrdinal  /*= int*/,
class GlobalOrdinal /*= LocalOrdinal*/,
class Node          /*= KokkosClassic::DefaultNode::DefaultNodeType*/>
class IteratorOps {

private:
#undef XPETRA_ITERATOROPS_SHORT
#include "Xpetra_UseShortNames.hpp"

public:

  static RCP<Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> >
  Jacobi(Scalar omega,
      const Vector& Dinv,
      const Matrix& A,
      const Matrix& B,
      RCP<Matrix> C_in,
      Teuchos::FancyOStream &fos,
      const std::string & label) {
    // Sanity checks
    if (!A.isFillComplete())
      throw Exceptions::RuntimeError("A is not fill-completed");
    if (!B.isFillComplete())
      throw Exceptions::RuntimeError("B is not fill-completed");

    // Default case: Xpetra Jacobi
    RCP<Matrix> C = C_in;
    if (C == Teuchos::null)
      C = MatrixFactory::Build(B.getRowMap(),Teuchos::OrdinalTraits<LO>::zero());

    Xpetra::Jacobi<Scalar, LocalOrdinal, GlobalOrdinal, Node>(omega, Dinv, A, B, *C, true,true,label);
    C->CreateView("stridedMaps", rcpFromRef(A),false, rcpFromRef(B), false);
    return C;
  } //Jacobi


};

/*!
    @class IteratorOps
    @brief Xpetra utility class containing iteration operators.

    Currently it only contains routines to generate the Jacobi iteration operator

    Specialization for SC=double and LO=int
 */
template <class GlobalOrdinal, class Node>
class IteratorOps<double,int,GlobalOrdinal,Node> {
public:
  typedef double Scalar;
  typedef int LocalOrdinal;

  static RCP<Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> >
  Jacobi(Scalar omega,
      const Xpetra::Vector<Scalar,LocalOrdinal,GlobalOrdinal,Node>& Dinv,
      const Xpetra::Matrix<Scalar,LocalOrdinal,GlobalOrdinal,Node>& A,
      const Xpetra::Matrix<Scalar,LocalOrdinal,GlobalOrdinal,Node>& B,
      RCP<Xpetra::Matrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> > C_in,
      Teuchos::FancyOStream &fos,
      const std::string & label) {
    // Sanity checks
    if (!A.isFillComplete())
      throw Exceptions::RuntimeError("A is not fill-completed");
    if (!B.isFillComplete())
      throw Exceptions::RuntimeError("B is not fill-completed");

    // Default case: Xpetra Jacobi
    RCP<Xpetra::Matrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> > C = C_in;
    if (C == Teuchos::null)
      C = Xpetra::MatrixFactory<Scalar,LocalOrdinal,GlobalOrdinal,Node>::Build(B.getRowMap(),Teuchos::OrdinalTraits<LocalOrdinal>::zero());

    JacobiT<GlobalOrdinal>(omega, Dinv, A, B, *C, true,true,label);
    C->CreateView("stridedMaps", rcpFromRef(A),false, rcpFromRef(B), false);
    return C;
  } //Jacobi
};

/*!
    @class IteratorOps
    @brief Xpetra utility class containing iteration operators.

    Currently it only contains routines to generate the Jacobi iteration operator

    Specialization for LO=GO=int and SC=double
 */
#ifndef HAVE_XPETRA_TPETRA_INST_INT_INT
template <class Node>
class IteratorOps<double,int,int,Node> {
public:
  typedef int LocalOrdinal;
  typedef int GlobalOrdinal;
  typedef double Scalar;

  static RCP<Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> >
  Jacobi(Scalar omega,
      const Xpetra::Vector<Scalar,LocalOrdinal,GlobalOrdinal,Node>& Dinv,
      const Xpetra::Matrix<Scalar,LocalOrdinal,GlobalOrdinal,Node>& A,
      const Xpetra::Matrix<Scalar,LocalOrdinal,GlobalOrdinal,Node>& B,
      RCP<Xpetra::Matrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> > C_in,
      Teuchos::FancyOStream &fos,
      const std::string & label) {
    // Sanity checks
    if (!A.isFillComplete())
      throw Exceptions::RuntimeError("A is not fill-completed");
    if (!B.isFillComplete())
      throw Exceptions::RuntimeError("B is not fill-completed");

    // Default case: Xpetra Jacobi
    RCP<Xpetra::Matrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> > C = C_in;
    if (C == Teuchos::null)
      C = Xpetra::MatrixFactory<Scalar,LocalOrdinal,GlobalOrdinal,Node>::Build(B.getRowMap(),Teuchos::OrdinalTraits<LocalOrdinal>::zero());

    Xpetra::JacobiInt(omega, Dinv, A, B, *C, true,true,label);
    C->CreateView("stridedMaps", rcpFromRef(A),false, rcpFromRef(B), false);
    return C;
  } //Jacobi
};
#endif




} // end namespace Xpetra

#define XPETRA_ITERATOROPS_SHORT

#endif /* PACKAGES_XPETRA_SUP_UTILS_XPETRA_ITERATOROPS_HPP_ */
