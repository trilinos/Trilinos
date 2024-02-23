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
#ifndef THYRA_MUELU_PRECONDITIONER_FACTORY_DEF_HPP
#define THYRA_MUELU_PRECONDITIONER_FACTORY_DEF_HPP

#include "Thyra_MueLuPreconditionerFactory_decl.hpp"
#include "MueLu_TpetraOperatorAsRowMatrix.hpp"

#if defined(HAVE_MUELU_STRATIMIKOS) && defined(HAVE_MUELU_THYRA)

// This is not as general as possible, but should be good enough for most builds.
#if ((defined(HAVE_TPETRA_INST_DOUBLE) && defined(HAVE_TPETRA_INST_FLOAT) && !defined(HAVE_TPETRA_INST_COMPLEX_DOUBLE) && !defined(HAVE_TPETRA_INST_COMPLEX_FLOAT)) || \
     (!defined(HAVE_TPETRA_INST_DOUBLE) && !defined(HAVE_TPETRA_INST_FLOAT) && defined(HAVE_TPETRA_INST_COMPLEX_DOUBLE) && defined(HAVE_TPETRA_INST_COMPLEX_FLOAT)) || \
     (defined(HAVE_TPETRA_INST_DOUBLE) && defined(HAVE_TPETRA_INST_FLOAT) && defined(HAVE_TPETRA_INST_COMPLEX_DOUBLE) && defined(HAVE_TPETRA_INST_COMPLEX_FLOAT)))
#define MUELU_CAN_USE_MIXED_PRECISION
#endif

namespace Thyra {

using Teuchos::ParameterList;
using Teuchos::RCP;
using Teuchos::rcp;
using Teuchos::rcp_const_cast;
using Teuchos::rcp_dynamic_cast;

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
bool Converters<Scalar, LocalOrdinal, GlobalOrdinal, Node>::replaceWithXpetra(ParameterList& paramList, std::string parameterName) {
  typedef typename Teuchos::ScalarTraits<Scalar>::magnitudeType Magnitude;
  typedef Xpetra::Operator<Scalar, LocalOrdinal, GlobalOrdinal, Node> XpOp;
  typedef Xpetra::ThyraUtils<Scalar, LocalOrdinal, GlobalOrdinal, Node> XpThyUtils;
  // typedef Xpetra::CrsMatrixWrap<Scalar,LocalOrdinal,GlobalOrdinal,Node>    XpCrsMatWrap;
  // typedef Xpetra::CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node>        XpCrsMat;
  typedef Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> XpMat;
  typedef Xpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node> XpMultVec;
  typedef Xpetra::MultiVector<Magnitude, LocalOrdinal, GlobalOrdinal, Node> XpMagMultVec;
  typedef Xpetra::Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node> XpVec;

  typedef Thyra::LinearOpBase<Scalar> ThyLinOpBase;
  typedef Thyra::DiagonalLinearOpBase<Scalar> ThyDiagLinOpBase;
  // typedef Thyra::XpetraLinearOp<Scalar, LocalOrdinal, GlobalOrdinal, Node> ThyXpOp;
  // typedef Thyra::SpmdVectorSpaceBase<Scalar>                               ThyVSBase;

  typedef Tpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> TpCrsMat;
  typedef Tpetra::Operator<Scalar, LocalOrdinal, GlobalOrdinal, Node> tOp;
  typedef Tpetra::Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node> tV;
  typedef Thyra::TpetraVector<Scalar, LocalOrdinal, GlobalOrdinal, Node> thyTpV;
  typedef Tpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node> tMV;
  typedef Tpetra::MultiVector<Magnitude, LocalOrdinal, GlobalOrdinal, Node> tMagMV;
#if defined(MUELU_CAN_USE_MIXED_PRECISION)
  typedef typename Teuchos::ScalarTraits<Magnitude>::halfPrecision HalfMagnitude;
  typedef Tpetra::MultiVector<HalfMagnitude, LocalOrdinal, GlobalOrdinal, Node> tHalfMagMV;
#endif

  if (paramList.isParameter(parameterName)) {
    if (paramList.isType<RCP<XpMat> >(parameterName))
      return true;
    else if (paramList.isType<RCP<const XpMat> >(parameterName)) {
      RCP<const XpMat> constM = paramList.get<RCP<const XpMat> >(parameterName);
      paramList.remove(parameterName);
      RCP<XpMat> M = rcp_const_cast<XpMat>(constM);
      paramList.set<RCP<XpMat> >(parameterName, M);
      return true;
    } else if (paramList.isType<RCP<XpMultVec> >(parameterName))
      return true;
    else if (paramList.isType<RCP<const XpMultVec> >(parameterName)) {
      RCP<const XpMultVec> constX = paramList.get<RCP<const XpMultVec> >(parameterName);
      paramList.remove(parameterName);
      RCP<XpMultVec> X = rcp_const_cast<XpMultVec>(constX);
      paramList.set<RCP<XpMultVec> >(parameterName, X);
      return true;
    } else if (paramList.isType<RCP<XpMagMultVec> >(parameterName))
      return true;
    else if (paramList.isType<RCP<const XpMagMultVec> >(parameterName)) {
      RCP<const XpMagMultVec> constX = paramList.get<RCP<const XpMagMultVec> >(parameterName);
      paramList.remove(parameterName);
      RCP<XpMagMultVec> X = rcp_const_cast<XpMagMultVec>(constX);
      paramList.set<RCP<XpMagMultVec> >(parameterName, X);
      return true;
    } else if (paramList.isType<RCP<TpCrsMat> >(parameterName)) {
      RCP<TpCrsMat> tM = paramList.get<RCP<TpCrsMat> >(parameterName);
      paramList.remove(parameterName);
      RCP<XpMat> xM = MueLu::TpetraCrs_To_XpetraMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>(tM);
      paramList.set<RCP<XpMat> >(parameterName, xM);
      return true;
    } else if (paramList.isType<RCP<tMV> >(parameterName)) {
      RCP<tMV> tpetra_X = paramList.get<RCP<tMV> >(parameterName);
      paramList.remove(parameterName);
      RCP<XpMultVec> X = MueLu::TpetraMultiVector_To_XpetraMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>(tpetra_X);
      paramList.set<RCP<XpMultVec> >(parameterName, X);
      TEUCHOS_TEST_FOR_EXCEPT(Teuchos::is_null(X));
      return true;
    } else if (paramList.isType<RCP<tMagMV> >(parameterName)) {
      RCP<tMagMV> tpetra_X = paramList.get<RCP<tMagMV> >(parameterName);
      paramList.remove(parameterName);
      RCP<XpMagMultVec> X = MueLu::TpetraMultiVector_To_XpetraMultiVector<Magnitude, LocalOrdinal, GlobalOrdinal, Node>(tpetra_X);
      paramList.set<RCP<XpMagMultVec> >(parameterName, X);
      TEUCHOS_TEST_FOR_EXCEPT(Teuchos::is_null(X));
      return true;
    }
#if defined(MUELU_CAN_USE_MIXED_PRECISION)
    else if (paramList.isType<RCP<tHalfMagMV> >(parameterName)) {
      RCP<tHalfMagMV> tpetra_hX = paramList.get<RCP<tHalfMagMV> >(parameterName);
      paramList.remove(parameterName);
      RCP<tMagMV> tpetra_X = rcp(new tMagMV(tpetra_hX->getMap(), tpetra_hX->getNumVectors()));
      Tpetra::deep_copy(*tpetra_X, *tpetra_hX);
      RCP<XpMagMultVec> X = MueLu::TpetraMultiVector_To_XpetraMultiVector<Magnitude, LocalOrdinal, GlobalOrdinal, Node>(tpetra_X);
      paramList.set<RCP<XpMagMultVec> >(parameterName, X);
      TEUCHOS_TEST_FOR_EXCEPT(Teuchos::is_null(X));
      return true;
    }
#endif
    else if (paramList.isType<RCP<const ThyDiagLinOpBase> >(parameterName) ||
             (paramList.isType<RCP<const ThyLinOpBase> >(parameterName) && !rcp_dynamic_cast<const ThyDiagLinOpBase>(paramList.get<RCP<const ThyLinOpBase> >(parameterName)).is_null())) {
      RCP<const ThyDiagLinOpBase> thyM;
      if (paramList.isType<RCP<const ThyDiagLinOpBase> >(parameterName))
        thyM = paramList.get<RCP<const ThyDiagLinOpBase> >(parameterName);
      else
        thyM = rcp_dynamic_cast<const ThyDiagLinOpBase>(paramList.get<RCP<const ThyLinOpBase> >(parameterName), true);
      paramList.remove(parameterName);
      RCP<const Thyra::VectorBase<Scalar> > diag = thyM->getDiag();

      RCP<const XpVec> xpDiag;
      if (!rcp_dynamic_cast<const thyTpV>(diag).is_null()) {
        RCP<const tV> tDiag = Thyra::TpetraOperatorVectorExtraction<Scalar, LocalOrdinal, GlobalOrdinal, Node>::getConstTpetraVector(diag);
        if (!tDiag.is_null())
          xpDiag = Xpetra::toXpetra(tDiag);
      }
      TEUCHOS_ASSERT(!xpDiag.is_null());
      RCP<XpMat> M = Xpetra::MatrixFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Build(xpDiag);
      paramList.set<RCP<XpMat> >(parameterName, M);
      return true;
    } else if (paramList.isType<RCP<const ThyLinOpBase> >(parameterName)) {
      RCP<const ThyLinOpBase> thyM = paramList.get<RCP<const ThyLinOpBase> >(parameterName);
      paramList.remove(parameterName);
      try {
        RCP<XpMat> M = XpThyUtils::toXpetra(Teuchos::rcp_const_cast<ThyLinOpBase>(thyM));
        paramList.set<RCP<XpMat> >(parameterName, M);
      } catch (std::exception& e) {
        RCP<XpOp> M                                                                  = XpThyUtils::toXpetraOperator(Teuchos::rcp_const_cast<ThyLinOpBase>(thyM));
        RCP<Xpetra::TpetraOperator<Scalar, LocalOrdinal, GlobalOrdinal, Node> > tpOp = rcp_dynamic_cast<Xpetra::TpetraOperator<Scalar, LocalOrdinal, GlobalOrdinal, Node> >(M, true);
        RCP<tOp> tO                                                                  = tpOp->getOperator();
        RCP<tV> diag;
        if (tO->hasDiagonal()) {
          diag = rcp(new tV(tO->getRangeMap()));
          tO->getLocalDiagCopy(*diag);
        }
        auto fTpRow                                                                   = rcp(new MueLu::TpetraOperatorAsRowMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>(tO, diag));
        RCP<Xpetra::TpetraOperator<Scalar, LocalOrdinal, GlobalOrdinal, Node> > tpFOp = rcp(new Xpetra::TpetraOperator<Scalar, LocalOrdinal, GlobalOrdinal, Node>(fTpRow));
        auto op                                                                       = rcp_dynamic_cast<XpOp>(tpFOp);
        paramList.set<RCP<XpOp> >(parameterName, op);
      }
      return true;
    } else {
      TEUCHOS_TEST_FOR_EXCEPTION(true, MueLu::Exceptions::RuntimeError, "Parameter " << parameterName << " has wrong type.");
      return false;
    }
  } else
    return false;
}

#ifdef HAVE_MUELU_EPETRA
template <class GlobalOrdinal>
bool Converters<double, int, GlobalOrdinal, Tpetra::KokkosCompat::KokkosSerialWrapperNode>::replaceWithXpetra(ParameterList& paramList, std::string parameterName) {
  typedef double Scalar;
  typedef int LocalOrdinal;
  typedef Tpetra::KokkosCompat::KokkosSerialWrapperNode Node;
  typedef typename Teuchos::ScalarTraits<Scalar>::magnitudeType Magnitude;
  typedef Xpetra::Operator<Scalar, LocalOrdinal, GlobalOrdinal, Node> XpOp;
  typedef Xpetra::ThyraUtils<Scalar, LocalOrdinal, GlobalOrdinal, Node> XpThyUtils;
  typedef Xpetra::CrsMatrixWrap<Scalar, LocalOrdinal, GlobalOrdinal, Node> XpCrsMatWrap;
  typedef Xpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> XpCrsMat;
  typedef Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> XpMat;
  typedef Xpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node> XpMultVec;
  typedef Xpetra::MultiVector<Magnitude, LocalOrdinal, GlobalOrdinal, Node> XpMagMultVec;
  typedef Xpetra::Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node> XpVec;

  typedef Thyra::LinearOpBase<Scalar> ThyLinOpBase;
  typedef Thyra::DiagonalLinearOpBase<Scalar> ThyDiagLinOpBase;
  typedef Thyra::SpmdVectorSpaceBase<Scalar> ThyVSBase;

  typedef Tpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> TpCrsMat;
  typedef Tpetra::Operator<Scalar, LocalOrdinal, GlobalOrdinal, Node> tOp;
  typedef Tpetra::Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node> tV;
  typedef Thyra::TpetraVector<Scalar, LocalOrdinal, GlobalOrdinal, Node> thyTpV;
  typedef Tpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node> tMV;
  typedef Tpetra::MultiVector<Magnitude, LocalOrdinal, GlobalOrdinal, Node> tMagMV;
#if defined(MUELU_CAN_USE_MIXED_PRECISION)
  typedef typename Teuchos::ScalarTraits<Magnitude>::halfPrecision HalfMagnitude;
  typedef Tpetra::MultiVector<HalfMagnitude, LocalOrdinal, GlobalOrdinal, Node> tHalfMagMV;
#endif
#if defined(HAVE_MUELU_EPETRA)
  typedef Xpetra::EpetraCrsMatrixT<GlobalOrdinal, Node> XpEpCrsMat;
#endif

  if (paramList.isParameter(parameterName)) {
    if (paramList.isType<RCP<XpMat> >(parameterName))
      return true;
    else if (paramList.isType<RCP<const XpMat> >(parameterName)) {
      RCP<const XpMat> constM = paramList.get<RCP<const XpMat> >(parameterName);
      paramList.remove(parameterName);
      RCP<XpMat> M = rcp_const_cast<XpMat>(constM);
      paramList.set<RCP<XpMat> >(parameterName, M);
      return true;
    } else if (paramList.isType<RCP<XpMultVec> >(parameterName))
      return true;
    else if (paramList.isType<RCP<const XpMultVec> >(parameterName)) {
      RCP<const XpMultVec> constX = paramList.get<RCP<const XpMultVec> >(parameterName);
      paramList.remove(parameterName);
      RCP<XpMultVec> X = rcp_const_cast<XpMultVec>(constX);
      paramList.set<RCP<XpMultVec> >(parameterName, X);
      return true;
    } else if (paramList.isType<RCP<XpMagMultVec> >(parameterName))
      return true;
    else if (paramList.isType<RCP<const XpMagMultVec> >(parameterName)) {
      RCP<const XpMagMultVec> constX = paramList.get<RCP<const XpMagMultVec> >(parameterName);
      paramList.remove(parameterName);
      RCP<XpMagMultVec> X = rcp_const_cast<XpMagMultVec>(constX);
      paramList.set<RCP<XpMagMultVec> >(parameterName, X);
      return true;
    } else if (paramList.isType<RCP<TpCrsMat> >(parameterName)) {
      RCP<TpCrsMat> tM = paramList.get<RCP<TpCrsMat> >(parameterName);
      paramList.remove(parameterName);
      RCP<XpMat> xM = MueLu::TpetraCrs_To_XpetraMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>(tM);
      paramList.set<RCP<XpMat> >(parameterName, xM);
      return true;
    } else if (paramList.isType<RCP<tMV> >(parameterName)) {
      RCP<tMV> tpetra_X = paramList.get<RCP<tMV> >(parameterName);
      paramList.remove(parameterName);
      RCP<XpMultVec> X = MueLu::TpetraMultiVector_To_XpetraMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>(tpetra_X);
      paramList.set<RCP<XpMultVec> >(parameterName, X);
      TEUCHOS_TEST_FOR_EXCEPT(Teuchos::is_null(X));
      return true;
    } else if (paramList.isType<RCP<tMagMV> >(parameterName)) {
      RCP<tMagMV> tpetra_X = paramList.get<RCP<tMagMV> >(parameterName);
      paramList.remove(parameterName);
      RCP<XpMagMultVec> X = MueLu::TpetraMultiVector_To_XpetraMultiVector<Magnitude, LocalOrdinal, GlobalOrdinal, Node>(tpetra_X);
      paramList.set<RCP<XpMagMultVec> >(parameterName, X);
      TEUCHOS_TEST_FOR_EXCEPT(Teuchos::is_null(X));
      return true;
    }
#if defined(MUELU_CAN_USE_MIXED_PRECISION)
    else if (paramList.isType<RCP<tHalfMagMV> >(parameterName)) {
      RCP<tHalfMagMV> tpetra_hX = paramList.get<RCP<tHalfMagMV> >(parameterName);
      paramList.remove(parameterName);
      RCP<tMagMV> tpetra_X = rcp(new tMagMV(tpetra_hX->getMap(), tpetra_hX->getNumVectors()));
      Tpetra::deep_copy(*tpetra_X, *tpetra_hX);
      RCP<XpMagMultVec> X = MueLu::TpetraMultiVector_To_XpetraMultiVector<Magnitude, LocalOrdinal, GlobalOrdinal, Node>(tpetra_X);
      paramList.set<RCP<XpMagMultVec> >(parameterName, X);
      TEUCHOS_TEST_FOR_EXCEPT(Teuchos::is_null(X));
      return true;
    }
#endif
#ifdef HAVE_MUELU_EPETRA
    else if (paramList.isType<RCP<Epetra_CrsMatrix> >(parameterName)) {
      RCP<Epetra_CrsMatrix> eM = paramList.get<RCP<Epetra_CrsMatrix> >(parameterName);
      paramList.remove(parameterName);
      RCP<XpEpCrsMat> xeM   = rcp(new XpEpCrsMat(eM));
      RCP<XpCrsMat> xCrsM   = rcp_dynamic_cast<XpCrsMat>(xeM, true);
      RCP<XpCrsMatWrap> xwM = rcp(new XpCrsMatWrap(xCrsM));
      RCP<XpMat> xM         = rcp_dynamic_cast<XpMat>(xwM);
      paramList.set<RCP<XpMat> >(parameterName, xM);
      return true;
    } else if (paramList.isType<RCP<Epetra_MultiVector> >(parameterName)) {
      RCP<Epetra_MultiVector> epetra_X = Teuchos::null;
      epetra_X                         = paramList.get<RCP<Epetra_MultiVector> >(parameterName);
      paramList.remove(parameterName);
      RCP<Xpetra::EpetraMultiVectorT<int, Node> > xpEpX           = rcp(new Xpetra::EpetraMultiVectorT<int, Node>(epetra_X));
      RCP<Xpetra::MultiVector<Scalar, int, int, Node> > xpEpXMult = rcp_dynamic_cast<Xpetra::MultiVector<Scalar, int, int, Node> >(xpEpX, true);
      RCP<XpMultVec> X                                            = rcp_dynamic_cast<XpMultVec>(xpEpXMult, true);
      paramList.set<RCP<XpMultVec> >(parameterName, X);
      return true;
    }
#endif
    else if (paramList.isType<RCP<const ThyDiagLinOpBase> >(parameterName) ||
             (paramList.isType<RCP<const ThyLinOpBase> >(parameterName) && !rcp_dynamic_cast<const ThyDiagLinOpBase>(paramList.get<RCP<const ThyLinOpBase> >(parameterName)).is_null())) {
      RCP<const ThyDiagLinOpBase> thyM;
      if (paramList.isType<RCP<const ThyDiagLinOpBase> >(parameterName))
        thyM = paramList.get<RCP<const ThyDiagLinOpBase> >(parameterName);
      else
        thyM = rcp_dynamic_cast<const ThyDiagLinOpBase>(paramList.get<RCP<const ThyLinOpBase> >(parameterName), true);
      paramList.remove(parameterName);
      RCP<const Thyra::VectorBase<Scalar> > diag = thyM->getDiag();

      RCP<const XpVec> xpDiag;
      if (!rcp_dynamic_cast<const thyTpV>(diag).is_null()) {
        RCP<const tV> tDiag = Thyra::TpetraOperatorVectorExtraction<Scalar, LocalOrdinal, GlobalOrdinal, Node>::getConstTpetraVector(diag);
        if (!tDiag.is_null())
          xpDiag = Xpetra::toXpetra(tDiag);
      }
#ifdef HAVE_MUELU_EPETRA
      if (xpDiag.is_null()) {
        RCP<const Epetra_Comm> comm = Thyra::get_Epetra_Comm(*rcp_dynamic_cast<const ThyVSBase>(thyM->range())->getComm());
        RCP<const Epetra_Map> map   = Thyra::get_Epetra_Map(*(thyM->range()), comm);
        if (!map.is_null()) {
          RCP<const Epetra_Vector> eDiag                  = Thyra::get_Epetra_Vector(*map, diag);
          RCP<Epetra_Vector> nceDiag                      = rcp_const_cast<Epetra_Vector>(eDiag);
          RCP<Xpetra::EpetraVectorT<int, Node> > xpEpDiag = rcp(new Xpetra::EpetraVectorT<int, Node>(nceDiag));
          xpDiag                                          = rcp_dynamic_cast<XpVec>(xpEpDiag, true);
        }
      }
#endif
      TEUCHOS_ASSERT(!xpDiag.is_null());
      RCP<XpMat> M = Xpetra::MatrixFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Build(xpDiag);
      paramList.set<RCP<XpMat> >(parameterName, M);
      return true;
    } else if (paramList.isType<RCP<const ThyLinOpBase> >(parameterName)) {
      RCP<const ThyLinOpBase> thyM = paramList.get<RCP<const ThyLinOpBase> >(parameterName);
      paramList.remove(parameterName);
      try {
        RCP<XpMat> M = XpThyUtils::toXpetra(Teuchos::rcp_const_cast<ThyLinOpBase>(thyM));
        paramList.set<RCP<XpMat> >(parameterName, M);
      } catch (std::exception& e) {
        RCP<XpOp> M                                                                  = XpThyUtils::toXpetraOperator(Teuchos::rcp_const_cast<ThyLinOpBase>(thyM));
        RCP<Xpetra::TpetraOperator<Scalar, LocalOrdinal, GlobalOrdinal, Node> > tpOp = rcp_dynamic_cast<Xpetra::TpetraOperator<Scalar, LocalOrdinal, GlobalOrdinal, Node> >(M, true);
        RCP<tOp> tO                                                                  = tpOp->getOperator();
        RCP<tV> diag;
        if (tO->hasDiagonal()) {
          diag = rcp(new tV(tO->getRangeMap()));
          tO->getLocalDiagCopy(*diag);
        }
        auto fTpRow                                                                   = rcp(new MueLu::TpetraOperatorAsRowMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>(tO, diag));
        RCP<Xpetra::TpetraOperator<Scalar, LocalOrdinal, GlobalOrdinal, Node> > tpFOp = rcp(new Xpetra::TpetraOperator<Scalar, LocalOrdinal, GlobalOrdinal, Node>(fTpRow));
        auto op                                                                       = rcp_dynamic_cast<XpOp>(tpFOp);
        paramList.set<RCP<XpOp> >(parameterName, op);
      }
      return true;
    } else {
      TEUCHOS_TEST_FOR_EXCEPTION(true, MueLu::Exceptions::RuntimeError, "Parameter " << parameterName << " has wrong type.");
      return false;
    }
  } else
    return false;
}
#endif

// Constructors/initializers/accessors

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
MueLuPreconditionerFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::MueLuPreconditionerFactory()
  : paramList_(rcp(new ParameterList())) {}

// Overridden from PreconditionerFactoryBase

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
bool MueLuPreconditionerFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::isCompatible(const LinearOpSourceBase<Scalar>& fwdOpSrc) const {
  const RCP<const LinearOpBase<Scalar> > fwdOp = fwdOpSrc.getOp();

  if (Xpetra::ThyraUtils<Scalar, LocalOrdinal, GlobalOrdinal, Node>::isTpetra(fwdOp)) return true;

#ifdef HAVE_MUELU_EPETRA
  if (Xpetra::ThyraUtils<Scalar, LocalOrdinal, GlobalOrdinal, Node>::isEpetra(fwdOp)) return true;
#endif

  if (Xpetra::ThyraUtils<Scalar, LocalOrdinal, GlobalOrdinal, Node>::isBlockedOperator(fwdOp)) return true;

  return false;
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
RCP<PreconditionerBase<Scalar> > MueLuPreconditionerFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::createPrec() const {
  return Teuchos::rcp(new DefaultPreconditioner<Scalar>);
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void MueLuPreconditionerFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    initializePrec(const RCP<const LinearOpSourceBase<Scalar> >& fwdOpSrc, PreconditionerBase<Scalar>* prec, const ESupportSolveUse supportSolveUse) const {
  Teuchos::TimeMonitor tM(*Teuchos::TimeMonitor::getNewTimer(std::string("ThyraMueLu::initializePrec")));

  using Teuchos::rcp_dynamic_cast;

  // we are using typedefs here, since we are using objects from different packages (Xpetra, Thyra,...)
  typedef Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node> XpMap;
  typedef Xpetra::Operator<Scalar, LocalOrdinal, GlobalOrdinal, Node> XpOp;
  typedef MueLu::XpetraOperator<Scalar, LocalOrdinal, GlobalOrdinal, Node> MueLuXpOp;
  typedef Xpetra::ThyraUtils<Scalar, LocalOrdinal, GlobalOrdinal, Node> XpThyUtils;
  // typedef Xpetra::CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node>        XpCrsMat;
  typedef Xpetra::BlockedCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> XpBlockedCrsMat;
  typedef Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> XpMat;
  // typedef Xpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>      XpMultVec;
  // typedef Xpetra::MultiVector<typename Teuchos::ScalarTraits<Scalar>::coordinateType,LocalOrdinal,GlobalOrdinal,Node>      XpMultVecDouble;
  typedef Thyra::LinearOpBase<Scalar> ThyLinOpBase;
  typedef Thyra::XpetraLinearOp<Scalar, LocalOrdinal, GlobalOrdinal, Node> ThyXpOp;
  typedef Xpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node> XpMV;
  typedef typename Teuchos::ScalarTraits<Scalar>::magnitudeType Magnitude;
  typedef Xpetra::MultiVector<Magnitude, LocalOrdinal, GlobalOrdinal, Node> XpmMV;
#if defined(MUELU_CAN_USE_MIXED_PRECISION)
  typedef Xpetra::TpetraHalfPrecisionOperator<Scalar, LocalOrdinal, GlobalOrdinal, Node> XpHalfPrecOp;
  typedef typename XpHalfPrecOp::HalfScalar HalfScalar;
  typedef Xpetra::Operator<HalfScalar, LocalOrdinal, GlobalOrdinal, Node> XpHalfOp;
  typedef MueLu::XpetraOperator<HalfScalar, LocalOrdinal, GlobalOrdinal, Node> MueLuHalfXpOp;
  typedef typename Teuchos::ScalarTraits<Magnitude>::halfPrecision HalfMagnitude;
  typedef Xpetra::MultiVector<HalfScalar, LocalOrdinal, GlobalOrdinal, Node> XphMV;
  typedef Xpetra::MultiVector<HalfMagnitude, LocalOrdinal, GlobalOrdinal, Node> XphmMV;
  typedef Xpetra::Matrix<HalfScalar, LocalOrdinal, GlobalOrdinal, Node> XphMat;
#endif

  // Check precondition
  TEUCHOS_ASSERT(Teuchos::nonnull(fwdOpSrc));
  TEUCHOS_ASSERT(this->isCompatible(*fwdOpSrc));
  TEUCHOS_ASSERT(prec);

  // Create a copy, as we may remove some things from the list
  ParameterList paramList = *paramList_;

  // Retrieve wrapped concrete Xpetra matrix from FwdOp
  const RCP<const ThyLinOpBase> fwdOp = fwdOpSrc->getOp();
  TEUCHOS_TEST_FOR_EXCEPT(Teuchos::is_null(fwdOp));

  // Check whether it is Epetra/Tpetra
  bool bIsEpetra  = XpThyUtils::isEpetra(fwdOp);
  bool bIsTpetra  = XpThyUtils::isTpetra(fwdOp);
  bool bIsBlocked = XpThyUtils::isBlockedOperator(fwdOp);
  TEUCHOS_TEST_FOR_EXCEPT((bIsEpetra == true && bIsTpetra == true));
  TEUCHOS_TEST_FOR_EXCEPT((bIsEpetra == bIsTpetra) && bIsBlocked == false);
  TEUCHOS_TEST_FOR_EXCEPT((bIsEpetra != bIsTpetra) && bIsBlocked == true);

  RCP<XpMat> A = Teuchos::null;
  if (bIsBlocked) {
    Teuchos::RCP<const Thyra::BlockedLinearOpBase<Scalar> > ThyBlockedOp =
        Teuchos::rcp_dynamic_cast<const Thyra::BlockedLinearOpBase<Scalar> >(fwdOp);
    TEUCHOS_TEST_FOR_EXCEPT(Teuchos::is_null(ThyBlockedOp));

    TEUCHOS_TEST_FOR_EXCEPT(ThyBlockedOp->blockExists(0, 0) == false);

    Teuchos::RCP<const LinearOpBase<Scalar> > b00 = ThyBlockedOp->getBlock(0, 0);
    TEUCHOS_TEST_FOR_EXCEPT(Teuchos::is_null(b00));

    // wrap the forward operator as an Xpetra::Matrix that MueLu can work with
    // MueLu needs a non-const object as input
    RCP<XpMat> A00 = XpThyUtils::toXpetra(Teuchos::rcp_const_cast<LinearOpBase<Scalar> >(b00));
    TEUCHOS_TEST_FOR_EXCEPT(Teuchos::is_null(A00));

    RCP<const XpMap> rowmap00           = A00->getRowMap();
    RCP<const Teuchos::Comm<int> > comm = rowmap00->getComm();

    // create a Xpetra::BlockedCrsMatrix which derives from Xpetra::Matrix that MueLu can work with
    RCP<XpBlockedCrsMat> bMat = Teuchos::rcp(new XpBlockedCrsMat(ThyBlockedOp, comm));
    TEUCHOS_TEST_FOR_EXCEPT(Teuchos::is_null(bMat));

    // save blocked matrix
    A = bMat;
  } else {
    // wrap the forward operator as an Xpetra::Matrix that MueLu can work with
    // MueLu needs a non-const object as input
    A = XpThyUtils::toXpetra(Teuchos::rcp_const_cast<ThyLinOpBase>(fwdOp));
  }
  TEUCHOS_TEST_FOR_EXCEPT(Teuchos::is_null(A));

  // Retrieve concrete preconditioner object
  const Teuchos::Ptr<DefaultPreconditioner<Scalar> > defaultPrec = Teuchos::ptr(dynamic_cast<DefaultPreconditioner<Scalar>*>(prec));
  TEUCHOS_TEST_FOR_EXCEPT(Teuchos::is_null(defaultPrec));

  // extract preconditioner operator
  RCP<ThyLinOpBase> thyra_precOp = Teuchos::null;
  thyra_precOp                   = rcp_dynamic_cast<Thyra::LinearOpBase<Scalar> >(defaultPrec->getNonconstUnspecifiedPrecOp(), true);

  // make a decision whether to (re)build the multigrid preconditioner or reuse the old one
  // rebuild preconditioner if startingOver == true
  // reuse preconditioner if startingOver == false
  const bool startingOver = (thyra_precOp.is_null() || !paramList.isParameter("reuse: type") || paramList.get<std::string>("reuse: type") == "none");
  bool useHalfPrecision   = false;
  if (paramList.isParameter("half precision"))
    useHalfPrecision = paramList.get<bool>("half precision");
  else if (paramList.isSublist("Hierarchy") && paramList.sublist("Hierarchy").isParameter("half precision"))
    useHalfPrecision = paramList.sublist("Hierarchy").get<bool>("half precision");
  if (useHalfPrecision)
    TEUCHOS_TEST_FOR_EXCEPTION(!bIsTpetra, MueLu::Exceptions::RuntimeError, "The only scalar type Epetra knows is double, so a half precision preconditioner cannot be constructed.");

  RCP<XpOp> xpPrecOp;
  if (startingOver == true) {
    // Convert to Xpetra
    std::list<std::string> convertXpetra = {"Coordinates", "Nullspace"};
    for (auto it = convertXpetra.begin(); it != convertXpetra.end(); ++it)
      Converters<Scalar, LocalOrdinal, GlobalOrdinal, Node>::replaceWithXpetra(paramList, *it);

    for (int lvlNo = 0; lvlNo < 10; ++lvlNo) {
      if (paramList.isSublist("level " + std::to_string(lvlNo) + " user data")) {
        ParameterList& lvlList = paramList.sublist("level " + std::to_string(lvlNo) + " user data");
        std::list<std::string> convertKeys;
        for (auto it = lvlList.begin(); it != lvlList.end(); ++it)
          convertKeys.push_back(lvlList.name(it));
        for (auto it = convertKeys.begin(); it != convertKeys.end(); ++it)
          Converters<Scalar, LocalOrdinal, GlobalOrdinal, Node>::replaceWithXpetra(lvlList, *it);
      }
    }

    if (useHalfPrecision) {
#if defined(MUELU_CAN_USE_MIXED_PRECISION)

      // CAG: There is nothing special about the combination double-float,
      //      except that I feel somewhat confident that Trilinos builds
      //      with both scalar types.

      // convert to half precision
      RCP<XphMat> halfA                     = Xpetra::convertToHalfPrecision(A);
      const std::string userName            = "user data";
      Teuchos::ParameterList& userParamList = paramList.sublist(userName);
      if (userParamList.isType<RCP<XpmMV> >("Coordinates")) {
        RCP<XpmMV> coords = userParamList.get<RCP<XpmMV> >("Coordinates");
        userParamList.remove("Coordinates");
        RCP<XphmMV> halfCoords = Xpetra::convertToHalfPrecision(coords);
        userParamList.set("Coordinates", halfCoords);
      }
      if (userParamList.isType<RCP<XpMV> >("Nullspace")) {
        RCP<XpMV> nullspace = userParamList.get<RCP<XpMV> >("Nullspace");
        userParamList.remove("Nullspace");
        RCP<XphMV> halfNullspace = Xpetra::convertToHalfPrecision(nullspace);
        userParamList.set("Nullspace", halfNullspace);
      }
      if (paramList.isType<RCP<XpmMV> >("Coordinates")) {
        RCP<XpmMV> coords = paramList.get<RCP<XpmMV> >("Coordinates");
        paramList.remove("Coordinates");
        RCP<XphmMV> halfCoords = Xpetra::convertToHalfPrecision(coords);
        userParamList.set("Coordinates", halfCoords);
      }
      if (paramList.isType<RCP<XpMV> >("Nullspace")) {
        RCP<XpMV> nullspace = paramList.get<RCP<XpMV> >("Nullspace");
        paramList.remove("Nullspace");
        RCP<XphMV> halfNullspace = Xpetra::convertToHalfPrecision(nullspace);
        userParamList.set("Nullspace", halfNullspace);
      }

      // build a new half-precision MueLu preconditioner

      RCP<MueLu::Hierarchy<HalfScalar, LocalOrdinal, GlobalOrdinal, Node> > H = MueLu::CreateXpetraPreconditioner(halfA, paramList);
      RCP<MueLuHalfXpOp> xpOp                                                 = rcp(new MueLuHalfXpOp(H));
      xpPrecOp                                                                = rcp(new XpHalfPrecOp(xpOp));
#else
      TEUCHOS_TEST_FOR_EXCEPT(true);
#endif
    } else {
      const std::string userName            = "user data";
      Teuchos::ParameterList& userParamList = paramList.sublist(userName);
      if (paramList.isType<RCP<XpmMV> >("Coordinates")) {
        RCP<XpmMV> coords = paramList.get<RCP<XpmMV> >("Coordinates");
        paramList.remove("Coordinates");
        userParamList.set("Coordinates", coords);
      }
      if (paramList.isType<RCP<XpMV> >("Nullspace")) {
        RCP<XpMV> nullspace = paramList.get<RCP<XpMV> >("Nullspace");
        paramList.remove("Nullspace");
        userParamList.set("Nullspace", nullspace);
      }

      // build a new MueLu RefMaxwell preconditioner
      RCP<MueLu::Hierarchy<Scalar, LocalOrdinal, GlobalOrdinal, Node> > H = MueLu::CreateXpetraPreconditioner(A, paramList);
      xpPrecOp                                                            = rcp(new MueLuXpOp(H));
    }
  } else {
    // reuse old MueLu hierarchy stored in MueLu Xpetra operator and put in new matrix
    RCP<ThyXpOp> thyXpOp = rcp_dynamic_cast<ThyXpOp>(thyra_precOp, true);
    xpPrecOp             = rcp_dynamic_cast<XpOp>(thyXpOp->getXpetraOperator(), true);
#if defined(MUELU_CAN_USE_MIXED_PRECISION)
    RCP<XpHalfPrecOp> xpHalfPrecOp = rcp_dynamic_cast<XpHalfPrecOp>(xpPrecOp);
    if (!xpHalfPrecOp.is_null()) {
      RCP<MueLu::Hierarchy<HalfScalar, LocalOrdinal, GlobalOrdinal, Node> > H = rcp_dynamic_cast<MueLuHalfXpOp>(xpHalfPrecOp->GetHalfPrecisionOperator(), true)->GetHierarchy();
      RCP<XphMat> halfA                                                       = Xpetra::convertToHalfPrecision(A);

      TEUCHOS_TEST_FOR_EXCEPTION(!H->GetNumLevels(), MueLu::Exceptions::RuntimeError,
                                 "Thyra::MueLuPreconditionerFactory: Hierarchy has no levels in it");
      TEUCHOS_TEST_FOR_EXCEPTION(!H->GetLevel(0)->IsAvailable("A"), MueLu::Exceptions::RuntimeError,
                                 "Thyra::MueLuPreconditionerFactory: Hierarchy has no fine level operator");
      RCP<MueLu::Level> level0 = H->GetLevel(0);
      RCP<XpHalfOp> O0         = level0->Get<RCP<XpHalfOp> >("A");
      RCP<XphMat> A0           = rcp_dynamic_cast<XphMat>(O0, true);

      if (!A0.is_null()) {
        // If a user provided a "number of equations" argument in a parameter list
        // during the initial setup, we must honor that settings and reuse it for
        // all consequent setups.
        halfA->SetFixedBlockSize(A0->GetFixedBlockSize());
      }

      // set new matrix
      level0->Set("A", halfA);

      H->SetupRe();
    } else
#endif
    {
      // get old MueLu hierarchy
      RCP<MueLuXpOp> xpOp                                                 = rcp_dynamic_cast<MueLuXpOp>(thyXpOp->getXpetraOperator(), true);
      RCP<MueLu::Hierarchy<Scalar, LocalOrdinal, GlobalOrdinal, Node> > H = xpOp->GetHierarchy();
      ;

      TEUCHOS_TEST_FOR_EXCEPTION(!H->GetNumLevels(), MueLu::Exceptions::RuntimeError,
                                 "Thyra::MueLuPreconditionerFactory: Hierarchy has no levels in it");
      TEUCHOS_TEST_FOR_EXCEPTION(!H->GetLevel(0)->IsAvailable("A"), MueLu::Exceptions::RuntimeError,
                                 "Thyra::MueLuPreconditionerFactory: Hierarchy has no fine level operator");
      RCP<MueLu::Level> level0 = H->GetLevel(0);
      RCP<XpOp> O0             = level0->Get<RCP<XpOp> >("A");
      RCP<XpMat> A0            = rcp_dynamic_cast<XpMat>(O0);

      if (!A0.is_null()) {
        // If a user provided a "number of equations" argument in a parameter list
        // during the initial setup, we must honor that settings and reuse it for
        // all consequent setups.
        A->SetFixedBlockSize(A0->GetFixedBlockSize());
      }

      // set new matrix
      level0->Set("A", A);

      H->SetupRe();
    }
  }

  // wrap preconditioner in thyraPrecOp
  RCP<const VectorSpaceBase<Scalar> > thyraRangeSpace  = Xpetra::ThyraUtils<Scalar, LocalOrdinal, GlobalOrdinal, Node>::toThyra(xpPrecOp->getRangeMap());
  RCP<const VectorSpaceBase<Scalar> > thyraDomainSpace = Xpetra::ThyraUtils<Scalar, LocalOrdinal, GlobalOrdinal, Node>::toThyra(xpPrecOp->getDomainMap());

  RCP<ThyLinOpBase> thyraPrecOp = Thyra::xpetraLinearOp<Scalar, LocalOrdinal, GlobalOrdinal, Node>(thyraRangeSpace, thyraDomainSpace, xpPrecOp);
  TEUCHOS_TEST_FOR_EXCEPT(Teuchos::is_null(thyraPrecOp));

  defaultPrec->initializeUnspecified(thyraPrecOp);
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void MueLuPreconditionerFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    uninitializePrec(PreconditionerBase<Scalar>* prec, RCP<const LinearOpSourceBase<Scalar> >* fwdOp, ESupportSolveUse* supportSolveUse) const {
  TEUCHOS_ASSERT(prec);

  // Retrieve concrete preconditioner object
  const Teuchos::Ptr<DefaultPreconditioner<Scalar> > defaultPrec = Teuchos::ptr(dynamic_cast<DefaultPreconditioner<Scalar>*>(prec));
  TEUCHOS_TEST_FOR_EXCEPT(Teuchos::is_null(defaultPrec));

  if (fwdOp) {
    // TODO: Implement properly instead of returning default value
    *fwdOp = Teuchos::null;
  }

  if (supportSolveUse) {
    // TODO: Implement properly instead of returning default value
    *supportSolveUse = Thyra::SUPPORT_SOLVE_UNSPECIFIED;
  }

  defaultPrec->uninitialize();
}

// Overridden from ParameterListAcceptor
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void MueLuPreconditionerFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::setParameterList(RCP<ParameterList> const& paramList) {
  TEUCHOS_TEST_FOR_EXCEPT(Teuchos::is_null(paramList));
  paramList_ = paramList;
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
RCP<ParameterList> MueLuPreconditionerFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::getNonconstParameterList() {
  return paramList_;
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
RCP<ParameterList> MueLuPreconditionerFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::unsetParameterList() {
  RCP<ParameterList> savedParamList = paramList_;
  paramList_                        = Teuchos::null;
  return savedParamList;
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
RCP<const ParameterList> MueLuPreconditionerFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::getParameterList() const {
  return paramList_;
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
RCP<const ParameterList> MueLuPreconditionerFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::getValidParameters() const {
  static RCP<const ParameterList> validPL;

  if (Teuchos::is_null(validPL))
    validPL = rcp(new ParameterList());

  return validPL;
}

// Public functions overridden from Teuchos::Describable
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
std::string MueLuPreconditionerFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::description() const {
  return "Thyra::MueLuPreconditionerFactory";
}
}  // namespace Thyra

#endif  // HAVE_MUELU_STRATIMIKOS

#endif  // ifdef THYRA_MUELU_PRECONDITIONER_FACTORY_DEF_HPP
