// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef _MiniEM_Utils_cpp_
#define _MiniEM_Utils_cpp_

#include "Teuchos_CompilerCodeTweakMacros.hpp"
#include "MiniEM_Utils.hpp"
#include "Thyra_DiagonalLinearOpBase.hpp"
#include "Thyra_DefaultDiagonalLinearOp.hpp"

namespace mini_em {

  void writeOut(const std::string & s,const Thyra::LinearOpBase<double> & op)
  {
    using Teuchos::RCP;
    using Teuchos::rcp_dynamic_cast;
    using NT = panzer::TpetraNodeType;
    const RCP<const Thyra::TpetraLinearOp<double,int,panzer::GlobalOrdinal,NT> > tOp = rcp_dynamic_cast<const Thyra::TpetraLinearOp<double,int,panzer::GlobalOrdinal,NT> >(Teuchos::rcpFromRef(op));

    if(tOp != Teuchos::null) {
      *Teko::getOutputStream() << "Dumping matrix \'" << s << "\'" << std::endl;
      const RCP<const Tpetra::CrsMatrix<double,int,panzer::GlobalOrdinal,NT> > crsOp = rcp_dynamic_cast<const Tpetra::CrsMatrix<double,int,panzer::GlobalOrdinal,NT> >(tOp->getConstTpetraOperator());
      if (crsOp != Teuchos::null) {
        Tpetra::MatrixMarket::Writer<Tpetra::CrsMatrix<double,int,panzer::GlobalOrdinal,NT> >::writeMapFile(("rowmap_"+s).c_str(),*(crsOp->getRowMap()));
        Tpetra::MatrixMarket::Writer<Tpetra::CrsMatrix<double,int,panzer::GlobalOrdinal,NT> >::writeMapFile(("colmap_"+s).c_str(),*(crsOp->getColMap()));
        Tpetra::MatrixMarket::Writer<Tpetra::CrsMatrix<double,int,panzer::GlobalOrdinal,NT> >::writeMapFile(("domainmap_"+s).c_str(),*(crsOp->getDomainMap()));
        Tpetra::MatrixMarket::Writer<Tpetra::CrsMatrix<double,int,panzer::GlobalOrdinal,NT> >::writeMapFile(("rangemap_"+s).c_str(),*(crsOp->getRangeMap()));
        Tpetra::MatrixMarket::Writer<Tpetra::CrsMatrix<double,int,panzer::GlobalOrdinal,NT> >::writeSparseFile(s.c_str(),crsOp);
      } else {
        auto tO = tOp->getConstTpetraOperator();
        if (tO->hasDiagonal()) {
          typedef Tpetra::Vector<double,int,panzer::GlobalOrdinal,NT> tV;
          RCP<tV> diag = rcp(new tV(tO->getRangeMap()));
          tO->getLocalDiagCopy(*diag);
          Tpetra::MatrixMarket::Writer<Tpetra::CrsMatrix<double,int,panzer::GlobalOrdinal,NT> >::writeDenseFile(("diag_"+s).c_str(),*diag);
        }
        *Teko::getOutputStream() << "Cannot dump operator \'" << s << "\'" << std::endl;
      }
    } else
      TEUCHOS_ASSERT(false);
  }


  void describeMatrix(const std::string & s,const Thyra::LinearOpBase<double> & op,Teuchos::RCP<Teuchos::FancyOStream> out)
  {
    using Teuchos::RCP;
    using Teuchos::rcp_dynamic_cast;
    using NT = Tpetra::Map<>::node_type;
    if (out!=Teuchos::null) {
      const RCP<const Thyra::TpetraLinearOp<double,int,panzer::GlobalOrdinal,NT> > tOp = rcp_dynamic_cast<const Thyra::TpetraLinearOp<double,int,panzer::GlobalOrdinal,NT> >(Teuchos::rcpFromRef(op));
      if(tOp != Teuchos::null) {
        const RCP<const Tpetra::CrsMatrix<double,int,panzer::GlobalOrdinal,NT> > crsOp = rcp_dynamic_cast<const Tpetra::CrsMatrix<double,int,panzer::GlobalOrdinal,NT> >(tOp->getConstTpetraOperator(),true);
        *out << "\nDebug: " << s << std::endl;
        crsOp->describe(*out,Teuchos::VERB_MEDIUM);
      } else
        TEUCHOS_ASSERT(false);
    }
  }

  void describeAndWriteMatrix(const std::string & s, const Thyra::LinearOpBase<double> & op, Teuchos::RCP<Teuchos::FancyOStream> out, const bool doWrite) {
    describeMatrix(s, op, out);
    if (doWrite)
      writeOut(s+".mm", op);
  }


  template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  Teuchos::RCP<const Tpetra::CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> > get_Tpetra_CrsMatrix(const Thyra::LinearOpBase<double> & op) {
    using Teuchos::RCP;
    using Teuchos::rcp_dynamic_cast;
    const RCP<const Thyra::TpetraLinearOp<Scalar,LocalOrdinal,GlobalOrdinal,Node> > tOp = rcp_dynamic_cast<const Thyra::TpetraLinearOp<Scalar,LocalOrdinal,GlobalOrdinal,Node> >(Teuchos::rcpFromRef(op),true);
    RCP<const Tpetra::CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> > crsOp = rcp_dynamic_cast<const Tpetra::CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> >(tOp->getConstTpetraOperator(),true);
    return crsOp;
  }

  template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  Teuchos::RCP<Tpetra::CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> >
  getIdentityMatrixTpetra (Teuchos::RCP<const Tpetra::Map<LocalOrdinal,GlobalOrdinal,Node> >& rowMap,
                           Scalar scaling)
  {
    typedef Tpetra::CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> Matrix;

    Teuchos::RCP<Matrix> identityMatrix = Teuchos::rcp(new Matrix(rowMap, rowMap, 1));
    Teuchos::Array<LocalOrdinal> col(1);
    Teuchos::Array<Scalar> val(1, scaling);
    for (LocalOrdinal rowLID = 0; rowLID < Teuchos::as<LocalOrdinal>(rowMap->getLocalNumElements()); rowLID++) {
      col[0] = rowLID;
      identityMatrix->insertLocalValues(rowLID, col(), val());
    }
    identityMatrix->fillComplete();
    return identityMatrix;
  }

  Teko::LinearOp getIdentityMatrix(const Teko::LinearOp& op, double scaling)
  {
    using Teuchos::RCP;
    using Teuchos::rcp_dynamic_cast;
    using Node = panzer::TpetraNodeType;

    const RCP<const Thyra::TpetraLinearOp<double,int,panzer::GlobalOrdinal,Node> > tOp = rcp_dynamic_cast<const Thyra::TpetraLinearOp<double,int,panzer::GlobalOrdinal,Node> >(op);
    if(tOp != Teuchos::null) {
      using Scalar = double;
      using LocalOrdinal = int;
      using GlobalOrdinal = panzer::GlobalOrdinal;
      RCP<Thyra::TpetraLinearOp<Scalar,LocalOrdinal,GlobalOrdinal,Node> > tOp2 = Teuchos::rcp_const_cast<Thyra::TpetraLinearOp<Scalar,LocalOrdinal,GlobalOrdinal,Node>>(tOp);
      RCP<Tpetra::CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> > crsOp = rcp_dynamic_cast<Tpetra::CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> >(tOp2->getTpetraOperator(),true);
      auto tpMap = crsOp->getRowMap();
      auto tpId = getIdentityMatrixTpetra(tpMap, scaling);
      Teko::LinearOp thyId = Thyra::tpetraLinearOp<Scalar,LocalOrdinal,GlobalOrdinal,Node>(Teko::rangeSpace(op),
                                                                                           Teko::domainSpace(op),
                                                                                           tpId);
      return thyId;
    } else {
      TEUCHOS_ASSERT(false);
      TEUCHOS_UNREACHABLE_RETURN(Teuchos::null);
    }
  }

  bool isMatrixFreeOperator(const Teko::LinearOp& op) {
    using Teuchos::RCP;
    using Teuchos::rcp_dynamic_cast;
    using Node = panzer::TpetraNodeType;

    const RCP<const Thyra::TpetraLinearOp<double,int,panzer::GlobalOrdinal,Node> > tOp = rcp_dynamic_cast<const Thyra::TpetraLinearOp<double,int,panzer::GlobalOrdinal,Node> >(op);
    if(tOp != Teuchos::null) {
      using Scalar = double;
      using LocalOrdinal = int;
      using GlobalOrdinal = panzer::GlobalOrdinal;
      RCP<Thyra::TpetraLinearOp<Scalar,LocalOrdinal,GlobalOrdinal,Node> > tOp2 = Teuchos::rcp_const_cast<Thyra::TpetraLinearOp<Scalar,LocalOrdinal,GlobalOrdinal,Node>>(tOp);
      RCP<Tpetra::CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> > crsOp = rcp_dynamic_cast<Tpetra::CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> >(tOp2->getTpetraOperator());
      return crsOp.is_null();
    } else
      return false;
  }

  Teko::LinearOp getLumpedInverseDiagonal(const Teko::LinearOp& op) {
    using Teuchos::RCP;
    using Scalar = double;
    RCP<Thyra::VectorBase<Scalar> > ones = Thyra::createMember(op->domain());
    RCP<Thyra::VectorBase<Scalar> > diagonal = Thyra::createMember(op->range());
    Thyra::assign(ones.ptr(),1.0);
    // compute lumped diagonal
    Thyra::apply(*op,Thyra::NOTRANS,*ones,diagonal.ptr());
    Thyra::reciprocal(*diagonal,diagonal.ptr());
    return Teuchos::rcp(new Thyra::DefaultDiagonalLinearOp<Scalar>(diagonal));
  }
}

#endif
