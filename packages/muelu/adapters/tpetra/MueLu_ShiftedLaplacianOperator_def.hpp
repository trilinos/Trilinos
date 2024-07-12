// @HEADER
// *****************************************************************************
//        MueLu: A package for multigrid based preconditioning
//
// Copyright 2012 NTESS and the MueLu contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef MUELU_SHIFTEDLAPLACIANOPERATOR_DEF_HPP
#define MUELU_SHIFTEDLAPLACIANOPERATOR_DEF_HPP

#include "MueLu_ConfigDefs.hpp"

#include <Xpetra_Matrix.hpp>
#include <Xpetra_CrsMatrixWrap.hpp>
#include <Xpetra_BlockedCrsMatrix.hpp>
#include <Xpetra_TpetraMultiVector.hpp>
#include <Xpetra_MultiVectorFactory.hpp>

#include "MueLu_ShiftedLaplacianOperator_decl.hpp"
#include "MueLu_Hierarchy.hpp"
#include "MueLu_Utilities.hpp"

namespace MueLu {

// ------------- getDomainMap -----------------------

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
Teuchos::RCP<const Tpetra::Map<LocalOrdinal, GlobalOrdinal, Node> >
ShiftedLaplacianOperator<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    getDomainMap() const {
  typedef Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> XMatrix;

  RCP<MueLu::Level> L0 = Hierarchy_->GetLevel(0);
  RCP<XMatrix> A       = L0->Get<RCP<XMatrix> >("A");

  RCP<Xpetra::BlockedCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> > tpbA =
      Teuchos::rcp_dynamic_cast<Xpetra::BlockedCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> >(A);
  if (tpbA != Teuchos::null) {
    return Xpetra::toTpetraNonZero(tpbA->getDomainMap());
  }

  RCP<Tpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> > tpA =
      Utilities<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Op2NonConstTpetraCrs(A);
  return tpA->getDomainMap();
}

// ------------- getRangeMap -----------------------

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
Teuchos::RCP<const Tpetra::Map<LocalOrdinal, GlobalOrdinal, Node> >
ShiftedLaplacianOperator<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    getRangeMap() const {
  typedef Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> XMatrix;

  RCP<MueLu::Level> L0 = Hierarchy_->GetLevel(0);
  RCP<XMatrix> A       = L0->Get<RCP<XMatrix> >("A");

  RCP<Xpetra::BlockedCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> > tpbA =
      Teuchos::rcp_dynamic_cast<Xpetra::BlockedCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> >(A);
  if (tpbA != Teuchos::null)
    return Xpetra::toTpetraNonZero(tpbA->getRangeMap());

  RCP<Tpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> > tpA =
      Utilities<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Op2NonConstTpetraCrs(A);
  return tpA->getRangeMap();
}

// ------------- apply -----------------------

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void ShiftedLaplacianOperator<Scalar, LocalOrdinal, GlobalOrdinal, Node>::apply(const Tpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>& X,
                                                                                Tpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>& Y,
                                                                                Teuchos::ETransp /* mode */, Scalar /* alpha */, Scalar /* beta */) const {
  typedef Tpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node> TMV;
  typedef Xpetra::TpetraMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node> XTMV;
  // typedef Xpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>        XMV; // unused

  TMV& temp_x = const_cast<TMV&>(X);
  const XTMV tX(rcpFromRef(temp_x));
  XTMV tY(rcpFromRef(Y));

  try {
    tY.putScalar(0.0);
    Hierarchy_->Iterate(tX, tY, cycles_, true);
  }

  catch (std::exception& e) {
    // FIXME add message and rethrow
    std::cerr << "Caught an exception in MueLu::ShiftedLaplacianOperator::ApplyInverse():" << std::endl
              << e.what() << std::endl;
  }

  // update solution with 2-grid error correction
  /*if(option_==1) {
    for(int j=0; j<cycles_; j++) {
      RCP<XMV> residual       = MueLu::Utilities<Scalar,LocalOrdinal,GlobalOrdinal,Node>::Residual(*A_, tY, tX);
      RCP<XMV> coarseResidual = Xpetra::MultiVectorFactory<Scalar,LocalOrdinal,GlobalOrdinal,Node>::Build(R_->getRangeMap(), tX.getNumVectors());
      RCP<XMV> coarseError    = Xpetra::MultiVectorFactory<Scalar,LocalOrdinal,GlobalOrdinal,Node>::Build(R_->getRangeMap(), tX.getNumVectors());
      R_ -> apply(*residual, *coarseResidual, Teuchos::NO_TRANS, (Scalar) 1.0, (Scalar) 0.0);
      RCP<TMV> tcoarseR = MueLu::Utilities<Scalar,LocalOrdinal,GlobalOrdinal,Node>::MV2NonConstTpetraMV(coarseResidual);
      RCP<TMV> tcoarseE = MueLu::Utilities<Scalar,LocalOrdinal,GlobalOrdinal,Node>::MV2NonConstTpetraMV(coarseError);
      BelosLP_ -> setProblem(tcoarseE,tcoarseR);
      BelosSM_ -> solve();
      RCP<XMV> fineError      = Xpetra::MultiVectorFactory<Scalar,LocalOrdinal,GlobalOrdinal,Node>::Build(P_->getRangeMap(), tX.getNumVectors());
      XTMV tmpcoarseE(rcpFromRef(*tcoarseE));
      P_ -> apply(tmpcoarseE, *fineError, Teuchos::NO_TRANS, (Scalar) 1.0, (Scalar) 0.0);
      tY.update((Scalar) 1.0, *fineError, (Scalar) 1.0);
    }
  }

  try {
    Hierarchy_->Iterate(tX, tY, 1, false);
  }

  catch(std::exception& e) {
    //FIXME add message and rethrow
    std::cerr << "Caught an exception in MueLu::ShiftedLaplacianOperator::ApplyInverse():" << std::endl
    << e.what() << std::endl;
    }*/
}

// ------------- apply -----------------------
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
bool ShiftedLaplacianOperator<Scalar, LocalOrdinal, GlobalOrdinal, Node>::hasTransposeApply() const {
  return false;
}

}  // namespace MueLu

#endif  // ifdef MUELU_SHIFTEDLAPLACIANOPERATOR_DEF_HPP
