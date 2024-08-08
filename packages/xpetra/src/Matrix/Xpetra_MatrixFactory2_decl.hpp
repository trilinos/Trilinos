// @HEADER
// *****************************************************************************
//             Xpetra: A linear algebra interface package
//
// Copyright 2012 NTESS and the Xpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

// WARNING: This code is experimental. Backwards compatibility should not be expected.

#ifndef XPETRA_MATRIXFACTORY2_DECL_HPP
#define XPETRA_MATRIXFACTORY2_DECL_HPP

#include "Xpetra_ConfigDefs.hpp"
#include "Xpetra_MapExtractor_fwd.hpp"
#include "Xpetra_Matrix.hpp"
#include "Xpetra_CrsMatrixWrap.hpp"
#include "Xpetra_BlockedCrsMatrix_fwd.hpp"
#include "Xpetra_Map.hpp"
#include "Xpetra_BlockedMap.hpp"
#include "Xpetra_Vector.hpp"
#include "Xpetra_BlockedVector.hpp"
#include "Xpetra_Exceptions.hpp"

namespace Xpetra {

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node = Tpetra::KokkosClassic::DefaultNode::DefaultNodeType>
class MatrixFactory2 {
#undef XPETRA_MATRIXFACTORY2_SHORT
#include "Xpetra_UseShortNames.hpp"

 public:
  static RCP<Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> > BuildCopy(const RCP<const Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> > A, bool setFixedBlockSize = true);
};
#define XPETRA_MATRIXFACTORY2_SHORT

// template<>
// class MatrixFactory2<double,int,int,typename Xpetra::Matrix<double, int, int>::node_type> {
template <class Node>
class MatrixFactory2<double, int, int, Node> {
  typedef double Scalar;
  typedef int LocalOrdinal;
  typedef int GlobalOrdinal;
  // typedef Matrix<double, int, GlobalOrdinal>::node_type Node;
#undef XPETRA_MATRIXFACTORY2_SHORT
#include "Xpetra_UseShortNames.hpp"
 public:
  static RCP<Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> > BuildCopy(const RCP<const Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> > A, bool setFixedBlockSize = true) {
    RCP<const CrsMatrixWrap> oldOp = Teuchos::rcp_dynamic_cast<const CrsMatrixWrap>(A);
    if (oldOp == Teuchos::null)
      throw Exceptions::BadCast("Cast from Xpetra::Matrix to Xpetra::CrsMatrixWrap failed");

    RCP<const CrsMatrix> oldCrsOp = oldOp->getCrsMatrix();

#ifdef HAVE_XPETRA_EPETRA
#ifndef XPETRA_EPETRA_NO_32BIT_GLOBAL_INDICES
    RCP<const EpetraCrsMatrixT<GlobalOrdinal, Node> > oldECrsOp = Teuchos::rcp_dynamic_cast<const EpetraCrsMatrixT<GlobalOrdinal, Node> >(oldCrsOp);
    if (oldECrsOp != Teuchos::null) {
      // Underlying matrix is Epetra
      RCP<CrsMatrix> newECrsOp(new EpetraCrsMatrixT<GlobalOrdinal, Node>(*oldECrsOp));
      RCP<CrsMatrixWrap> newOp(new CrsMatrixWrap(newECrsOp));
      if (setFixedBlockSize)
        newOp->SetFixedBlockSize(A->GetFixedBlockSize());
      return newOp;
    }
#endif
#endif

#ifdef HAVE_XPETRA_TPETRA
    // Underlying matrix is Tpetra
    RCP<const TpetraCrsMatrix> oldTCrsOp = Teuchos::rcp_dynamic_cast<const TpetraCrsMatrix>(oldCrsOp);
    if (oldTCrsOp != Teuchos::null) {
      RCP<CrsMatrix> newTCrsOp(new TpetraCrsMatrix(*oldTCrsOp));
      RCP<CrsMatrixWrap> newOp(new CrsMatrixWrap(newTCrsOp));
      if (setFixedBlockSize)
        newOp->SetFixedBlockSize(A->GetFixedBlockSize());
      return newOp;
    }
    return Teuchos::null;
#else
    throw Exceptions::BadCast("Cast from Xpetra::Matrix to Xpetra::EpetraCrsMatrix or Xpetra::TpetraCrsMatrix failed");
    TEUCHOS_UNREACHABLE_RETURN(Teuchos::null);  // make compiler happy
#endif

  }  // BuildCopy
};

#define XPETRA_MATRIXFACTORY2_SHORT

#ifdef HAVE_XPETRA_INT_LONG_LONG
template <class Node>
class MatrixFactory2<double, int, long long, Node> {
  typedef double Scalar;
  typedef int LocalOrdinal;
  typedef long long GlobalOrdinal;
  // typedef Matrix<double, int, GlobalOrdinal>::node_type Node;
#undef XPETRA_MATRIXFACTORY2_SHORT
#include "Xpetra_UseShortNames.hpp"
 public:
  static RCP<Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> > BuildCopy(const RCP<const Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> > A, bool setFixedBlockSize = true) {
    RCP<const CrsMatrixWrap> oldOp = Teuchos::rcp_dynamic_cast<const CrsMatrixWrap>(A);
    if (oldOp == Teuchos::null)
      throw Exceptions::BadCast("Cast from Xpetra::Matrix to Xpetra::CrsMatrixWrap failed");

    RCP<const CrsMatrix> oldCrsOp = oldOp->getCrsMatrix();

#ifdef HAVE_XPETRA_EPETRA
#ifndef XPETRA_EPETRA_NO_64BIT_GLOBAL_INDICES
    RCP<const EpetraCrsMatrixT<GlobalOrdinal, Node> > oldECrsOp = Teuchos::rcp_dynamic_cast<const EpetraCrsMatrixT<GlobalOrdinal, Node> >(oldCrsOp);
    if (oldECrsOp != Teuchos::null) {
      // Underlying matrix is Epetra
      RCP<CrsMatrix> newECrsOp(new EpetraCrsMatrixT<GlobalOrdinal, Node>(*oldECrsOp));
      RCP<CrsMatrixWrap> newOp(new CrsMatrixWrap(newECrsOp));
      if (setFixedBlockSize)
        newOp->SetFixedBlockSize(A->GetFixedBlockSize());
      return newOp;
    }
#endif
#endif

#ifdef HAVE_XPETRA_TPETRA
    // Underlying matrix is Tpetra
    RCP<const TpetraCrsMatrix> oldTCrsOp = Teuchos::rcp_dynamic_cast<const TpetraCrsMatrix>(oldCrsOp);
    if (oldTCrsOp != Teuchos::null) {
      RCP<CrsMatrix> newTCrsOp(new TpetraCrsMatrix(*oldTCrsOp));
      RCP<CrsMatrixWrap> newOp(new CrsMatrixWrap(newTCrsOp));
      if (setFixedBlockSize)
        newOp->SetFixedBlockSize(A->GetFixedBlockSize());
      return newOp;
    }
#else
    throw Exceptions::BadCast("Cast from Xpetra::Matrix to Xpetra::EpetraCrsMatrix or Xpetra::TpetraCrsMatrix failed");
#endif

    return Teuchos::null;  // make compiler happy
  }
};
#endif  // HAVE_XPETRA_INT_LONG_LONG

#define XPETRA_MATRIXFACTORY2_SHORT

}  // namespace Xpetra

#define XPETRA_MATRIXFACTORY2_SHORT

#endif  // ifndef XPETRA_MATRIXFACTORY2_DECL_HPP
