// @HEADER
// *****************************************************************************
//        MueLu: A package for multigrid based preconditioning
//
// Copyright 2012 NTESS and the MueLu contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef MUELU_UTILITIES_DEF_HPP
#define MUELU_UTILITIES_DEF_HPP

#include <Teuchos_DefaultComm.hpp>
#include <Teuchos_ParameterList.hpp>
#include <iostream>

#include "MueLu_ConfigDefs.hpp"
#include "Xpetra_TpetraRowMatrix.hpp"

#ifdef HAVE_MUELU_EPETRA
#ifdef HAVE_MPI
#include "Epetra_MpiComm.h"
#endif
#endif

#if defined(HAVE_MUELU_EPETRA) && defined(HAVE_MUELU_EPETRAEXT)
#include <EpetraExt_MatrixMatrix.h>
#include <EpetraExt_RowMatrixOut.h>
#include <EpetraExt_MultiVectorOut.h>
#include <EpetraExt_CrsMatrixIn.h>
#include <EpetraExt_MultiVectorIn.h>
#include <EpetraExt_BlockMapIn.h>
#include <Xpetra_EpetraUtils.hpp>
#include <Xpetra_EpetraMultiVector.hpp>
#include <EpetraExt_BlockMapOut.h>
#endif

#include <MatrixMarket_Tpetra.hpp>
#include <Tpetra_RowMatrixTransposer.hpp>
#include <TpetraExt_MatrixMatrix.hpp>
#include <Xpetra_TpetraMultiVector.hpp>
#include <Xpetra_TpetraOperator.hpp>
#include <Xpetra_TpetraCrsMatrix.hpp>
#include <Xpetra_TpetraBlockCrsMatrix.hpp>

#ifdef HAVE_MUELU_EPETRA
#include <Xpetra_EpetraMap.hpp>
#endif

#include <Xpetra_BlockedCrsMatrix.hpp>
//#include <Xpetra_DefaultPlatform.hpp>
#include <Xpetra_IO.hpp>
#include <Xpetra_Map.hpp>
#include <Xpetra_MapFactory.hpp>
#include <Xpetra_Matrix.hpp>
#include <Xpetra_MatrixFactory.hpp>
#include <Xpetra_MultiVector.hpp>
#include <Xpetra_MultiVectorFactory.hpp>
#include <Xpetra_Operator.hpp>
#include <Xpetra_Vector.hpp>
#include <Xpetra_VectorFactory.hpp>

#include <Xpetra_MatrixMatrix.hpp>

#include <MueLu_Utilities_decl.hpp>
#if defined(HAVE_MUELU_EPETRA) && defined(HAVE_MUELU_ML)
#include <ml_operator.h>
#include <ml_epetra_utils.h>
#endif

namespace MueLu {

#ifdef HAVE_MUELU_EPETRA
// using Xpetra::EpetraCrsMatrix;   // TODO: mv in Xpetra_UseShortNamesScalar
// using Xpetra::EpetraMultiVector;
#endif

#ifdef HAVE_MUELU_EPETRA
template <typename SC, typename LO, typename GO, typename NO>
RCP<Xpetra::CrsMatrixWrap<SC, LO, GO, NO>> Convert_Epetra_CrsMatrix_ToXpetra_CrsMatrixWrap(RCP<Epetra_CrsMatrix>& epAB) {
  return Xpetra::Convert_Epetra_CrsMatrix_ToXpetra_CrsMatrixWrap<SC, LO, GO, NO>(epAB);
}
#endif

#ifdef HAVE_MUELU_EPETRA
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
RCP<const Epetra_MultiVector> Utilities<Scalar, LocalOrdinal, GlobalOrdinal, Node>::MV2EpetraMV(const RCP<Xpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>> vec) {
  RCP<const Xpetra::EpetraMultiVectorT<GlobalOrdinal, Node>> tmpVec = rcp_dynamic_cast<Xpetra::EpetraMultiVectorT<GlobalOrdinal, Node>>(vec);
  if (tmpVec == Teuchos::null)
    throw Exceptions::BadCast("Cast from Xpetra::MultiVector to Xpetra::EpetraMultiVector failed");
  return tmpVec->getEpetra_MultiVector();
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
RCP<Epetra_MultiVector> Utilities<Scalar, LocalOrdinal, GlobalOrdinal, Node>::MV2NonConstEpetraMV(RCP<Xpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>> vec) {
  RCP<const Xpetra::EpetraMultiVectorT<GlobalOrdinal, Node>> tmpVec = rcp_dynamic_cast<Xpetra::EpetraMultiVectorT<GlobalOrdinal, Node>>(vec);
  if (tmpVec == Teuchos::null)
    throw Exceptions::BadCast("Cast from Xpetra::MultiVector to Xpetra::EpetraMultiVector failed");
  return tmpVec->getEpetra_MultiVector();
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
Epetra_MultiVector& Utilities<Scalar, LocalOrdinal, GlobalOrdinal, Node>::MV2NonConstEpetraMV(Xpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>& vec) {
  const Xpetra::EpetraMultiVectorT<GlobalOrdinal, Node>& tmpVec = dynamic_cast<const Xpetra::EpetraMultiVectorT<GlobalOrdinal, Node>&>(vec);
  return *(tmpVec.getEpetra_MultiVector());
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
const Epetra_MultiVector& Utilities<Scalar, LocalOrdinal, GlobalOrdinal, Node>::MV2EpetraMV(const Xpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>& vec) {
  const Xpetra::EpetraMultiVectorT<GlobalOrdinal, Node>& tmpVec = dynamic_cast<const Xpetra::EpetraMultiVectorT<GlobalOrdinal, Node>&>(vec);
  return *(tmpVec.getEpetra_MultiVector());
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
RCP<const Epetra_CrsMatrix> Utilities<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Op2EpetraCrs(RCP<const Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>> Op) {
  RCP<const Xpetra::CrsMatrixWrap<Scalar, LocalOrdinal, GlobalOrdinal, Node>> crsOp = rcp_dynamic_cast<const Xpetra::CrsMatrixWrap<Scalar, LocalOrdinal, GlobalOrdinal, Node>>(Op);
  if (crsOp == Teuchos::null)
    throw Exceptions::BadCast("Cast from Xpetra::Matrix to Xpetra::CrsMatrixWrap failed");
  const RCP<const Xpetra::EpetraCrsMatrixT<GlobalOrdinal, Node>>& tmp_ECrsMtx = rcp_dynamic_cast<const Xpetra::EpetraCrsMatrixT<GlobalOrdinal, Node>>(crsOp->getCrsMatrix());
  if (tmp_ECrsMtx == Teuchos::null)
    throw Exceptions::BadCast("Cast from Xpetra::CrsMatrix to Xpetra::EpetraCrsMatrix failed");
  return tmp_ECrsMtx->getEpetra_CrsMatrix();
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
RCP<Epetra_CrsMatrix> Utilities<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Op2NonConstEpetraCrs(RCP<Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>> Op) {
  RCP<const Xpetra::CrsMatrixWrap<Scalar, LocalOrdinal, GlobalOrdinal, Node>> crsOp = rcp_dynamic_cast<const Xpetra::CrsMatrixWrap<Scalar, LocalOrdinal, GlobalOrdinal, Node>>(Op);
  if (crsOp == Teuchos::null)
    throw Exceptions::BadCast("Cast from Xpetra::Matrix to Xpetra::CrsMatrixWrap failed");
  const RCP<const Xpetra::EpetraCrsMatrixT<GlobalOrdinal, Node>>& tmp_ECrsMtx = rcp_dynamic_cast<const Xpetra::EpetraCrsMatrixT<GlobalOrdinal, Node>>(crsOp->getCrsMatrix());
  if (tmp_ECrsMtx == Teuchos::null)
    throw Exceptions::BadCast("Cast from Xpetra::CrsMatrix to Xpetra::EpetraCrsMatrix failed");
  return tmp_ECrsMtx->getEpetra_CrsMatrixNonConst();
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
const Epetra_CrsMatrix& Utilities<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Op2EpetraCrs(const Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>& Op) {
  try {
    const Xpetra::CrsMatrixWrap<Scalar, LocalOrdinal, GlobalOrdinal, Node>& crsOp = dynamic_cast<const Xpetra::CrsMatrixWrap<Scalar, LocalOrdinal, GlobalOrdinal, Node>&>(Op);
    try {
      const Xpetra::EpetraCrsMatrixT<GlobalOrdinal, Node>& tmp_ECrsMtx = dynamic_cast<const Xpetra::EpetraCrsMatrixT<GlobalOrdinal, Node>&>(*crsOp.getCrsMatrix());
      return *tmp_ECrsMtx.getEpetra_CrsMatrix();
    } catch (std::bad_cast&) {
      throw Exceptions::BadCast("Cast from Xpetra::CrsMatrix to Xpetra::EpetraCrsMatrix failed");
    }
  } catch (std::bad_cast&) {
    throw Exceptions::BadCast("Cast from Xpetra::Matrix to Xpetra::CrsMatrixWrap failed");
  }
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
Epetra_CrsMatrix& Utilities<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Op2NonConstEpetraCrs(Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>& Op) {
  try {
    Xpetra::CrsMatrixWrap<Scalar, LocalOrdinal, GlobalOrdinal, Node>& crsOp = dynamic_cast<Xpetra::CrsMatrixWrap<Scalar, LocalOrdinal, GlobalOrdinal, Node>&>(Op);
    try {
      Xpetra::EpetraCrsMatrixT<GlobalOrdinal, Node>& tmp_ECrsMtx = dynamic_cast<Xpetra::EpetraCrsMatrixT<GlobalOrdinal, Node>&>(*crsOp.getCrsMatrix());
      return *tmp_ECrsMtx.getEpetra_CrsMatrixNonConst();
    } catch (std::bad_cast&) {
      throw Exceptions::BadCast("Cast from Xpetra::CrsMatrix to Xpetra::EpetraCrsMatrix failed");
    }
  } catch (std::bad_cast&) {
    throw Exceptions::BadCast("Cast from Xpetra::Matrix to Xpetra::CrsMatrixWrap failed");
  }
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
const Epetra_Map& Utilities<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Map2EpetraMap(const Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node>& map) {
  RCP<const Xpetra::EpetraMapT<GlobalOrdinal, Node>> xeMap = rcp_dynamic_cast<const Xpetra::EpetraMapT<GlobalOrdinal, Node>>(rcpFromRef(map));
  if (xeMap == Teuchos::null)
    throw Exceptions::BadCast("Utilities::Map2EpetraMap : Cast from Xpetra::Map to Xpetra::EpetraMap failed");
  return xeMap->getEpetra_Map();
}
#endif

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
RCP<const Tpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>>
Utilities<Scalar, LocalOrdinal, GlobalOrdinal, Node>::MV2TpetraMV(RCP<Xpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>> const vec) {
  RCP<const Xpetra::TpetraMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>> tmpVec = rcp_dynamic_cast<Xpetra::TpetraMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>>(vec);
  if (tmpVec == Teuchos::null)
    throw Exceptions::BadCast("Cast from Xpetra::MultiVector to Xpetra::TpetraMultiVector failed");
  return tmpVec->getTpetra_MultiVector();
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
RCP<Tpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>> Utilities<Scalar, LocalOrdinal, GlobalOrdinal, Node>::MV2NonConstTpetraMV(RCP<Xpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>> vec) {
  RCP<const Xpetra::TpetraMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>> tmpVec = rcp_dynamic_cast<Xpetra::TpetraMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>>(vec);
  if (tmpVec == Teuchos::null)
    throw Exceptions::BadCast("Cast from Xpetra::MultiVector to Xpetra::TpetraMultiVector failed");
  return tmpVec->getTpetra_MultiVector();
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
Tpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>& Utilities<Scalar, LocalOrdinal, GlobalOrdinal, Node>::MV2NonConstTpetraMV(Xpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>& vec) {
  const Xpetra::TpetraMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>& tmpVec = dynamic_cast<const Xpetra::TpetraMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>&>(vec);
  return *(tmpVec.getTpetra_MultiVector());
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
RCP<Tpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>> Utilities<Scalar, LocalOrdinal, GlobalOrdinal, Node>::MV2NonConstTpetraMV2(Xpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>& vec) {
  const Xpetra::TpetraMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>& tmpVec = dynamic_cast<const Xpetra::TpetraMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>&>(vec);
  return tmpVec.getTpetra_MultiVector();
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
const Tpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>&
Utilities<Scalar, LocalOrdinal, GlobalOrdinal, Node>::MV2TpetraMV(const Xpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>& vec) {
  const Xpetra::TpetraMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>& tmpVec = dynamic_cast<const Xpetra::TpetraMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>&>(vec);
  return *(tmpVec.getTpetra_MultiVector());
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
RCP<const Tpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>> Utilities<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Op2TpetraCrs(RCP<const Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>> Op) {
  // Get the underlying Tpetra Mtx
  RCP<const Xpetra::CrsMatrixWrap<Scalar, LocalOrdinal, GlobalOrdinal, Node>> crsOp = rcp_dynamic_cast<const Xpetra::CrsMatrixWrap<Scalar, LocalOrdinal, GlobalOrdinal, Node>>(Op);
  if (crsOp == Teuchos::null)
    throw Exceptions::BadCast("Cast from Xpetra::Matrix to Xpetra::CrsMatrixWrap failed");
  const RCP<const Xpetra::TpetraCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>>& tmp_ECrsMtx = rcp_dynamic_cast<const Xpetra::TpetraCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>>(crsOp->getCrsMatrix());
  if (tmp_ECrsMtx == Teuchos::null)
    throw Exceptions::BadCast("Cast from Xpetra::CrsMatrix to Xpetra::TpetraCrsMatrix failed");
  return tmp_ECrsMtx->getTpetra_CrsMatrix();
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
RCP<Tpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>> Utilities<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Op2NonConstTpetraCrs(RCP<Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>> Op) {
  RCP<const Xpetra::CrsMatrixWrap<Scalar, LocalOrdinal, GlobalOrdinal, Node>> crsOp = rcp_dynamic_cast<const Xpetra::CrsMatrixWrap<Scalar, LocalOrdinal, GlobalOrdinal, Node>>(Op);
  if (crsOp == Teuchos::null)
    throw Exceptions::BadCast("Cast from Xpetra::Matrix to Xpetra::CrsMatrixWrap failed");
  const RCP<const Xpetra::TpetraCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>>& tmp_ECrsMtx = rcp_dynamic_cast<const Xpetra::TpetraCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>>(crsOp->getCrsMatrix());
  if (tmp_ECrsMtx == Teuchos::null)
    throw Exceptions::BadCast("Cast from Xpetra::CrsMatrix to Xpetra::TpetraCrsMatrix failed");
  return tmp_ECrsMtx->getTpetra_CrsMatrixNonConst();
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
const Tpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>& Utilities<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Op2TpetraCrs(const Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>& Op) {
  try {
    const Xpetra::CrsMatrixWrap<Scalar, LocalOrdinal, GlobalOrdinal, Node>& crsOp = dynamic_cast<const Xpetra::CrsMatrixWrap<Scalar, LocalOrdinal, GlobalOrdinal, Node>&>(Op);
    try {
      const Xpetra::TpetraCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>& tmp_ECrsMtx = dynamic_cast<const Xpetra::TpetraCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>&>(*crsOp.getCrsMatrix());
      return *tmp_ECrsMtx.getTpetra_CrsMatrix();
    } catch (std::bad_cast&) {
      throw Exceptions::BadCast("Cast from Xpetra::CrsMatrix to Xpetra::TpetraCrsMatrix failed");
    }
  } catch (std::bad_cast&) {
    throw Exceptions::BadCast("Cast from Xpetra::Matrix to Xpetra::CrsMatrixWrap failed");
  }
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
Tpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>& Utilities<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Op2NonConstTpetraCrs(Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>& Op) {
  try {
    Xpetra::CrsMatrixWrap<Scalar, LocalOrdinal, GlobalOrdinal, Node>& crsOp = dynamic_cast<Xpetra::CrsMatrixWrap<Scalar, LocalOrdinal, GlobalOrdinal, Node>&>(Op);
    try {
      Xpetra::TpetraCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>& tmp_ECrsMtx = dynamic_cast<Xpetra::TpetraCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>&>(*crsOp.getCrsMatrix());
      return *tmp_ECrsMtx.getTpetra_CrsMatrixNonConst();
    } catch (std::bad_cast&) {
      throw Exceptions::BadCast("Cast from Xpetra::CrsMatrix to Xpetra::TpetraCrsMatrix failed");
    }
  } catch (std::bad_cast&) {
    throw Exceptions::BadCast("Cast from Xpetra::Matrix to Xpetra::CrsMatrixWrap failed");
  }
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
RCP<const Tpetra::BlockCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>> Utilities<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Op2TpetraBlockCrs(RCP<const Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>> Op) {
  using XCrsMatrixWrap = Xpetra::CrsMatrixWrap<Scalar, LocalOrdinal, GlobalOrdinal, Node>;
  // Get the underlying Tpetra Mtx
  RCP<const XCrsMatrixWrap> crsOp = rcp_dynamic_cast<const XCrsMatrixWrap>(Op);
  if (crsOp == Teuchos::null)
    throw Exceptions::BadCast("Cast from Xpetra::Matrix to Xpetra::CrsMatrixWrap failed");
  const RCP<const Xpetra::TpetraBlockCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>>& tmp_ECrsMtx = rcp_dynamic_cast<const Xpetra::TpetraBlockCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>>(crsOp->getCrsMatrix());
  if (tmp_ECrsMtx == Teuchos::null)
    throw Exceptions::BadCast("Cast from Xpetra::CrsMatrix to Xpetra::TpetraBlockCrsMatrix failed");
  return tmp_ECrsMtx->getTpetra_BlockCrsMatrix();
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
RCP<Tpetra::BlockCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>> Utilities<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Op2NonConstTpetraBlockCrs(RCP<Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>> Op) {
  using XCrsMatrixWrap            = Xpetra::CrsMatrixWrap<Scalar, LocalOrdinal, GlobalOrdinal, Node>;
  RCP<const XCrsMatrixWrap> crsOp = rcp_dynamic_cast<const XCrsMatrixWrap>(Op);
  if (crsOp == Teuchos::null)
    throw Exceptions::BadCast("Cast from Xpetra::Matrix to Xpetra::CrsMatrixWrap failed");
  const RCP<const Xpetra::TpetraBlockCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>>& tmp_ECrsMtx = rcp_dynamic_cast<const Xpetra::TpetraBlockCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>>(crsOp->getCrsMatrix());
  if (tmp_ECrsMtx == Teuchos::null)
    throw Exceptions::BadCast("Cast from Xpetra::CrsMatrix to Xpetra::TpetraBlockCrsMatrix failed");
  return tmp_ECrsMtx->getTpetra_BlockCrsMatrixNonConst();
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
const Tpetra::BlockCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>& Utilities<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Op2TpetraBlockCrs(const Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>& Op) {
  try {
    using XCrsMatrixWrap        = Xpetra::CrsMatrixWrap<Scalar, LocalOrdinal, GlobalOrdinal, Node>;
    const XCrsMatrixWrap& crsOp = dynamic_cast<const XCrsMatrixWrap&>(Op);
    try {
      const Xpetra::TpetraBlockCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>& tmp_ECrsMtx = dynamic_cast<const Xpetra::TpetraBlockCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>&>(*crsOp.getCrsMatrix());
      return *tmp_ECrsMtx.getTpetra_BlockCrsMatrix();
    } catch (std::bad_cast&) {
      throw Exceptions::BadCast("Cast from Xpetra::CrsMatrix to Xpetra::TpetraBlockCrsMatrix failed");
    }
  } catch (std::bad_cast&) {
    throw Exceptions::BadCast("Cast from Xpetra::Matrix to Xpetra::CrsMatrixWrap failed");
  }
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
Tpetra::BlockCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>& Utilities<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Op2NonConstTpetraBlockCrs(Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>& Op) {
  try {
    using XCrsMatrixWrap  = Xpetra::CrsMatrixWrap<Scalar, LocalOrdinal, GlobalOrdinal, Node>;
    XCrsMatrixWrap& crsOp = dynamic_cast<XCrsMatrixWrap&>(Op);
    try {
      Xpetra::TpetraBlockCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>& tmp_ECrsMtx = dynamic_cast<Xpetra::TpetraBlockCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>&>(*crsOp.getCrsMatrix());
      return *tmp_ECrsMtx.getTpetra_BlockCrsMatrixNonConst();
    } catch (std::bad_cast&) {
      throw Exceptions::BadCast("Cast from Xpetra::CrsMatrix to Xpetra::TpetraBlockCrsMatrix failed");
    }
  } catch (std::bad_cast&) {
    throw Exceptions::BadCast("Cast from Xpetra::Matrix to Xpetra::CrsMatrixWrap failed");
  }
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
RCP<const Tpetra::RowMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>> Utilities<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Op2TpetraRow(RCP<const Xpetra::Operator<Scalar, LocalOrdinal, GlobalOrdinal, Node>> Op) {
  RCP<const Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>> mat           = rcp_dynamic_cast<const Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>>(Op);
  RCP<const Xpetra::TpetraRowMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>> rmat = rcp_dynamic_cast<const Xpetra::TpetraRowMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>>(Op);
  if (!mat.is_null()) {
    RCP<const Xpetra::CrsMatrixWrap<Scalar, LocalOrdinal, GlobalOrdinal, Node>> crsOp = rcp_dynamic_cast<const Xpetra::CrsMatrixWrap<Scalar, LocalOrdinal, GlobalOrdinal, Node>>(mat);
    if (crsOp == Teuchos::null)
      throw Exceptions::BadCast("Cast from Xpetra::Matrix to Xpetra::CrsMatrixWrap failed");

    RCP<const Xpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>> crsMat              = crsOp->getCrsMatrix();
    const RCP<const Xpetra::TpetraCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>> tmp_Crs = rcp_dynamic_cast<const Xpetra::TpetraCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>>(crsMat);
    RCP<const Xpetra::TpetraBlockCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>> tmp_BlockCrs;
    if (!tmp_Crs.is_null()) {
      return tmp_Crs->getTpetra_CrsMatrixNonConst();
    } else {
      tmp_BlockCrs = rcp_dynamic_cast<const Xpetra::TpetraBlockCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>>(crsMat);
      if (tmp_BlockCrs.is_null())
        throw Exceptions::BadCast("Cast from Xpetra::CrsMatrix to Xpetra::TpetraCrsMatrix and Xpetra::TpetraBlockCrsMatrix failed");
      return tmp_BlockCrs->getTpetra_BlockCrsMatrixNonConst();
    }
  } else if (!rmat.is_null()) {
    return rmat->getTpetra_RowMatrix();
  } else {
    RCP<const Xpetra::TpetraOperator<Scalar, LocalOrdinal, GlobalOrdinal, Node>> tpOp = rcp_dynamic_cast<const Xpetra::TpetraOperator<Scalar, LocalOrdinal, GlobalOrdinal, Node>>(Op, true);
    RCP<const Tpetra::Operator<Scalar, LocalOrdinal, GlobalOrdinal, Node>> tOp        = tpOp->getOperatorConst();
    RCP<const Tpetra::RowMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>> tRow      = rcp_dynamic_cast<const Tpetra::RowMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>>(tOp, true);
    return tRow;
  }
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
RCP<Tpetra::RowMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>> Utilities<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Op2NonConstTpetraRow(RCP<Xpetra::Operator<Scalar, LocalOrdinal, GlobalOrdinal, Node>> Op) {
  RCP<Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>> mat           = rcp_dynamic_cast<Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>>(Op);
  RCP<Xpetra::TpetraRowMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>> rmat = rcp_dynamic_cast<Xpetra::TpetraRowMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>>(Op);
  if (!mat.is_null()) {
    RCP<const Xpetra::CrsMatrixWrap<Scalar, LocalOrdinal, GlobalOrdinal, Node>> crsOp = rcp_dynamic_cast<const Xpetra::CrsMatrixWrap<Scalar, LocalOrdinal, GlobalOrdinal, Node>>(mat);
    if (crsOp == Teuchos::null)
      throw Exceptions::BadCast("Cast from Xpetra::Matrix to Xpetra::CrsMatrixWrap failed");

    RCP<const Xpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>> crsMat              = crsOp->getCrsMatrix();
    const RCP<const Xpetra::TpetraCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>> tmp_Crs = rcp_dynamic_cast<const Xpetra::TpetraCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>>(crsMat);
    RCP<const Xpetra::TpetraBlockCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>> tmp_BlockCrs;
    if (!tmp_Crs.is_null()) {
      return tmp_Crs->getTpetra_CrsMatrixNonConst();
    } else {
      tmp_BlockCrs = rcp_dynamic_cast<const Xpetra::TpetraBlockCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>>(crsMat);
      if (tmp_BlockCrs.is_null())
        throw Exceptions::BadCast("Cast from Xpetra::CrsMatrix to Xpetra::TpetraCrsMatrix and Xpetra::TpetraBlockCrsMatrix failed");
      return tmp_BlockCrs->getTpetra_BlockCrsMatrixNonConst();
    }
  } else if (!rmat.is_null()) {
    return rmat->getTpetra_RowMatrixNonConst();
  } else {
    RCP<Xpetra::TpetraOperator<Scalar, LocalOrdinal, GlobalOrdinal, Node>> tpOp = rcp_dynamic_cast<Xpetra::TpetraOperator<Scalar, LocalOrdinal, GlobalOrdinal, Node>>(Op, true);
    RCP<Tpetra::Operator<Scalar, LocalOrdinal, GlobalOrdinal, Node>> tOp        = tpOp->getOperator();
    RCP<Tpetra::RowMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>> tRow      = rcp_dynamic_cast<Tpetra::RowMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>>(tOp, true);
    return tRow;
  }
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
const RCP<const Tpetra::Map<LocalOrdinal, GlobalOrdinal, Node>> Utilities<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Map2TpetraMap(const Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node>& map) {
  const RCP<const Xpetra::TpetraMap<LocalOrdinal, GlobalOrdinal, Node>>& tmp_TMap = rcp_dynamic_cast<const Xpetra::TpetraMap<LocalOrdinal, GlobalOrdinal, Node>>(rcpFromRef(map));
  if (tmp_TMap == Teuchos::null)
    throw Exceptions::BadCast("Utilities::Map2TpetraMap : Cast from Xpetra::Map to Xpetra::TpetraMap failed");
  return tmp_TMap->getTpetra_Map();
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void Utilities<Scalar, LocalOrdinal, GlobalOrdinal, Node>::MyOldScaleMatrix(Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>& Op, const Teuchos::ArrayRCP<const Scalar>& scalingVector, bool doInverse,
                                                                            bool doFillComplete,
                                                                            bool doOptimizeStorage) {
  Scalar one = Teuchos::ScalarTraits<Scalar>::one();
  Teuchos::ArrayRCP<Scalar> sv(scalingVector.size());
  if (doInverse) {
    for (int i = 0; i < scalingVector.size(); ++i)
      sv[i] = one / scalingVector[i];
  } else {
    for (int i = 0; i < scalingVector.size(); ++i)
      sv[i] = scalingVector[i];
  }

  switch (Op.getRowMap()->lib()) {
    case Xpetra::UseTpetra:
      MyOldScaleMatrix_Tpetra(Op, sv, doFillComplete, doOptimizeStorage);
      break;

    case Xpetra::UseEpetra:
      MyOldScaleMatrix_Epetra(Op, sv, doFillComplete, doOptimizeStorage);
      break;

    default:
      throw Exceptions::RuntimeError("Only Epetra and Tpetra matrices can be scaled.");
  }
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void Utilities<Scalar, LocalOrdinal, GlobalOrdinal, Node>::MyOldScaleMatrix_Epetra(Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>& /* Op */, const Teuchos::ArrayRCP<Scalar>& /* scalingVector */, bool /* doFillComplete */, bool /* doOptimizeStorage */) {
  throw Exceptions::RuntimeError("MyOldScaleMatrix_Epetra: Epetra needs SC=double and LO=GO=int.");
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void Utilities<Scalar, LocalOrdinal, GlobalOrdinal, Node>::MyOldScaleMatrix_Tpetra(Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>& Op, const Teuchos::ArrayRCP<Scalar>& scalingVector,
                                                                                   bool doFillComplete,
                                                                                   bool doOptimizeStorage) {
  try {
    Tpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>& tpOp = Op2NonConstTpetraCrs(Op);

    const RCP<const Tpetra::Map<LocalOrdinal, GlobalOrdinal, Node>> rowMap    = tpOp.getRowMap();
    const RCP<const Tpetra::Map<LocalOrdinal, GlobalOrdinal, Node>> domainMap = tpOp.getDomainMap();
    const RCP<const Tpetra::Map<LocalOrdinal, GlobalOrdinal, Node>> rangeMap  = tpOp.getRangeMap();

    size_t maxRowSize = tpOp.getLocalMaxNumRowEntries();
    if (maxRowSize == Teuchos::as<size_t>(-1))  // hasn't been determined yet
      maxRowSize = 20;

    if (tpOp.isFillComplete())
      tpOp.resumeFill();

    if (Op.isLocallyIndexed() == true) {
      typename Tpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::local_inds_host_view_type cols;
      typename Tpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::values_host_view_type vals;

      for (size_t i = 0; i < rowMap->getLocalNumElements(); ++i) {
        tpOp.getLocalRowView(i, cols, vals);

        size_t nnz = tpOp.getNumEntriesInLocalRow(i);
        typename Tpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::nonconst_values_host_view_type scaledVals("ScaledVals", nnz);

        for (size_t j = 0; j < nnz; ++j) {
          scaledVals[j] = scalingVector[i] * vals[j];
        }

        if (nnz > 0) {
          tpOp.replaceLocalValues(i, cols, scaledVals);
        }
      }  // for (size_t i=0; ...

    } else {
      typename Tpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::global_inds_host_view_type cols;
      typename Tpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::values_host_view_type vals;

      for (size_t i = 0; i < rowMap->getLocalNumElements(); ++i) {
        GlobalOrdinal gid = rowMap->getGlobalElement(i);
        tpOp.getGlobalRowView(gid, cols, vals);
        size_t nnz = tpOp.getNumEntriesInGlobalRow(gid);
        typename Tpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::nonconst_values_host_view_type scaledVals("ScaledVals", nnz);

        // FIXME FIXME FIXME FIXME FIXME FIXME
        for (size_t j = 0; j < nnz; ++j)
          scaledVals[j] = scalingVector[i] * vals[j];  // FIXME i or gid?

        if (nnz > 0) {
          tpOp.replaceGlobalValues(gid, cols, scaledVals);
        }
      }  // for (size_t i=0; ...
    }

    if (doFillComplete) {
      if (domainMap == Teuchos::null || rangeMap == Teuchos::null)
        throw Exceptions::RuntimeError("In Utilities::Scaling: cannot fillComplete because the domain and/or range map hasn't been defined");

      RCP<Teuchos::ParameterList> params = rcp(new Teuchos::ParameterList());
      params->set("Optimize Storage", doOptimizeStorage);
      params->set("No Nonlocal Changes", true);
      Op.fillComplete(Op.getDomainMap(), Op.getRangeMap(), params);
    }
  } catch (...) {
    throw Exceptions::RuntimeError("Only Tpetra::CrsMatrix types can be scaled (Err.1)");
  }
}  // MyOldScaleMatrix_Tpetra()

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
RCP<Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>>
Utilities<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    Transpose(Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>& Op, bool /* optimizeTranspose */, const std::string& label, const Teuchos::RCP<Teuchos::ParameterList>& params) {
#if defined(HAVE_MUELU_EPETRA) && defined(HAVE_MUELU_EPETRAEXT)
  std::string TorE = "epetra";
#else
  std::string TorE = "tpetra";
#endif

#if defined(HAVE_MUELU_EPETRA) && defined(HAVE_MUELU_EPETRAEXT)
  try {
    const Epetra_CrsMatrix& epetraOp = Utilities<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Op2NonConstEpetraCrs(Op);
    (void)epetraOp;  // silence "unused variable" compiler warning
  } catch (...) {
    TorE = "tpetra";
  }
#endif

  if (TorE == "tpetra") {
    using Helpers = Xpetra::Helpers<Scalar, LocalOrdinal, GlobalOrdinal, Node>;
    /***************************************************************/
    if (Helpers::isTpetraCrs(Op)) {
      const Tpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>& tpetraOp = Utilities<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Op2TpetraCrs(Op);

      RCP<Tpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>> A;
      Tpetra::RowMatrixTransposer<Scalar, LocalOrdinal, GlobalOrdinal, Node> transposer(rcpFromRef(tpetraOp), label);  // more than meets the eye

      {
        using Teuchos::ParameterList;
        using Teuchos::rcp;
        RCP<ParameterList> transposeParams = params.is_null() ? rcp(new ParameterList) : rcp(new ParameterList(*params));
        transposeParams->set("sort", false);
        A = transposer.createTranspose(transposeParams);
      }

      RCP<Xpetra::TpetraCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>> AA = rcp(new Xpetra::TpetraCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>(A));
      RCP<Xpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>> AAA      = rcp_implicit_cast<Xpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>>(AA);
      RCP<Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>> AAAA        = rcp(new Xpetra::CrsMatrixWrap<Scalar, LocalOrdinal, GlobalOrdinal, Node>(AAA));
      if (!AAAA->isFillComplete())
        AAAA->fillComplete(Op.getRangeMap(), Op.getDomainMap());

      if (Op.IsView("stridedMaps"))
        AAAA->CreateView("stridedMaps", Teuchos::rcpFromRef(Op), true /*doTranspose*/);

      return AAAA;
    } else if (Helpers::isTpetraBlockCrs(Op)) {
      using XMatrix        = Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>;
      using XCrsMatrix     = Xpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>;
      using XCrsMatrixWrap = Xpetra::CrsMatrixWrap<Scalar, LocalOrdinal, GlobalOrdinal, Node>;
      using BCRS           = Tpetra::BlockCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>;
      // using CRS  = Tpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>;
      const BCRS& tpetraOp = Utilities<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Op2TpetraBlockCrs(Op);

      RCP<BCRS> At;
      {
        Tpetra::BlockCrsMatrixTransposer<Scalar, LocalOrdinal, GlobalOrdinal, Node> transposer(rcpFromRef(tpetraOp), label);

        using Teuchos::ParameterList;
        using Teuchos::rcp;
        RCP<ParameterList> transposeParams = params.is_null() ? rcp(new ParameterList) : rcp(new ParameterList(*params));
        transposeParams->set("sort", false);
        At = transposer.createTranspose(transposeParams);
      }

      RCP<Xpetra::TpetraBlockCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>> AA = rcp(new Xpetra::TpetraBlockCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>(At));
      RCP<XCrsMatrix> AAA                                                             = rcp_implicit_cast<XCrsMatrix>(AA);
      RCP<XMatrix> AAAA                                                               = rcp(new XCrsMatrixWrap(AAA));

      if (Op.IsView("stridedMaps"))
        AAAA->CreateView("stridedMaps", Teuchos::rcpFromRef(Op), true /*doTranspose*/);

      return AAAA;
    } else {
      throw Exceptions::RuntimeError("Utilities::Transpose failed, perhaps because matrix is not a Crs matrix");
    }
  }  // if

  // Epetra case
  std::cout << "Utilities::Transpose() not implemented for Epetra" << std::endl;
  return Teuchos::null;

}  // Transpose

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
RCP<Xpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>>
Utilities<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    RealValuedToScalarMultiVector(RCP<Xpetra::MultiVector<typename Teuchos::ScalarTraits<Scalar>::coordinateType, LocalOrdinal, GlobalOrdinal, Node>> X) {
  RCP<Xpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>> Xscalar;
#if defined(HAVE_XPETRA_TPETRA) && (defined(HAVE_TPETRA_INST_COMPLEX_DOUBLE) || defined(HAVE_TPETRA_INST_COMPLEX_FLOAT))
  using range_type = Kokkos::RangePolicy<LocalOrdinal, typename Node::execution_space>;

  // Need to cast the real-valued multivector to Scalar=complex
  if ((typeid(Scalar).name() == typeid(std::complex<double>).name()) ||
      (typeid(Scalar).name() == typeid(std::complex<float>).name())) {
    size_t numVecs  = X->getNumVectors();
    Xscalar         = Xpetra::MultiVectorFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Build(X->getMap(), numVecs);
    auto XVec       = X->getDeviceLocalView(Xpetra::Access::ReadOnly);
    auto XVecScalar = Xscalar->getDeviceLocalView(Xpetra::Access::ReadWrite);

    Kokkos::parallel_for(
        "MueLu:Utils::RealValuedToScalarMultiVector", range_type(0, X->getLocalLength()),
        KOKKOS_LAMBDA(const size_t i) {
          for (size_t j = 0; j < numVecs; j++)
            XVecScalar(i, j) = XVec(i, j);
        });
  } else
#endif
    Xscalar = rcp_dynamic_cast<Xpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>>(X);
  return Xscalar;
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
RCP<Xpetra::MultiVector<typename Teuchos::ScalarTraits<Scalar>::magnitudeType, LocalOrdinal, GlobalOrdinal, Node>>
Utilities<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    ExtractCoordinatesFromParameterList(ParameterList& paramList) {
  RCP<Xpetra::MultiVector<typename Teuchos::ScalarTraits<Scalar>::magnitudeType, LocalOrdinal, GlobalOrdinal, Node>> coordinates = Teuchos::null;

  // check whether coordinates are contained in parameter list
  if (paramList.isParameter("Coordinates") == false)
    return coordinates;

    // define Tpetra::MultiVector type with Scalar=float only if
    // * ETI is turned off, since then the compiler will instantiate it automatically OR
    // * Tpetra is instantiated on Scalar=float
#if !defined(HAVE_TPETRA_EXPLICIT_INSTANTIATION) || defined(HAVE_TPETRA_INST_FLOAT)
  typedef Tpetra::MultiVector<float, LocalOrdinal, GlobalOrdinal, Node> tfMV;
  RCP<tfMV> floatCoords = Teuchos::null;
#endif

  // define Tpetra::MultiVector type with Scalar=double only if
  // * ETI is turned off, since then the compiler will instantiate it automatically OR
  // * Tpetra is instantiated on Scalar=double
#if !defined(HAVE_TPETRA_EXPLICIT_INSTANTIATION) || defined(HAVE_TPETRA_INST_DOUBLE)
  typedef Tpetra::MultiVector<typename Teuchos::ScalarTraits<Scalar>::magnitudeType, LocalOrdinal, GlobalOrdinal, Node> tdMV;
  RCP<tdMV> doubleCoords = Teuchos::null;
  if (paramList.isType<RCP<tdMV>>("Coordinates")) {
    // Coordinates are stored as a double vector
    doubleCoords = paramList.get<RCP<tdMV>>("Coordinates");
    paramList.remove("Coordinates");
  }
#if !defined(HAVE_TPETRA_EXPLICIT_INSTANTIATION) || defined(HAVE_TPETRA_INST_FLOAT)
  else if (paramList.isType<RCP<tfMV>>("Coordinates")) {
    // check if coordinates are stored as a float vector
    floatCoords = paramList.get<RCP<tfMV>>("Coordinates");
    paramList.remove("Coordinates");
    doubleCoords = rcp(new tdMV(floatCoords->getMap(), floatCoords->getNumVectors()));
    deep_copy(*doubleCoords, *floatCoords);
  }
#endif
  // We have the coordinates in a Tpetra double vector
  if (doubleCoords != Teuchos::null) {
    // rcp(new Xpetra::TpetraMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>(Vtpetra));
    coordinates = Teuchos::rcp_dynamic_cast<Xpetra::MultiVector<typename Teuchos::ScalarTraits<Scalar>::magnitudeType, LocalOrdinal, GlobalOrdinal, Node>>(MueLu::TpetraMultiVector_To_XpetraMultiVector<typename Teuchos::ScalarTraits<Scalar>::magnitudeType, LocalOrdinal, GlobalOrdinal, Node>(doubleCoords));
    TEUCHOS_TEST_FOR_EXCEPT(Teuchos::is_null(coordinates));
    TEUCHOS_TEST_FOR_EXCEPT(doubleCoords->getNumVectors() != coordinates->getNumVectors());
  }
#else
  // coordinates usually are stored as double vector
  // Tpetra is not instantiated on scalar=double
  throw Exceptions::RuntimeError("ExtractCoordinatesFromParameterList: The coordinates vector in parameter list is expected to be a Tpetra multivector with SC=double or float.");
#endif

  // check for Xpetra coordinates vector
  if (paramList.isType<decltype(coordinates)>("Coordinates")) {
    coordinates = paramList.get<decltype(coordinates)>("Coordinates");
  }

  return coordinates;
}  // ExtractCoordinatesFromParameterList

}  // namespace MueLu

#define MUELU_UTILITIES_SHORT
#endif  // MUELU_UTILITIES_DEF_HPP

//  LocalWords:  LocalOrdinal
