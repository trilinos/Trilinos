#ifndef XPETRA_HELPERS_DEF_HPP
#define XPETRA_HELPERS_DEF_HPP

#include "Xpetra_Helpers_decl.hpp"

namespace Xpetra {

#ifdef HAVE_XPETRA_EPETRA
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
RCP<const Epetra_CrsMatrix> Helpers<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Op2EpetraCrs(RCP<Matrix> Op) {
  // Get the underlying Epetra Mtx
  RCP<const CrsMatrixWrap> crsOp = Teuchos::rcp_dynamic_cast<const CrsMatrixWrap>(Op);
  TEUCHOS_TEST_FOR_EXCEPTION(crsOp == Teuchos::null, Xpetra::Exceptions::BadCast,
                             "Cast from Xpetra::Matrix to Xpetra::CrsMatrixWrap failed");

  RCP<const CrsMatrix> tmp_CrsMtx                        = crsOp->getCrsMatrix();
  const RCP<const EpetraCrsMatrixT<GO, NO>>& tmp_ECrsMtx = Teuchos::rcp_dynamic_cast<const EpetraCrsMatrixT<GO, NO>>(tmp_CrsMtx);
  TEUCHOS_TEST_FOR_EXCEPTION(tmp_ECrsMtx == Teuchos::null, Exceptions::BadCast,
                             "Cast from Xpetra::CrsMatrix to Xpetra::EpetraCrsMatrix failed");

  return tmp_ECrsMtx->getEpetra_CrsMatrix();
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
RCP<Epetra_CrsMatrix> Helpers<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Op2NonConstEpetraCrs(RCP<Matrix> Op) {
  RCP<Epetra_CrsMatrix> A;
  // Get the underlying Epetra Mtx
  RCP<const CrsMatrixWrap> crsOp = Teuchos::rcp_dynamic_cast<const CrsMatrixWrap>(Op);
  TEUCHOS_TEST_FOR_EXCEPTION(crsOp == Teuchos::null, Exceptions::BadCast,
                             "Cast from Xpetra::Matrix to Xpetra::CrsMatrixWrap failed");

  RCP<const CrsMatrix> tmp_CrsMtx                        = crsOp->getCrsMatrix();
  const RCP<const EpetraCrsMatrixT<GO, NO>>& tmp_ECrsMtx = Teuchos::rcp_dynamic_cast<const EpetraCrsMatrixT<GO, NO>>(tmp_CrsMtx);
  TEUCHOS_TEST_FOR_EXCEPTION(tmp_ECrsMtx == Teuchos::null, Exceptions::BadCast, "Cast from Xpetra::CrsMatrix to Xpetra::EpetraCrsMatrix failed");

  return tmp_ECrsMtx->getEpetra_CrsMatrixNonConst();
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
const Epetra_CrsMatrix& Helpers<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Op2EpetraCrs(const Matrix& Op) {
  // Get the underlying Epetra Mtx
  try {
    const CrsMatrixWrap& crsOp                             = dynamic_cast<const CrsMatrixWrap&>(Op);
    RCP<const CrsMatrix> tmp_CrsMtx                        = crsOp.getCrsMatrix();
    const RCP<const EpetraCrsMatrixT<GO, NO>>& tmp_ECrsMtx = Teuchos::rcp_dynamic_cast<const EpetraCrsMatrixT<GO, NO>>(tmp_CrsMtx);
    TEUCHOS_TEST_FOR_EXCEPTION(tmp_ECrsMtx == Teuchos::null, Xpetra::Exceptions::BadCast,
                               "Cast from Xpetra::CrsMatrix to Xpetra::EpetraCrsMatrix failed");

    return *tmp_ECrsMtx->getEpetra_CrsMatrix();

  } catch (...) {
    throw(Xpetra::Exceptions::BadCast("Cast from Xpetra::Matrix to Xpetra::CrsMatrixWrap failed"));
  }
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
Epetra_CrsMatrix& Helpers<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Op2NonConstEpetraCrs(const Matrix& Op) {
  RCP<Epetra_CrsMatrix> A;
  // Get the underlying Epetra Mtx
  try {
    const CrsMatrixWrap& crsOp                             = dynamic_cast<const CrsMatrixWrap&>(Op);
    RCP<const CrsMatrix> tmp_CrsMtx                        = crsOp.getCrsMatrix();
    const RCP<const EpetraCrsMatrixT<GO, NO>>& tmp_ECrsMtx = Teuchos::rcp_dynamic_cast<const EpetraCrsMatrixT<GO, NO>>(tmp_CrsMtx);
    TEUCHOS_TEST_FOR_EXCEPTION(tmp_ECrsMtx == Teuchos::null, Xpetra::Exceptions::BadCast, "Cast from Xpetra::CrsMatrix to Xpetra::EpetraCrsMatrix failed");

    return *Teuchos::rcp_const_cast<Epetra_CrsMatrix>(tmp_ECrsMtx->getEpetra_CrsMatrix());

  } catch (...) {
    throw(Xpetra::Exceptions::BadCast("Cast from Xpetra::Matrix to Xpetra::CrsMatrixWrap failed"));
  }
}
#endif  // HAVE_XPETRA_EPETRA

#ifdef HAVE_XPETRA_TPETRA
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
RCP<const Tpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>> Helpers<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Op2TpetraCrs(RCP<Matrix> Op) {
  // Get the underlying Tpetra Mtx
  RCP<const CrsMatrixWrap> crsOp = Teuchos::rcp_dynamic_cast<const CrsMatrixWrap>(Op);
  TEUCHOS_TEST_FOR_EXCEPTION(crsOp == Teuchos::null, Xpetra::Exceptions::BadCast, "Cast from Xpetra::Matrix to Xpetra::CrsMatrixWrap failed");

  RCP<const CrsMatrix> tmp_CrsMtx                                       = crsOp->getCrsMatrix();
  const RCP<const Xpetra::TpetraCrsMatrix<SC, LO, GO, NO>>& tmp_ECrsMtx = Teuchos::rcp_dynamic_cast<const Xpetra::TpetraCrsMatrix<SC, LO, GO, NO>>(tmp_CrsMtx);
  TEUCHOS_TEST_FOR_EXCEPTION(tmp_ECrsMtx == Teuchos::null, Xpetra::Exceptions::BadCast, "Cast from Xpetra::CrsMatrix to Xpetra::TpetraCrsMatrix failed");

  return tmp_ECrsMtx->getTpetra_CrsMatrix();
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
RCP<Tpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>> Helpers<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Op2NonConstTpetraCrs(RCP<Matrix> Op) {
  // Get the underlying Tpetra Mtx
  RCP<const CrsMatrixWrap> crsOp = Teuchos::rcp_dynamic_cast<const CrsMatrixWrap>(Op);
  TEUCHOS_TEST_FOR_EXCEPTION(crsOp == Teuchos::null, Xpetra::Exceptions::BadCast, "Cast from Xpetra::Matrix to Xpetra::CrsMatrixWrap failed");
  RCP<const CrsMatrix> tmp_CrsMtx                                       = crsOp->getCrsMatrix();
  const RCP<const Xpetra::TpetraCrsMatrix<SC, LO, GO, NO>>& tmp_ECrsMtx = Teuchos::rcp_dynamic_cast<const Xpetra::TpetraCrsMatrix<SC, LO, GO, NO>>(tmp_CrsMtx);
  TEUCHOS_TEST_FOR_EXCEPTION(tmp_ECrsMtx == Teuchos::null, Xpetra::Exceptions::BadCast, "Cast from Xpetra::CrsMatrix to Xpetra::TpetraCrsMatrix failed");

  return tmp_ECrsMtx->getTpetra_CrsMatrixNonConst();
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
const Tpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>& Helpers<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Op2TpetraCrs(const Matrix& Op) {
  // Get the underlying Tpetra Mtx
  try {
    const CrsMatrixWrap& crsOp                                            = dynamic_cast<const CrsMatrixWrap&>(Op);
    RCP<const CrsMatrix> tmp_CrsMtx                                       = crsOp.getCrsMatrix();
    const RCP<const Xpetra::TpetraCrsMatrix<SC, LO, GO, NO>>& tmp_TCrsMtx = Teuchos::rcp_dynamic_cast<const Xpetra::TpetraCrsMatrix<SC, LO, GO, NO>>(tmp_CrsMtx);
    TEUCHOS_TEST_FOR_EXCEPTION(tmp_TCrsMtx == Teuchos::null, Xpetra::Exceptions::BadCast, "Cast from Xpetra::CrsMatrix to Xpetra::TpetraCrsMatrix failed");

    return *tmp_TCrsMtx->getTpetra_CrsMatrix();

  } catch (...) {
    throw(Xpetra::Exceptions::BadCast("Cast from Xpetra::Matrix to Xpetra::CrsMatrixWrap failed"));
  }
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
Tpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>& Helpers<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Op2NonConstTpetraCrs(const Matrix& Op) {
  // Get the underlying Tpetra Mtx
  try {
    const CrsMatrixWrap& crsOp                                            = dynamic_cast<const CrsMatrixWrap&>(Op);
    RCP<const CrsMatrix> tmp_CrsMtx                                       = crsOp.getCrsMatrix();
    const RCP<const Xpetra::TpetraCrsMatrix<SC, LO, GO, NO>>& tmp_TCrsMtx = Teuchos::rcp_dynamic_cast<const Xpetra::TpetraCrsMatrix<SC, LO, GO, NO>>(tmp_CrsMtx);
    TEUCHOS_TEST_FOR_EXCEPTION(tmp_TCrsMtx == Teuchos::null, Xpetra::Exceptions::BadCast, "Cast from Xpetra::CrsMatrix to Xpetra::TpetraCrsMatrix failed");

    return *Teuchos::rcp_const_cast<Tpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>>(tmp_TCrsMtx->getTpetra_CrsMatrix());

  } catch (...) {
    throw(Xpetra::Exceptions::BadCast("Cast from Xpetra::Matrix to Xpetra::CrsMatrixWrap failed"));
  }
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
bool Helpers<Scalar, LocalOrdinal, GlobalOrdinal, Node>::isTpetraCrs(RCP<Matrix> Op) {
  RCP<const CrsMatrixWrap> crsOp = Teuchos::rcp_dynamic_cast<const CrsMatrixWrap>(Op);
  if (crsOp == Teuchos::null) return false;
  RCP<const CrsMatrix> tmp_CrsMtx                                       = crsOp->getCrsMatrix();
  const RCP<const Xpetra::TpetraCrsMatrix<SC, LO, GO, NO>>& tmp_ECrsMtx = Teuchos::rcp_dynamic_cast<const Xpetra::TpetraCrsMatrix<SC, LO, GO, NO>>(tmp_CrsMtx);
  if (tmp_ECrsMtx == Teuchos::null)
    return false;
  else
    return true;
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
bool Helpers<Scalar, LocalOrdinal, GlobalOrdinal, Node>::isTpetraCrs(const Matrix& Op) {
  try {
    const CrsMatrixWrap& crsOp                                            = dynamic_cast<const CrsMatrixWrap&>(Op);
    RCP<const CrsMatrix> tmp_CrsMtx                                       = crsOp.getCrsMatrix();
    const RCP<const Xpetra::TpetraCrsMatrix<SC, LO, GO, NO>>& tmp_ECrsMtx = Teuchos::rcp_dynamic_cast<const Xpetra::TpetraCrsMatrix<SC, LO, GO, NO>>(tmp_CrsMtx);
    if (tmp_ECrsMtx == Teuchos::null)
      return false;
    else
      return true;
  } catch (...) {
    return false;
  }
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
RCP<const Tpetra::BlockCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>> Helpers<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Op2TpetraBlockCrs(RCP<Matrix> Op) {
  RCP<const CrsMatrixWrap> crsOp = Teuchos::rcp_dynamic_cast<const CrsMatrixWrap>(Op);
  TEUCHOS_TEST_FOR_EXCEPTION(crsOp == Teuchos::null, Xpetra::Exceptions::BadCast, "Cast from Xpetra::Matrix to Xpetra::CrsMatrixWrap failed");

  RCP<const CrsMatrix> tmp_CrsMtx              = crsOp->getCrsMatrix();
  RCP<const TpetraBlockCrsMatrix> tmp_BlockCrs = Teuchos::rcp_dynamic_cast<const TpetraBlockCrsMatrix>(tmp_CrsMtx);
  TEUCHOS_TEST_FOR_EXCEPTION(tmp_BlockCrs == Teuchos::null, Xpetra::Exceptions::BadCast, "Cast from Xpetra::CrsMatrix to Xpetra::TpetraBlockCrsMatrix failed");
  return tmp_BlockCrs->getTpetra_BlockCrsMatrix();
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
RCP<Tpetra::BlockCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>> Helpers<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Op2NonConstTpetraBlockCrs(RCP<Matrix> Op) {
  RCP<const CrsMatrixWrap> crsOp = Teuchos::rcp_dynamic_cast<const CrsMatrixWrap>(Op);
  TEUCHOS_TEST_FOR_EXCEPTION(crsOp == Teuchos::null, Xpetra::Exceptions::BadCast, "Cast from Xpetra::Matrix to Xpetra::CrsMatrixWrap failed");

  RCP<const CrsMatrix> tmp_CrsMtx              = crsOp->getCrsMatrix();
  RCP<const TpetraBlockCrsMatrix> tmp_BlockCrs = Teuchos::rcp_dynamic_cast<const TpetraBlockCrsMatrix>(tmp_CrsMtx);
  TEUCHOS_TEST_FOR_EXCEPTION(tmp_BlockCrs == Teuchos::null, Xpetra::Exceptions::BadCast, "Cast from Xpetra::CrsMatrix to Xpetra::TpetraBlockCrsMatrix failed");
  return tmp_BlockCrs->getTpetra_BlockCrsMatrixNonConst();
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
const Tpetra::BlockCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>& Helpers<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Op2TpetraBlockCrs(const Matrix& Op) {
  try {
    const CrsMatrixWrap& crsOp                   = dynamic_cast<const CrsMatrixWrap&>(Op);
    RCP<const CrsMatrix> tmp_CrsMtx              = crsOp.getCrsMatrix();
    RCP<const TpetraBlockCrsMatrix> tmp_BlockCrs = Teuchos::rcp_dynamic_cast<const TpetraBlockCrsMatrix>(tmp_CrsMtx);
    TEUCHOS_TEST_FOR_EXCEPTION(tmp_BlockCrs == Teuchos::null, Xpetra::Exceptions::BadCast, "Cast from Xpetra::CrsMatrix to Xpetra::TpetraBlockCrsMatrix failed");
    return *tmp_BlockCrs->getTpetra_BlockCrsMatrix();
  } catch (...) {
    throw(Xpetra::Exceptions::BadCast("Cast from Xpetra::Matrix to Xpetra::CrsMatrixWrap failed"));
  }
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
Tpetra::BlockCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>& Helpers<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Op2NonConstTpetraBlockCrs(const Matrix& Op) {
  try {
    const CrsMatrixWrap& crsOp                   = dynamic_cast<const CrsMatrixWrap&>(Op);
    RCP<const CrsMatrix> tmp_CrsMtx              = crsOp.getCrsMatrix();
    RCP<const TpetraBlockCrsMatrix> tmp_BlockCrs = Teuchos::rcp_dynamic_cast<const TpetraBlockCrsMatrix>(tmp_CrsMtx);
    TEUCHOS_TEST_FOR_EXCEPTION(tmp_BlockCrs == Teuchos::null, Xpetra::Exceptions::BadCast, "Cast from Xpetra::CrsMatrix to Xpetra::TpetraBlockCrsMatrix failed");
    return *tmp_BlockCrs->getTpetra_BlockCrsMatrixNonConst();
  } catch (...) {
    throw(Xpetra::Exceptions::BadCast("Cast from Xpetra::Matrix to Xpetra::CrsMatrixWrap failed"));
  }
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
bool Helpers<Scalar, LocalOrdinal, GlobalOrdinal, Node>::isTpetraBlockCrs(RCP<Matrix> Op) {
  RCP<const CrsMatrixWrap> crsOp = Teuchos::rcp_dynamic_cast<const CrsMatrixWrap>(Op);
  if (crsOp == Teuchos::null) return false;
  RCP<const CrsMatrix> tmp_CrsMtx              = crsOp->getCrsMatrix();
  RCP<const TpetraBlockCrsMatrix> tmp_BlockCrs = Teuchos::rcp_dynamic_cast<const TpetraBlockCrsMatrix>(tmp_CrsMtx);
  if (tmp_BlockCrs == Teuchos::null)
    return false;
  else
    return true;
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
bool Helpers<Scalar, LocalOrdinal, GlobalOrdinal, Node>::isTpetraBlockCrs(const Matrix& Op) {
  try {
    const CrsMatrixWrap& crsOp                   = dynamic_cast<const CrsMatrixWrap&>(Op);
    RCP<const CrsMatrix> tmp_CrsMtx              = crsOp.getCrsMatrix();
    RCP<const TpetraBlockCrsMatrix> tmp_BlockCrs = Teuchos::rcp_dynamic_cast<const TpetraBlockCrsMatrix>(tmp_CrsMtx);
    if (tmp_BlockCrs == Teuchos::null)
      return false;
    else
      return true;
  } catch (...) {
    return false;
  }
}
#else  // HAVE_XPETRA_TPETRA
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
bool Helpers<Scalar, LocalOrdinal, GlobalOrdinal, Node>::isTpetraCrs(const Matrix& Op) {
  return false;
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
bool Helpers<Scalar, LocalOrdinal, GlobalOrdinal, Node>::isTpetraBlockCrs(const Matrix& Op) {
  return false;
}

#endif  // HAVE_XPETRA_TPETRA

#ifdef HAVE_XPETRA_TPETRA
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
Teuchos::RCP<Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>> Helpers<Scalar, LocalOrdinal, GlobalOrdinal, Node>::tpetraAdd(
    const tcrs_matrix_type& A, bool transposeA, const typename tcrs_matrix_type::scalar_type alpha,
    const tcrs_matrix_type& B, bool transposeB, const typename tcrs_matrix_type::scalar_type beta) {
  using Teuchos::Array;
  using Teuchos::ParameterList;
  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::rcp_implicit_cast;
  using Teuchos::rcpFromRef;
  using XTCrsType                    = Xpetra::TpetraCrsMatrix<SC, LO, GO, NO>;
  using CrsType                      = Xpetra::CrsMatrix<SC, LO, GO, NO>;
  using CrsWrap                      = Xpetra::CrsMatrixWrap<SC, LO, GO, NO>;
  using transposer_type              = Tpetra::RowMatrixTransposer<SC, LO, GO, NO>;
  using import_type                  = Tpetra::Import<LO, GO, NO>;
  RCP<const tcrs_matrix_type> Aprime = rcpFromRef(A);
  if (transposeA)
    Aprime = transposer_type(Aprime).createTranspose();
  // Decide whether the fast code path can be taken.
  if (A.isFillComplete() && B.isFillComplete()) {
    RCP<tcrs_matrix_type> C      = rcp(new tcrs_matrix_type(Aprime->getRowMap(), 0));
    RCP<ParameterList> addParams = rcp(new ParameterList);
    addParams->set("Call fillComplete", false);
    // passing A after B means C will have the same domain/range map as A (or A^T if transposeA)
    Tpetra::MatrixMatrix::add<SC, LO, GO, NO>(beta, transposeB, B, alpha, false, *Aprime, *C, Teuchos::null, Teuchos::null, addParams);
    return rcp_implicit_cast<Matrix>(rcp(new CrsWrap(rcp_implicit_cast<CrsType>(rcp(new XTCrsType(C))))));
  } else {
    // Slow case - one or both operands are non-fill complete.
    // TODO: deprecate this.
    // Need to compute the explicit transpose before add if transposeA and/or transposeB.
    // This is to count the number of entries to allocate per row in the final sum.
    RCP<const tcrs_matrix_type> Bprime = rcpFromRef(B);
    if (transposeB)
      Bprime = transposer_type(Bprime).createTranspose();
    // Use Aprime's row map as C's row map.
    if (!(Aprime->getRowMap()->isSameAs(*(Bprime->getRowMap())))) {
      auto import = rcp(new import_type(Bprime->getRowMap(), Aprime->getRowMap()));
      Bprime      = Tpetra::importAndFillCompleteCrsMatrix<tcrs_matrix_type>(Bprime, *import, Aprime->getDomainMap(), Aprime->getRangeMap());
    }
    // Count the entries to allocate for C in each local row.
    LO numLocalRows = Aprime->getLocalNumRows();
    Array<size_t> allocPerRow(numLocalRows);
    // 0 is always the minimum LID
    for (LO i = 0; i < numLocalRows; i++) {
      allocPerRow[i] = Aprime->getNumEntriesInLocalRow(i) + Bprime->getNumEntriesInLocalRow(i);
    }
    // Construct C
    RCP<tcrs_matrix_type> C = rcp(new tcrs_matrix_type(Aprime->getRowMap(), allocPerRow()));
    // Compute the sum in C (already took care of transposes)
    Tpetra::MatrixMatrix::Add<SC, LO, GO, NO>(
        *Aprime, false, alpha,
        *Bprime, false, beta,
        C);
    return rcp(new CrsWrap(rcp_implicit_cast<CrsType>(rcp(new XTCrsType(C)))));
  }
}
#endif

#ifdef HAVE_XPETRA_EPETRAEXT
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void Helpers<Scalar, LocalOrdinal, GlobalOrdinal, Node>::epetraExtMult(const Matrix& A, bool transposeA, const Matrix& B, bool transposeB, Matrix& C, bool fillCompleteResult) {
  Epetra_CrsMatrix& epA = Op2NonConstEpetraCrs(A);
  Epetra_CrsMatrix& epB = Op2NonConstEpetraCrs(B);
  Epetra_CrsMatrix& epC = Op2NonConstEpetraCrs(C);
  // Want a static profile (possibly fill complete) matrix as the result.
  // But, EpetraExt Multiply needs C to be dynamic profile, so compute the product in a temporary DynamicProfile matrix.
  const Epetra_Map& Crowmap = epC.RowMap();
  int errCode               = 0;
  Epetra_CrsMatrix Ctemp(::Copy, Crowmap, 0, false);
  if (fillCompleteResult) {
    errCode = EpetraExt::MatrixMatrix::Multiply(epA, transposeA, epB, transposeB, Ctemp, true);
    if (!errCode) {
      epC = Ctemp;
      C.fillComplete();
    }
  } else {
    errCode = EpetraExt::MatrixMatrix::Multiply(epA, transposeA, epB, transposeB, Ctemp, false);
    if (!errCode) {
      int numLocalRows             = Crowmap.NumMyElements();
      long long* globalElementList = nullptr;
      Crowmap.MyGlobalElementsPtr(globalElementList);
      Teuchos::Array<int> entriesPerRow(numLocalRows, 0);
      for (int i = 0; i < numLocalRows; i++) {
        entriesPerRow[i] = Ctemp.NumGlobalEntries(globalElementList[i]);
      }
      // know how many entries to allocate in epC (which must be static profile)
      epC = Epetra_CrsMatrix(::Copy, Crowmap, entriesPerRow.data(), true);
      for (int i = 0; i < numLocalRows; i++) {
        int gid = globalElementList[i];
        int numEntries;
        double* values;
        int* inds;
        Ctemp.ExtractGlobalRowView(gid, numEntries, values, inds);
        epC.InsertGlobalValues(gid, numEntries, values, inds);
      }
    }
  }
  if (errCode) {
    std::ostringstream buf;
    buf << errCode;
    std::string msg = "EpetraExt::MatrixMatrix::Multiply returned nonzero error code " + buf.str();
    throw(Exceptions::RuntimeError(msg));
  }
}
#endif

}  // namespace Xpetra

#endif
