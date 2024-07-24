#ifndef XPETRA_HELPERS_DECL_HPP
#define XPETRA_HELPERS_DECL_HPP

#include "Xpetra_ConfigDefs.hpp"

#include "Xpetra_BlockedCrsMatrix.hpp"
#include "Xpetra_CrsMatrixWrap.hpp"
#include "Xpetra_MapExtractor.hpp"
#include "Xpetra_Map.hpp"
#include "Xpetra_Matrix.hpp"

#ifdef HAVE_XPETRA_EPETRA
#include <Xpetra_EpetraCrsMatrix_fwd.hpp>
#endif

#ifdef HAVE_XPETRA_EPETRAEXT
#include <EpetraExt_MatrixMatrix.h>
#include <EpetraExt_RowMatrixOut.h>
#include <Epetra_RowMatrixTransposer.h>
#endif  // HAVE_XPETRA_EPETRAEXT

#ifdef HAVE_XPETRA_TPETRA
#include <Tpetra_CrsMatrix.hpp>
#include <Tpetra_BlockCrsMatrix.hpp>
#include <TpetraExt_MatrixMatrix.hpp>
#include <Tpetra_RowMatrixTransposer.hpp>
#include <MatrixMarket_Tpetra.hpp>
#include <Xpetra_TpetraCrsMatrix.hpp>
#include <Xpetra_TpetraBlockCrsMatrix.hpp>
#include <Tpetra_BlockCrsMatrix_Helpers.hpp>
#include <Xpetra_TpetraMultiVector.hpp>
#include <Xpetra_TpetraVector.hpp>
#endif  // HAVE_XPETRA_TPETRA

namespace Xpetra {

/*!
    @class Helpers
    @brief Xpetra utility class containing transformation routines between Xpetra::Matrix and Epetra/Tpetra objects

Note: this class is not in the Xpetra_UseShortNames.hpp
*/
template <class Scalar,
          class LocalOrdinal  = int,
          class GlobalOrdinal = LocalOrdinal,
          class Node          = Tpetra::KokkosClassic::DefaultNode::DefaultNodeType>
class Helpers {
#include "Xpetra_UseShortNames.hpp"

 public:
#ifdef HAVE_XPETRA_EPETRA
  static RCP<const Epetra_CrsMatrix> Op2EpetraCrs(RCP<Matrix> Op);

  static RCP<Epetra_CrsMatrix> Op2NonConstEpetraCrs(RCP<Matrix> Op);

  static const Epetra_CrsMatrix& Op2EpetraCrs(const Matrix& Op);

  static Epetra_CrsMatrix& Op2NonConstEpetraCrs(const Matrix& Op);
#endif  // HAVE_XPETRA_EPETRA

#ifdef HAVE_XPETRA_TPETRA
  static RCP<const Tpetra::CrsMatrix<SC, LO, GO, NO> > Op2TpetraCrs(RCP<Matrix> Op);

  static RCP<Tpetra::CrsMatrix<SC, LO, GO, NO> > Op2NonConstTpetraCrs(RCP<Matrix> Op);

  static const Tpetra::CrsMatrix<SC, LO, GO, NO>& Op2TpetraCrs(const Matrix& Op);

  static Tpetra::CrsMatrix<SC, LO, GO, NO>& Op2NonConstTpetraCrs(const Matrix& Op);

  static bool isTpetraCrs(RCP<Matrix> Op);

  static bool isTpetraCrs(const Matrix& Op);

  static RCP<const Tpetra::BlockCrsMatrix<SC, LO, GO, NO> > Op2TpetraBlockCrs(RCP<Matrix> Op);

  static RCP<Tpetra::BlockCrsMatrix<SC, LO, GO, NO> > Op2NonConstTpetraBlockCrs(RCP<Matrix> Op);

  static const Tpetra::BlockCrsMatrix<SC, LO, GO, NO>& Op2TpetraBlockCrs(const Matrix& Op);

  static Tpetra::BlockCrsMatrix<SC, LO, GO, NO>& Op2NonConstTpetraBlockCrs(const Matrix& Op);

  static bool isTpetraBlockCrs(RCP<Matrix> Op);

  static bool isTpetraBlockCrs(const Matrix& Op);
#else  // HAVE_XPETRA_TPETRA
  static bool isTpetraCrs(const Matrix& Op);

  static bool isTpetraBlockCrs(const Matrix& Op);

#endif  // HAVE_XPETRA_TPETRA

#ifdef HAVE_XPETRA_TPETRA
  using tcrs_matrix_type = Tpetra::CrsMatrix<SC, LO, GO, NO>;
  static Teuchos::RCP<Matrix> tpetraAdd(
      const tcrs_matrix_type& A, bool transposeA, const typename tcrs_matrix_type::scalar_type alpha,
      const tcrs_matrix_type& B, bool transposeB, const typename tcrs_matrix_type::scalar_type beta);
#endif

#ifdef HAVE_XPETRA_EPETRAEXT
  static void epetraExtMult(const Matrix& A, bool transposeA, const Matrix& B, bool transposeB, Matrix& C, bool fillCompleteResult);
#endif
};
}  // namespace Xpetra

#endif
