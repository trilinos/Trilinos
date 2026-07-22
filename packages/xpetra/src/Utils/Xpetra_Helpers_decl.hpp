// @HEADER
// *****************************************************************************
//             Xpetra: A linear algebra interface package
//
// Copyright 2012 NTESS and the Xpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef XPETRA_HELPERS_DECL_HPP
#define XPETRA_HELPERS_DECL_HPP

#include "Xpetra_ConfigDefs.hpp"

#include "Xpetra_BlockedCrsMatrix.hpp"
#include "Xpetra_CrsMatrixWrap.hpp"
#include "Xpetra_MapExtractor.hpp"
#include "Xpetra_Map.hpp"
#include "Xpetra_Matrix.hpp"

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

  using tcrs_matrix_type = Tpetra::CrsMatrix<SC, LO, GO, NO>;
  static Teuchos::RCP<Matrix> tpetraAdd(
      const tcrs_matrix_type& A, bool transposeA, const typename tcrs_matrix_type::scalar_type alpha,
      const tcrs_matrix_type& B, bool transposeB, const typename tcrs_matrix_type::scalar_type beta);
};
}  // namespace Xpetra

#endif
