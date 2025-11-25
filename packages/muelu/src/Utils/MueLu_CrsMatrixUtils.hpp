// @HEADER
// *****************************************************************************
//             Xpetra: A linear algebra interface package
//
// Copyright 2012 NTESS and the Xpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef PACKAGES_XPETRA_CRSMATRIX_UTILS_HPP_
#define PACKAGES_XPETRA_CRSMATRIX_UTILS_HPP_

#include "Xpetra_ConfigDefs.hpp"
#include "Xpetra_Exceptions.hpp"
#include "Xpetra_Map.hpp"  // definition of UnderlyingLib

#include "Tpetra_Import_Util2.hpp"

namespace Xpetra {

/*!
  @class CrsMatrixUtils
  @brief Xpetra utility class for CrsMatrix-related routines

  The routines should be independent from Epetra/Tpetra and be purely implemented in Xpetra.

*/
template <class Scalar,
          class LocalOrdinal,
          class GlobalOrdinal,
          class Node = Tpetra::KokkosClassic::DefaultNode::DefaultNodeType>
class CrsMatrixUtils {
#undef XPETRA_CRSMATRIXUTILS_SHORT

 public:
  /// \brief Sort the entries of the (raw CSR) matrix by column index
  ///   within each row.
  static void sortCrsEntries(const Teuchos::ArrayView<size_t>& CRS_rowptr,
                             const Teuchos::ArrayView<LocalOrdinal>& CRS_colind,
                             const Teuchos::ArrayView<Scalar>& CRS_vals,
                             const UnderlyingLib lib) {
    if (lib == Xpetra::UseTpetra) {
      Tpetra::Import_Util::sortCrsEntries(CRS_rowptr, CRS_colind, CRS_vals);
    }

    return;
  }

  /// \brief Sort and merge the entries of the (raw CSR) matrix by column index
  ///   within each row.
  static void sortAndMergeCrsEntries(const Teuchos::ArrayView<size_t>& CRS_rowptr,
                                     const Teuchos::ArrayView<LocalOrdinal>& CRS_colind,
                                     const Teuchos::ArrayView<Scalar>& CRS_vals,
                                     const UnderlyingLib lib) {
    if (lib == Xpetra::UseTpetra) {
      Tpetra::Import_Util::sortAndMergeCrsEntries(CRS_rowptr, CRS_colind, CRS_vals);
    }

    return;
  }

};  // end class CrsMatrixUtils

}  // end namespace Xpetra

#define XPETRA_CRSMATRIXUTILS_SHORT

#endif  // PACKAGES_XPETRA_CRSMATRIX_UTILS_HPP_
