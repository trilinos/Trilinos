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

#ifdef HAVE_XPETRA_EPETRA
#include "Epetra_Util.h"
#endif

#ifdef HAVE_XPETRA_TPETRA
#include "Tpetra_Import_Util2.hpp"
#endif

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
    if (lib == Xpetra::UseEpetra) {
#if defined(HAVE_XPETRA_EPETRA)
      throw(Xpetra::Exceptions::RuntimeError("Xpetra::CrsMatrixUtils::sortCrsEntries only available for GO=int or GO=long long with EpetraNode (Serial or OpenMP depending on configuration)"));
#endif  // HAVE_XPETRA_EPETRA
    } else if (lib == Xpetra::UseTpetra) {
#ifdef HAVE_XPETRA_TPETRA
      Tpetra::Import_Util::sortCrsEntries(CRS_rowptr, CRS_colind, CRS_vals);
#endif  // HAVE_XPETRA_TPETRA
    }

    return;
  }

  /// \brief Sort and merge the entries of the (raw CSR) matrix by column index
  ///   within each row.
  static void sortAndMergeCrsEntries(const Teuchos::ArrayView<size_t>& CRS_rowptr,
                                     const Teuchos::ArrayView<LocalOrdinal>& CRS_colind,
                                     const Teuchos::ArrayView<Scalar>& CRS_vals,
                                     const UnderlyingLib lib) {
    if (lib == Xpetra::UseEpetra) {
#if defined(HAVE_XPETRA_EPETRA)
      throw(Xpetra::Exceptions::RuntimeError("Xpetra::CrsMatrixUtils::sortAndMergeCrsEntries only available for GO=int or GO=long long with EpetraNode (Serial or OpenMP depending on configuration)"));
#endif  // HAVE_XPETRA_EPETRA
    } else if (lib == Xpetra::UseTpetra) {
#ifdef HAVE_XPETRA_TPETRA
      Tpetra::Import_Util::sortAndMergeCrsEntries(CRS_rowptr, CRS_colind, CRS_vals);
#endif  // HAVE_XPETRA_TPETRA
    }

    return;
  }

};  // end class CrsMatrixUtils

#ifdef HAVE_XPETRA_EPETRA
// Specialization for double, int, int, EpetraNode
template <>
class CrsMatrixUtils<double, int, int, EpetraNode> {
  typedef double Scalar;
  typedef int LocalOrdinal;
  typedef int GlobalOrdinal;
  typedef EpetraNode Node;
#undef XPETRA_CRSMATRIXUTILS_SHORT

 public:
  /// \brief Sort the entries of the (raw CSR) matrix by column index
  ///   within each row.
  static void sortCrsEntries(const Teuchos::ArrayView<size_t>& CRS_rowptr,
                             const Teuchos::ArrayView<LocalOrdinal>& CRS_colind,
                             const Teuchos::ArrayView<Scalar>& CRS_vals,
                             const UnderlyingLib lib) {
    if (lib == Xpetra::UseEpetra) {
#if defined(HAVE_XPETRA_EPETRA)
      int rv = Epetra_Util::SortCrsEntries(Teuchos::as<int>(CRS_rowptr.size() - 1),
                                           CRS_rowptr.getRawPtr(),
                                           CRS_colind.getRawPtr(),
                                           CRS_vals.getRawPtr());
      if (rv != 0) {
        throw Exceptions::RuntimeError("Epetra_Util::SortCrsEntries() return value of " + Teuchos::toString(rv));
      }
#endif  // HAVE_XPETRA_EPETRA
    } else if (lib == Xpetra::UseTpetra) {
#ifdef HAVE_XPETRA_TPETRA
      Tpetra::Import_Util::sortCrsEntries(CRS_rowptr, CRS_colind, CRS_vals);
#endif  // HAVE_XPETRA_TPETRA
    }

    return;
  }

  /// \brief Sort and merge the entries of the (raw CSR) matrix by column index
  ///   within each row.
  static void sortAndMergeCrsEntries(const Teuchos::ArrayView<size_t>& CRS_rowptr,
                                     const Teuchos::ArrayView<LocalOrdinal>& CRS_colind,
                                     const Teuchos::ArrayView<Scalar>& CRS_vals,
                                     const UnderlyingLib lib) {
    if (lib == Xpetra::UseEpetra) {
#if defined(HAVE_XPETRA_EPETRA)
      int rv = Epetra_Util::SortAndMergeCrsEntries(Teuchos::as<int>(CRS_rowptr.size() - 1),
                                                   CRS_rowptr.getRawPtr(),
                                                   CRS_colind.getRawPtr(),
                                                   CRS_vals.getRawPtr());
      if (rv != 0) {
        throw Exceptions::RuntimeError("Epetra_Util::SortCrsEntries() return value of " + Teuchos::toString(rv));
      }
#endif  // HAVE_XPETRA_EPETRA
    } else if (lib == Xpetra::UseTpetra) {
#ifdef HAVE_XPETRA_TPETRA
      Tpetra::Import_Util::sortAndMergeCrsEntries(CRS_rowptr, CRS_colind, CRS_vals);
#endif  // HAVE_XPETRA_TPETRA
    }

    return;
  }

};  // end class CrsMatrixUtils

// Specialization for double, int, long long, EpetraNode
template <>
class CrsMatrixUtils<double, int, long long, EpetraNode> {
  typedef double Scalar;
  typedef int LocalOrdinal;
  typedef long long GlobalOrdinal;
  typedef EpetraNode Node;
#undef XPETRA_CRSMATRIXUTILS_SHORT

 public:
  /// \brief Sort the entries of the (raw CSR) matrix by column index
  ///   within each row.
  static void sortCrsEntries(const Teuchos::ArrayView<size_t>& CRS_rowptr,
                             const Teuchos::ArrayView<LocalOrdinal>& CRS_colind,
                             const Teuchos::ArrayView<Scalar>& CRS_vals,
                             const UnderlyingLib lib) {
    if (lib == Xpetra::UseEpetra) {
#if defined(HAVE_XPETRA_EPETRA)
      int rv = Epetra_Util::SortCrsEntries(Teuchos::as<int>(CRS_rowptr.size() - 1),
                                           CRS_rowptr.getRawPtr(),
                                           CRS_colind.getRawPtr(),
                                           CRS_vals.getRawPtr());
      if (rv != 0) {
        throw Exceptions::RuntimeError("Epetra_Util::SortCrsEntries() return value of " + Teuchos::toString(rv));
      }
#endif  // HAVE_XPETRA_EPETRA
    } else if (lib == Xpetra::UseTpetra) {
#ifdef HAVE_XPETRA_TPETRA
      Tpetra::Import_Util::sortCrsEntries(CRS_rowptr, CRS_colind, CRS_vals);
#endif  // HAVE_XPETRA_TPETRA
    }

    return;
  }

  /// \brief Sort and merge the entries of the (raw CSR) matrix by column index
  ///   within each row.
  static void sortAndMergeCrsEntries(const Teuchos::ArrayView<size_t>& CRS_rowptr,
                                     const Teuchos::ArrayView<LocalOrdinal>& CRS_colind,
                                     const Teuchos::ArrayView<Scalar>& CRS_vals,
                                     const UnderlyingLib lib) {
    if (lib == Xpetra::UseEpetra) {
#if defined(HAVE_XPETRA_EPETRA)
      int rv = Epetra_Util::SortAndMergeCrsEntries(Teuchos::as<int>(CRS_rowptr.size() - 1),
                                                   CRS_rowptr.getRawPtr(),
                                                   CRS_colind.getRawPtr(),
                                                   CRS_vals.getRawPtr());
      if (rv != 0) {
        throw Exceptions::RuntimeError("Epetra_Util::SortCrsEntries() return value of " + Teuchos::toString(rv));
      }
#endif  // HAVE_XPETRA_EPETRA
    } else if (lib == Xpetra::UseTpetra) {
#ifdef HAVE_XPETRA_TPETRA
      Tpetra::Import_Util::sortAndMergeCrsEntries(CRS_rowptr, CRS_colind, CRS_vals);
#endif  // HAVE_XPETRA_TPETRA
    }

    return;
  }

};      // end class CrsMatrixUtils
#endif  // HAVE_XPETRA_EPETRA for Epetra scpecialization

}  // end namespace Xpetra

#define XPETRA_CRSMATRIXUTILS_SHORT

#endif  // PACKAGES_XPETRA_CRSMATRIX_UTILS_HPP_
