// @HEADER
// *****************************************************************************
//                           Stokhos Package
//
// Copyright 2009 NTESS and the Stokhos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef TPETRA_CRSMATRIX_UQ_PCE_DEF_HPP
#define TPETRA_CRSMATRIX_UQ_PCE_DEF_HPP

// Specializations of some Tpetra::CrsMatrix methods for UQ::PCE

#include "Tpetra_CrsMatrix_def.hpp"

// These are macros; the code isn't compiled here, so we don't need to
// say here what namespace the code is in.  We _do_ need to put the
// code in the correct namespace when we use the macro.

#define TPETRA_CRSMATRIX_UQ_PCE_SPEC(Scalar,LocalOrdinal,GlobalOrdinal,Node) \
  template<>                                                            \
  Teuchos::RCP<MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node> > \
  CrsMatrix<Scalar , LocalOrdinal, GlobalOrdinal, Node>::        \
  getColumnMapMultiVector (const MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>& X_domainMap, \
                           const bool force) const                      \
  {                                                                     \
    using Teuchos::null;                                                \
    using Teuchos::RCP;                                                 \
    using Teuchos::rcp;                                                 \
                                                                        \
    TEUCHOS_TEST_FOR_EXCEPTION(                                         \
      ! this->hasColMap (), std::runtime_error, "Tpetra::CrsMatrix::getColumn" \
      "MapMultiVector: You may only call this method if the matrix has a " \
      "column Map.  If the matrix does not yet have a column Map, you should " \
      "first call fillComplete (with domain and range Map if necessary)."); \
                                                                        \
    TEUCHOS_TEST_FOR_EXCEPTION(                                         \
      ! this->getGraph ()->isFillComplete (), std::runtime_error, "Tpetra::" \
      "CrsMatrix::getColumnMapMultiVector: You may only call this method if " \
      "this matrix's graph is fill complete.");                         \
                                                                        \
    const size_t numVecs = X_domainMap.getNumVectors ();                \
    RCP<const import_type> importer = this->getGraph ()->getImporter (); \
    RCP<const map_type> colMap = this->getColMap ();                    \
                                                                        \
    RCP<MV> X_colMap;                                                   \
                                                                        \
    if (! importer.is_null () || force) {                               \
      if (importMV_.is_null () ||                                       \
          importMV_->getNumVectors () != numVecs ||                     \
          (Kokkos::dimension_scalar(*importMV_) !=                      \
           Kokkos::dimension_scalar(X_domainMap))) {                    \
        X_colMap = rcp (new MV (colMap, numVecs));                      \
        importMV_ = X_colMap;                                           \
      }                                                                 \
      else {                                                            \
        X_colMap = importMV_;                                           \
      }                                                                 \
    }                                                                   \
    return X_colMap;                                                    \
  }                                                                     \
                                                                        \
  template <>                                                           \
  Teuchos::RCP<MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node> > \
  CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::       \
  getRowMapMultiVector (const MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>& Y_rangeMap, \
                        const bool force) const                         \
  {                                                                     \
    using Teuchos::null;                                                \
    using Teuchos::RCP;                                                 \
    using Teuchos::rcp;                                                 \
                                                                        \
    TEUCHOS_TEST_FOR_EXCEPTION(                                         \
      ! this->getGraph ()->isFillComplete (), std::runtime_error, "Tpetra::" \
      "CrsMatrix::getRowMapMultiVector: You may only call this method if this " \
      "matrix's graph is fill complete.");                              \
                                                                        \
    const size_t numVecs = Y_rangeMap.getNumVectors ();                 \
    RCP<const export_type> exporter = this->getGraph ()->getExporter (); \
    RCP<const map_type> rowMap = this->getRowMap ();                    \
                                                                        \
    RCP<MV> Y_rowMap;                                                   \
                                                                        \
    if (! exporter.is_null () || force) {                               \
      if (exportMV_.is_null () ||                                       \
          exportMV_->getNumVectors () != numVecs ||                     \
          (Kokkos::dimension_scalar(*exportMV_) !=                      \
           Kokkos::dimension_scalar(Y_rangeMap))) {                     \
        Y_rowMap = rcp (new MV (rowMap, numVecs));                      \
        exportMV_ = Y_rowMap;                                           \
      }                                                                 \
      else {                                                            \
        Y_rowMap = exportMV_;                                           \
      }                                                                 \
    }                                                                   \
    return Y_rowMap;                                                    \
  }

#endif // STOKHOS_TPETRA_UQ_PCE_DEF_HPP
