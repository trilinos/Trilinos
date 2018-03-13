// @HEADER
// ***********************************************************************
//
//                           Stokhos Package
//                 Copyright (2009) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Eric T. Phipps (etphipp@sandia.gov).
//
// ***********************************************************************
// @HEADER

#ifndef TPETRA_CRSMATRIX_UQ_PCE_DEF_HPP
#define TPETRA_CRSMATRIX_UQ_PCE_DEF_HPP

// Specializations of some Tpetra::CrsMatrix methods for UQ::PCE

#include "Tpetra_CrsMatrix_def.hpp"

namespace Tpetra {

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

} // namespace Tpetra

#endif // STOKHOS_TPETRA_UQ_PCE_DEF_HPP
