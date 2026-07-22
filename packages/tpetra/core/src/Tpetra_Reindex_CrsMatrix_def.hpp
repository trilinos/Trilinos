// @HEADER
// *****************************************************************************
//          Tpetra: Templated Linear Algebra Services Package
//
// Copyright 2008 NTESS and the Tpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef TPETRA_REINDEX_CRSMATRIX_DEF_HPP
#define TPETRA_REINDEX_CRSMATRIX_DEF_HPP

#include <Tpetra_Reindex_CrsMatrix_decl.hpp>
#include <Tpetra_Vector_decl.hpp>
#include <Tpetra_Import_decl.hpp>

#include <vector>

/// \file Tpetra_Reindex_CrsMatrix_def.hpp
/// \brief Definition of the Tpetra::Reindex_CrsMatrix class
///
/// If you want to use Tpetra::Reindex_CrsMatrix, include
/// "Tpetra_Reindex_CrsMatrix.hpp", a file which CMake generates
/// and installs for you.
///
/// If you only want the declaration of Tpetra::Reindex_CrsMatrix,
/// include "Tpetra_Reindex_CrsMatrix_decl.hpp".

namespace Tpetra {

template <class Scalar,
          class LocalOrdinal,
          class GlobalOrdinal,
          class Node>
Reindex_CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    Reindex_CrsMatrix(Teuchos::RCP<Map<LocalOrdinal, GlobalOrdinal, Node> const> newRowMap)
  : ViewTransform<CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> >()
  , newRowMap_(newRowMap)
  , newColMap_(Teuchos::null) {
  // Nothing to do
}

template <class Scalar,
          class LocalOrdinal,
          class GlobalOrdinal,
          class Node>
Reindex_CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::~Reindex_CrsMatrix() {
  // Nothing to do
}

template <class Scalar,
          class LocalOrdinal,
          class GlobalOrdinal,
          class Node>
typename Reindex_CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::NewType
Reindex_CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::operator()(OriginalType const& origMatrix) {
  using cm_t = CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>;

  this->origObj_ = origMatrix;

  assert(origMatrix->getRowMap()->getLocalNumElements() == newRowMap_->getLocalNumElements());

  if ((origMatrix->getDomainMap()->getGlobalNumElements() == 0) &&
      (origMatrix->getRowMap()->getGlobalNumElements() == 0)) {
    // Construct a zero matrix as a placeholder, don't do reindexing analysis.
    this->newObj_ = Teuchos::rcp<cm_t>(new cm_t(origMatrix->getRowMap(), origMatrix->getColMap(), 0));
  } else {
    using v_t   = Vector<GlobalOrdinal, LocalOrdinal, GlobalOrdinal, Node>;
    using map_t = Map<LocalOrdinal, GlobalOrdinal, Node>;

    // Construct new column map
    v_t cols(origMatrix->getDomainMap());
    {
      size_t origDomainMap_localSize = origMatrix->getDomainMap()->getLocalNumElements();
      map_t tmpColMap(origMatrix->getDomainMap()->getGlobalNumElements(), origDomainMap_localSize, 0, origMatrix->getDomainMap()->getComm());
      Kokkos::deep_copy(Kokkos::subview(cols.getLocalViewDevice(Tpetra::Access::OverwriteAll), Kokkos::ALL(), 0),
                        tmpColMap.getMyGlobalIndicesDevice());
    }

    v_t newCols(origMatrix->getColMap());
    {
      using imp_t                        = Import<LocalOrdinal, GlobalOrdinal, Node>;
      Teuchos::RCP<const imp_t> importer = origMatrix->getCrsGraph()->getImporter();
      if (importer.is_null()) {
        newCols = cols;
      } else {
        newCols.doImport(cols, *importer, INSERT, false);
      }
    }

    using kv_t = Kokkos::View<GlobalOrdinal*, typename Node::device_type>;
    kv_t newColIndices;
    {
      auto newColsView = newCols.getLocalViewDevice(Tpetra::Access::ReadOnly);
      newColIndices    = kv_t("newColIndices", newColsView.extent(0));
      Kokkos::deep_copy(newColIndices, Kokkos::subview(newColsView, Kokkos::ALL(), 0));
    }

    this->newColMap_ = Teuchos::RCP<map_t>(new map_t(origMatrix->getColMap()->getGlobalNumElements(),
                                                     newColIndices,
                                                     origMatrix->getColMap()->getIndexBase(),
                                                     origMatrix->getColMap()->getComm()));

    // Create the new matrix
    typename cm_t::local_matrix_device_type tmpLocalMatrix("", origMatrix->getLocalMatrixDevice());
    Teuchos::RCP<cm_t> newMatrix = Teuchos::rcp<cm_t>(new cm_t(tmpLocalMatrix, this->newRowMap_, this->newColMap_));

    this->newObj_ = newMatrix;
  }

  return this->newObj_;
}

//
// Explicit instantiation macro
//
// Must be expanded from within the Tpetra namespace!
//

#define TPETRA_REINDEXCRSMATRIX_INSTANT(SCALAR, LO, GO, NODE) \
  template class Reindex_CrsMatrix<SCALAR, LO, GO, NODE>;

}  // namespace Tpetra

#endif  // TPETRA_REINDEX_CRSMATRIX_DEF_HPP
