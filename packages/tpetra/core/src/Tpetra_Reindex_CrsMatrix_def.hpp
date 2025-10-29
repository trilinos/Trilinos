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
    using map_t = Map<LocalOrdinal, GlobalOrdinal, Node>;
    using imp_t = Import<LocalOrdinal, GlobalOrdinal, Node>;
    using v_t   = Vector<GlobalOrdinal, LocalOrdinal, GlobalOrdinal, Node>;

    // Construct new column map
    v_t cols(origMatrix->getDomainMap());
    {
      size_t origDomainMap_localSize = origMatrix->getDomainMap()->getLocalNumElements();
      map_t tmpColMap(origMatrix->getDomainMap()->getGlobalNumElements(), origDomainMap_localSize, 0, origMatrix->getDomainMap()->getComm());
      for (size_t i(0); i < origDomainMap_localSize; ++i) {
        cols.replaceLocalValue(i, tmpColMap.getGlobalElement(i));
      }
    }

    imp_t importer(origMatrix->getDomainMap(), origMatrix->getColMap());
    v_t newCols(origMatrix->getColMap());
    newCols.doImport(cols, importer, INSERT, false);

    auto newColsView = newCols.getLocalViewDevice(Tpetra::Access::ReadOnly);
    size_t newColsSize(newColsView.extent(0));
    Kokkos::View<GlobalOrdinal*, typename Node::device_type> newColIndices("newColIndices", newColsSize);
    {
      using exec_space = typename Node::device_type::execution_space;
      Kokkos::parallel_for(
          "Tpetra::Reindex_CrsMatrix::operator()",
          Kokkos::RangePolicy<exec_space, size_t>(0, newColsSize),
          KOKKOS_LAMBDA(size_t const i)->void {
	    newColIndices(i) = newColsView(i,0);
          });
    }

    this->newColMap_ = Teuchos::RCP<map_t>(new map_t(origMatrix->getColMap()->getGlobalNumElements(),
                                                     newColIndices,
                                                     origMatrix->getColMap()->getIndexBase(),
                                                     origMatrix->getColMap()->getComm()));

    // Create the new matrix
    auto origMatrixLocal = origMatrix->getLocalMatrixDevice();
    Teuchos::RCP<cm_t> newMatrix = Teuchos::rcp<cm_t>(new cm_t(origMatrixLocal, this->newRowMap_, this->newColMap_));

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
