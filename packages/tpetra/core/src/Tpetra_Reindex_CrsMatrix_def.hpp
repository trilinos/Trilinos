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
    newCols.doImport(cols  // const SrcDistObject& source
                     ,
                     importer  // const Import<LocalOrdinal, GlobalOrdinal, Node>& importer
                     ,
                     INSERT  // const CombineMode CM
                     ,
                     false  // const bool restrictedMode
    );

    Teuchos::ArrayRCP<GlobalOrdinal const> newColIndicesArray = newCols.getData();
    std::vector<GlobalOrdinal> newColIndicesVector(newColIndicesArray.size());
    for (size_t j(0); j < newColIndicesVector.size(); ++j) {
      newColIndicesVector[j] = newColIndicesArray[j];
    }
    this->newColMap_ = Teuchos::RCP<map_t>(new map_t(origMatrix->getColMap()->getGlobalNumElements()  // const global_size_t       numGlobalElements
                                                     ,
                                                     newColIndicesVector.data()  // const global_ordinal_type indexList[]
                                                     ,
                                                     origMatrix->getColMap()->getLocalNumElements()  // const local_ordinal_type  indexListSize
                                                     ,
                                                     origMatrix->getColMap()->getIndexBase()  // const global_ordinal_type indexBase
                                                     ,
                                                     origMatrix->getColMap()->getComm()  // const Teuchos::RCP<const Teuchos::Comm<int> >& comm
                                                     ));

    // Create the new matrix
    size_t const origMatrix_maxNumEntries = origMatrix->getGlobalMaxNumRowEntries();
    Teuchos::RCP<cm_t> newMatrix          = Teuchos::rcp<cm_t>(new cm_t(this->newRowMap_, this->newColMap_, origMatrix_maxNumEntries));

    std::vector<Scalar> newMatrix_localValues(origMatrix_maxNumEntries);
    std::vector<LocalOrdinal> newMatrix_localIndices(origMatrix_maxNumEntries);

    typename cm_t::local_inds_host_view_type origMatrix_localIndices;
    typename cm_t::values_host_view_type origMatrix_localValues;

    size_t const newMatrix_localNumRows = newMatrix->getLocalNumRows();
    for (size_t i(0); i < newMatrix_localNumRows; ++i) {
      origMatrix->getLocalRowView(i, origMatrix_localIndices, origMatrix_localValues);

      size_t const numEntries(origMatrix_localIndices.size());
      for (size_t j(0); j < numEntries; ++j) {
        newMatrix_localValues[j]  = origMatrix_localValues[j];
        newMatrix_localIndices[j] = origMatrix_localIndices[j];
      }

      newMatrix->insertLocalValues(i  // const LocalOrdinal localRow
                                   ,
                                   numEntries  // const LocalOrdinal numEnt
                                   ,
                                   newMatrix_localValues.data()  // const Scalar       vals[]
                                   ,
                                   newMatrix_localIndices.data()  // const LocalOrdinal cols[]
                                   ,
                                   INSERT  // const CombineMode  CM
      );
    }

    newMatrix->fillComplete();

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
