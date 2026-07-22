// @HEADER
// *****************************************************************************
//          Tpetra: Templated Linear Algebra Services Package
//
// Copyright 2008 NTESS and the Tpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef TPETRA_SOLVERMAP_CRSMATRIX_DEF_HPP
#define TPETRA_SOLVERMAP_CRSMATRIX_DEF_HPP

/// \file Tpetra_SolverMap_CrsMatrix_def.hpp
/// \brief Definition of the Tpetra::SolverMap_CrsMatrix class
///
/// If you want to use Tpetra::SolverMap_CrsMatrix, include
/// "Tpetra_SolverMap_CrsMatrix.hpp", a file which CMake generates
/// and installs for you.
///
/// If you only want the declaration of Tpetra::SolverMap_CrsMatrix,
/// include "Tpetra_SolverMap_CrsMatrix_decl.hpp".

#include <Tpetra_SolverMap_CrsMatrix_decl.hpp>

#include <vector>

namespace Tpetra {

template <class Scalar,
          class LocalOrdinal,
          class GlobalOrdinal,
          class Node>
SolverMap_CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::SolverMap_CrsMatrix()
  : StructuralSameTypeTransform<CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> >()
  , newColMap_(Teuchos::null)
  , newGraph_(Teuchos::null) {
  // Nothing to do
}

template <class Scalar,
          class LocalOrdinal,
          class GlobalOrdinal,
          class Node>
SolverMap_CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::~SolverMap_CrsMatrix() {
  // Nothing to do
}

template <class Scalar,
          class LocalOrdinal,
          class GlobalOrdinal,
          class Node>
typename SolverMap_CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::NewType
SolverMap_CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::operator()(OriginalType const& origMatrix) {
  using map_t = Map<LocalOrdinal, GlobalOrdinal, Node>;

  assert(!origMatrix->isGloballyIndexed());

  this->origObj_ = origMatrix;

  // *******************************************************************
  // Step 1/7: Check if domain map and col map are different
  // *******************************************************************
  Teuchos::RCP<map_t const> origDomainMap = origMatrix->getDomainMap();
  Teuchos::RCP<map_t const> origColMap    = origMatrix->getColMap();

  typename map_t::local_map_type localOrigDomainMap(origDomainMap->getLocalMap());
  typename map_t::local_map_type localOrigColMap(origColMap->getLocalMap());

  if (origDomainMap->isLocallyFitted(*origColMap)) {
    this->newObj_ = this->origObj_;
  } else {
    using cg_t = CrsGraph<LocalOrdinal, GlobalOrdinal, Node>;
    using cm_t = CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>;

    // *****************************************************************
    // Step 2/7: Fill newColMap_globalColIndices with the global indices
    //           of all entries in origDomainMap
    // *****************************************************************

    // Initialize the value of 'newColMap_localSize'
    size_t const origDomainMap_localSize = origDomainMap->getLocalNumElements();
    size_t newColMap_localSize(origDomainMap_localSize);

    // Increase the value of 'newColMap_localSize' as necessary
    size_t newColMap_extraSize(0);
    size_t const origColMap_localSize = origColMap->getLocalNumElements();
    {
      Kokkos::parallel_reduce(
          "Tpetra::SolverMap_CrsMatrix::construct::newColMap_extraSize", Kokkos::RangePolicy<typename Node::device_type::execution_space, size_t>(0, origColMap_localSize), KOKKOS_LAMBDA(size_t const i, size_t& sizeToUpdate)->void {
            GlobalOrdinal const globalColIndex(localOrigColMap.getGlobalElement(i));
            if (localOrigDomainMap.getLocalElement(globalColIndex) == ::Tpetra::Details::OrdinalTraits<LocalOrdinal>::invalid()) {
              sizeToUpdate += 1;
            }
          },
          newColMap_extraSize);
    }
    newColMap_localSize += newColMap_extraSize;

    // Instantiate newColMap_globalColIndices with the correct size 'newColMap_localSize'
    Kokkos::View<GlobalOrdinal*, typename Node::device_type> newColMap_globalColIndices("", newColMap_localSize);

    // Fill newColMap_globalColIndices with the global indices of all entries in origDomainMap
    {
      Kokkos::parallel_for(
          "Tpetra::SolverMap_CrsMatrix::construct::copyDomainMapToNewColMap", Kokkos::RangePolicy<typename Node::device_type::execution_space, size_t>(0, origDomainMap_localSize), KOKKOS_LAMBDA(size_t const i)->void {
            newColMap_globalColIndices(i) = localOrigDomainMap.getGlobalElement(i);
          });
    }

    // *****************************************************************
    // Step 3/7: Append newColMap_globalColIndices with those origColMap
    //           entries that are not in newColMap_globalColIndices yet
    // *****************************************************************
    {
      Kokkos::parallel_scan(
          "Tpetra::SolverMap_CrsMatrix::construct::appendNewColMap", Kokkos::RangePolicy<typename Node::device_type::execution_space, size_t>(0, origColMap_localSize), KOKKOS_LAMBDA(size_t const i, size_t& jToUpdate, bool const final)->void {
            GlobalOrdinal const globalColIndex(localOrigColMap.getGlobalElement(i));
            if (localOrigDomainMap.getLocalElement(globalColIndex) == ::Tpetra::Details::OrdinalTraits<LocalOrdinal>::invalid()) {
              if (final) {
                newColMap_globalColIndices(origDomainMap_localSize + jToUpdate) = globalColIndex;
              }
              jToUpdate += 1;
            }
          });
    }

    // *****************************************************************
    // Step 4/7: Create a new column map using newColMap_globalColIndices
    // *****************************************************************
    Teuchos::RCP<map_t const> origRowMap        = origMatrix->getRowMap();
    Teuchos::RCP<Teuchos::Comm<int> const> Comm = origRowMap->getComm();
    size_t const newColMap_localNumCols         = newColMap_globalColIndices.size();
    size_t newColMap_globalNumCols(0);
    Teuchos::reduceAll(*Comm, Teuchos::REDUCE_SUM, 1, &newColMap_localNumCols, &newColMap_globalNumCols);

    newColMap_ = Teuchos::rcp<map_t>(new map_t(newColMap_globalNumCols, newColMap_globalColIndices, origDomainMap->getIndexBase(), Comm));

    // *****************************************************************
    // Step 5/7: Create new graph
    // *****************************************************************
    cg_t const* origGraph = dynamic_cast<cg_t const*>(origMatrix->getGraph().get());
    newGraph_             = Teuchos::rcp<cg_t>(new cg_t(*origGraph));
    newGraph_->resumeFill();
    newGraph_->reindexColumns(newColMap_);
    newGraph_->fillComplete();

    // *****************************************************************
    // Step 6/7: Create new CRS matrix
    // *****************************************************************
    Teuchos::RCP<cm_t> newMatrix = Teuchos::rcp<cm_t>(new cm_t(*origMatrix, newGraph_));

    // *****************************************************************
    // Step 7/7: Update newObj_
    // *****************************************************************
    this->newObj_ = newMatrix;
  }

  return this->newObj_;
}

//
// Explicit instantiation macro
//
// Must be expanded from within the Tpetra namespace!
//

#define TPETRA_SOLVERMAPCRSMATRIX_INSTANT(SCALAR, LO, GO, NODE) \
  template class SolverMap_CrsMatrix<SCALAR, LO, GO, NODE>;

}  // namespace Tpetra

#endif  // TPETRA_SOLVERMAP_CRSMATRIX_DEF_HPP
