// @HEADER
// *****************************************************************************
//          Tpetra: Templated Linear Algebra Services Package
//
// Copyright 2008 NTESS and the Tpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef TPETRA_REPLACEDIAGONALCRSMATRIX_DEF_HPP
#define TPETRA_REPLACEDIAGONALCRSMATRIX_DEF_HPP

/// \file Tpetra_replaceDiagonalCrsMatrix_def.hpp
/// \brief Definition of Tpetra::repalceDiagonalCrsMatrix

#include "Tpetra_Details_Behavior.hpp"
#include "Tpetra_CrsMatrix.hpp"
#include "Tpetra_Vector.hpp"

namespace Tpetra {

template<class SC, class LO, class GO, class NT>
LO
replaceDiagonalCrsMatrix (CrsMatrix<SC, LO, GO, NT>& matrix,
    const Vector<SC, LO, GO, NT>& newDiag) {

  using map_type = Map<LO, GO, NT>;
  using crs_matrix_type = CrsMatrix<SC, LO, GO, NT>;
  // using vec_type = ::Tpetra::Vector<SC, LO, GO, NT>;

  const LO oneLO = Teuchos::OrdinalTraits<LO>::one();

  // Count number of successfully replaced diagonal entries
  LO numReplacedDiagEntries = 0;

  // Extract some useful data
  auto rowMapPtr = matrix.getRowMap();
  if (rowMapPtr.is_null() || rowMapPtr->getComm().is_null()) {
    // Processes on which the row Map or its communicator is null
    // don't participate.  Users shouldn't even call this method on
    // those processes.
    return numReplacedDiagEntries;
  }
  auto colMapPtr = matrix.getColMap();

  TEUCHOS_TEST_FOR_EXCEPTION
    (colMapPtr.get () == nullptr, std::invalid_argument,
     "replaceDiagonalCrsMatrix: "
     "Input matrix must have a nonnull column Map.");

  const map_type& rowMap = *rowMapPtr;
  const map_type& colMap = *colMapPtr;
  const LO myNumRows = static_cast<LO>(matrix.getLocalNumRows());
  const bool isFillCompleteOnInput = matrix.isFillComplete();

  if (Tpetra::Details::Behavior::debug()) {
    TEUCHOS_TEST_FOR_EXCEPTION(!rowMap.isSameAs(*newDiag.getMap()), std::runtime_error,
        "Row map of matrix and map of input vector do not match.");
  }

  // KJ: This fence is necessary for UVM. Views used in the row map and colmap
  // can use UVM and they are accessed in the following routine. So, we need to
  // make sure that the values are available for touching in host.
  typename crs_matrix_type::execution_space().fence();

  if (isFillCompleteOnInput)
    matrix.resumeFill();

  Teuchos::ArrayRCP<const SC> newDiagData = newDiag.getVector(0)->getData();
  LO numReplacedEntriesPerRow = 0;

  auto invalid = Teuchos::OrdinalTraits<LO>::invalid();

  // Loop over all local rows to replace the diagonal entry row by row
  for (LO lclRowInd = 0; lclRowInd < myNumRows; ++lclRowInd) {

    // Get row and column indices of the diagonal entry
    const GO gblInd = rowMap.getGlobalElement(lclRowInd);
    const LO lclColInd = colMap.getLocalElement(gblInd);

    // If the row map is not one-to-one, the diagonal may not be on this proc.
    // Skip this row; some processor will have the diagonal for this row.
    if (lclColInd == invalid) continue;

    const SC vals[] = {static_cast<SC>(newDiagData[lclRowInd])};
    const LO cols[] = {lclColInd};

    // Do the actual replacement of the diagonal element, if on this proc
    numReplacedEntriesPerRow = matrix.replaceLocalValues(lclRowInd, oneLO,
                                                         vals, cols);

    // Check for success of replacement. 
    // numReplacedEntriesPerRow is one if the diagonal was replaced.
    // numReplacedEntriesPerRow is zero if the diagonal is not on 
    // this processor.  For example, in a 2D matrix distribution, gblInd may
    // be in both the row and column map, but the diagonal may not be on 
    // this processor.
    if (numReplacedEntriesPerRow == oneLO) {
      ++numReplacedDiagEntries;
    }
  }

  if (isFillCompleteOnInput)
    matrix.fillComplete();

  return numReplacedDiagEntries;
}

} // namespace Tpetra
//
// Explicit instantiation macro
//
// Must be expanded from within the Tpetra namespace!
//

#define TPETRA_REPLACEDIAGONALCRSMATRIX_INSTANT(SCALAR,LO,GO,NODE)      \
                                                                        \
  template                                                              \
  LO replaceDiagonalCrsMatrix(                                          \
                              CrsMatrix<SCALAR, LO, GO, NODE>& matrix, \
                              const Vector<SCALAR, LO, GO, NODE>& newDiag); \
                                                                        \

#endif // #ifndef TPETRA_REPLACEDIAGONALCRSMATRIX_DEF_HPP
