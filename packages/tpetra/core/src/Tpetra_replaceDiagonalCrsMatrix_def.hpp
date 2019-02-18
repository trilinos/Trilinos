// @HEADER
// ***********************************************************************
//
//          Tpetra: Templated Linear Algebra Services Package
//                 Copyright (2008) Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
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
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
//
// ************************************************************************
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
     "Input matrix A must have a nonnull column Map.");

  const map_type& rowMap = *rowMapPtr;
  const map_type& colMap = *colMapPtr;
  const LO myNumRows = static_cast<LO>(matrix.getNodeNumRows());
  const bool isFillCompleteOnInput = matrix.isFillComplete();

  if (Tpetra::Details::Behavior::debug()) {
    TEUCHOS_TEST_FOR_EXCEPTION(!rowMap.isSameAs(*newDiag.getMap()), std::runtime_error,
        "Row map of matrix and map of input vector do not match.");
  }

  crs_matrix_type::execution_space::fence(); // for UVM's sake

  if (isFillCompleteOnInput)
    matrix.resumeFill();

  Teuchos::ArrayRCP<const SC> newDiagData = newDiag.getVector(0)->getData();
  LO numReplacedEntriesPerRow = 0;

  // Loop over all local rows to replace the diagonal entry row by row
  for (LO lclRowInd = 0; lclRowInd < myNumRows; ++lclRowInd) {

    // Get row and column indices of the diagonal entry
    const GO gblInd = rowMap.getGlobalElement(lclRowInd);
    const LO lclColInd = colMap.getLocalElement(gblInd);

    const SC vals[] = {static_cast<SC>(newDiagData[lclRowInd])};
    const LO cols[] = {lclColInd};

    // Do the actual replacement of the diagonal element
    numReplacedEntriesPerRow = matrix.replaceLocalValues(lclRowInd, oneLO, vals, cols);

    // Check for success of replacement
    if (numReplacedEntriesPerRow == oneLO) {
      ++numReplacedDiagEntries;
    } else {
      TEUCHOS_TEST_FOR_EXCEPTION(true, std::runtime_error,
        "Number of replaced entries in this row is not equal to one.  "
        "It has to be exactly one, since we want to replace the diagonal element.");
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
