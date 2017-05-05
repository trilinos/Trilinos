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

#ifndef TPETRA_DETAILS_UNPACKCRSMATRIX_HPP
#define TPETRA_DETAILS_UNPACKCRSMATRIX_HPP

#include "TpetraCore_config.h"
#include "Teuchos_Array.hpp"
#include "Teuchos_ArrayView.hpp"
#include "Tpetra_Details_OrdinalTraits.hpp"
#include <memory>
#include <string>

/// \file Tpetra_Details_unpackCrsMatrix.hpp
/// \brief Functions for unpacking the entries of a Tpetra::CrsMatrix
///   for communication, in the case where it is valid to go to the
///   KokkosSparse::CrsMatrix (local sparse matrix data structure)
///   directly.
/// \warning This file, and its contents, are implementation details
///   of Tpetra.  The file itself or its contents may disappear or
///   change at any time.

namespace Tpetra {

#ifndef DOXYGEN_SHOULD_SKIP_THIS
// Forward declaration of Distributor
class Distributor;
#endif // DOXYGEN_SHOULD_SKIP_THIS

//
// Users must never rely on anything in the Details namespace.
//
namespace Details {

template<class LocalMatrix, class LocalMap>
typename LocalMatrix::ordinal_type
combineCrsMatrixValues(LocalMatrix& lclMatrix,
                       LocalMap& lclColMap,
                       const typename LocalMatrix::ordinal_type lclRow,
                       const typename LocalMatrix::ordinal_type numEnt,
                       const typename LocalMatrix::value_type vals[],
                       const typename LocalMap::global_ordinal_type cols[],
                       const Tpetra::CombineMode combineMode,
                       const bool atomic)
{
  typedef typename LocalMatrix::ordinal_type LO;

  // INSERT doesn't make sense for a static graph, since you
  // aren't allowed to change the structure of the graph.
  // However, all the other combine modes work.

  LO numValid = 0; // number of valid input column indices
  if (combineMode == ADD) {
    for (int k=0; k<numEnt; k++) {
      LO lclColInd = lclColMap.getLocalElement(cols[k]);
      numValid += lclMatrix.sumIntoValues(lclRow, &lclColInd, 1, &vals[k], false, atomic);
    }
  }
  else if (combineMode == REPLACE) {
    for (int k=0; k<numEnt; k++) {
      LO lclColInd = lclColMap.getLocalElement(cols[k]);
      numValid += lclMatrix.replaceValues(lclRow, &lclColInd, 1, &vals[k], false, atomic);
    }
  }
  else if (combineMode == ABSMAX) {
    TEUCHOS_TEST_FOR_EXCEPTION(
        true, std::invalid_argument,
        "ABSMAX combine mode is not yet implemented for a matrix that has a "
        "static graph (i.e., was constructed with the CrsMatrix constructor "
        "that takes a const CrsGraph pointer).");
  }
  else if (combineMode == INSERT) {
    TEUCHOS_TEST_FOR_EXCEPTION(
        true, std::invalid_argument,
        "INSERT combine mode is not allowed if the matrix has a static graph "
        "(i.e., was constructed with the CrsMatrix constructor that takes a "
        "const CrsGraph pointer).");
  }
  else {
    TEUCHOS_TEST_FOR_EXCEPTION(
        true, std::logic_error, "Invalid combine mode; should never get "
        "here!  Please report this bug to the Tpetra developers.");
  }

  return numValid;

}

template<class LocalMatrix, class LocalMap>
bool
unpackCrsMatrixRow(LocalMatrix& lclMatrix,
                   LocalMap& lclColMap,
                   std::unique_ptr<std::string>& errStr,
                   typename LocalMatrix::value_type* const valInTmp,
                   typename LocalMap::global_ordinal_type* const indInTmp,
                   const size_t tmpSize,
                   const char* const valIn,
                   const char* const indIn,
                   const size_t numEnt,
                   const typename LocalMatrix::ordinal_type lclRow,
                   const Tpetra::CombineMode combineMode,
                   const bool atomic)
{
  if (tmpSize < numEnt || (numEnt != 0 && (valInTmp == NULL || indInTmp == NULL))) {
    return false;
  }

  typedef typename LocalMatrix::value_type LO;
  typedef typename LocalMap::global_ordinal_type GO;
  memcpy(valInTmp, valIn, numEnt * sizeof(LO));
  memcpy(indInTmp, indIn, numEnt * sizeof(GO));

  // FIXME (mfh 23 Mar 2017) It would make sense to use the return
  // value here as more than just a "did it succeed" Boolean test.

  // FIXME (mfh 23 Mar 2017) CrsMatrix_NonlocalSumInto_Ignore test
  // expects this method to ignore incoming entries that do not
  // exist on the process that owns those rows.  We would like to
  // distinguish between "errors" resulting from ignored entries,
  // vs. actual errors.

  //const LO numModified =
  combineCrsMatrixValues<LocalMatrix,LocalMap>(
      lclMatrix, lclColMap, lclRow, numEnt, valInTmp, indInTmp, combineMode, atomic);
  return true; // FIXME (mfh 23 Mar 2013) See above.
  //return numModified == numEnt;
}

/// \brief Unpack the imported column indices and values, and combine into matrix.
///
/// \warning The allowed \c combineMode are:
///   ADD, REPLACE, and ABSMAX. INSERT is not allowed.
template<class LocalMatrix, class LocalMap>
bool
unpackCrsMatrixAndCombine(
    LocalMatrix& lclMatrix,
    const LocalMap& lclColMap,
    std::unique_ptr<std::string>& errStr,
    const Teuchos::ArrayView<const typename LocalMatrix::ordinal_type>& importLIDs,
    const Teuchos::ArrayView<const char>& imports,
    const Teuchos::ArrayView<const size_t>& numPacketsPerLID,
    size_t constantNumPackets,
    Distributor & /* distor */,
    CombineMode combineMode,
    const bool atomic)
{
  //typedef LocalMatrix::ordinal_type LO;
  typedef typename LocalMatrix::value_type IST;
  typedef typename LocalMatrix::ordinal_type LO;
  typedef typename LocalMap::global_ordinal_type GO;
  typedef typename Teuchos::ArrayView<const LO>::size_type size_type;

  const size_type numImportLIDs = importLIDs.size();
  if (numImportLIDs != numPacketsPerLID.size()) {
    std::ostringstream os;
    os << "importLIDs.size() (" << numImportLIDs << ") != "
       << "numPacketsPerLID.size() (" << numPacketsPerLID.size() << ").";
    if (errStr.get() == NULL) {
      errStr = std::unique_ptr<std::string>(new std::string());
    }
    *errStr = os.str();
    return false;
  }

  // If a sanity check fails, keep track of some state at the
  // "first" place where it fails.  After the first failure, "run
  // through the motions" until the end of this method, then raise
  // an error with an informative message.
  size_type firstBadIndex = 0;
  size_t firstBadOffset = 0;
  size_t firstBadExpectedNumBytes = 0;
  size_t firstBadNumBytes = 0;
  LO firstBadNumEnt = 0;
  // We have sanity checks for three kinds of errors:
  //
  //   1. Offset into array of all the incoming data (for all rows)
  //      is out of bounds
  //   2. Too few bytes of incoming data for a row, given the
  //      reported number of entries in those incoming data
  //   3. Error in unpacking the row's incoming data
  //
  bool outOfBounds = false;
  bool wrongNumBytes = false;
  bool unpackErr = false;

  const size_t bufSize = static_cast<size_t>(imports.size());
  const char* const importsRawPtr = imports.getRawPtr();
  size_t offset = 0;

  // Temporary storage for incoming values and indices.  We need
  // this because the receive buffer does not align storage; it's
  // just contiguous bytes.  In order to avoid violating ANSI
  // aliasing rules, we memcpy each incoming row's data into these
  // temporary arrays.  We double their size every time we run out
  // of storage.
  std::vector<IST> valInTmp;
  std::vector<GO> indInTmp;
  for (size_type i = 0; i < numImportLIDs; ++i) {
    const LO lclRow = importLIDs[i];
    const size_t numBytes = numPacketsPerLID[i];

    if (numBytes > 0) { // there is actually something in the row
      const char* const numEntBeg = importsRawPtr + offset;
      const char* const numEntEnd = numEntBeg + sizeof(LO);

      // Now we know how many entries to expect in the received data
      // for this row.
      LO numEnt = 0;
      memcpy(&numEnt, numEntBeg, sizeof(LO));

      const char* const valBeg = numEntEnd;
      const char* const valEnd =
        valBeg + static_cast<size_t>(numEnt) * sizeof(IST);
      const char* const indBeg = valEnd;
      const size_t expectedNumBytes = sizeof(LO) +
        static_cast<size_t>(numEnt) * (sizeof(IST) + sizeof(GO));

      if (expectedNumBytes > numBytes) {
        firstBadIndex = i;
        firstBadOffset = offset;
        firstBadExpectedNumBytes = expectedNumBytes;
        firstBadNumBytes = numBytes;
        firstBadNumEnt = numEnt;
        wrongNumBytes = true;
        break;
      }
      if (offset > bufSize || offset + numBytes > bufSize) {
        firstBadIndex = i;
        firstBadOffset = offset;
        firstBadExpectedNumBytes = expectedNumBytes;
        firstBadNumBytes = numBytes;
        firstBadNumEnt = numEnt;
        outOfBounds = true;
        break;
      }
      size_t tmpNumEnt = static_cast<size_t>(valInTmp.size());
      if (tmpNumEnt < static_cast<size_t>(numEnt) ||
          static_cast<size_t>(indInTmp.size()) < static_cast<size_t>(numEnt)) {
        // Double the size of the temporary arrays for incoming data.
        tmpNumEnt = std::max(static_cast<size_t>(numEnt), tmpNumEnt * 2);
        valInTmp.resize(tmpNumEnt);
        indInTmp.resize(tmpNumEnt);
      }

      unpackErr = ! unpackCrsMatrixRow(lclMatrix, lclColMap, errStr,
          valInTmp.data(), indInTmp.data(), tmpNumEnt,
          valBeg, indBeg, numEnt, lclRow, combineMode, atomic);

      if (unpackErr) {
        firstBadIndex = i;
        firstBadOffset = offset;
        firstBadExpectedNumBytes = expectedNumBytes;
        firstBadNumBytes = numBytes;
        firstBadNumEnt = numEnt;
        break;
      }
      offset += numBytes;
    }
  }

  if (wrongNumBytes || outOfBounds || unpackErr) {
    std::ostringstream os;
    if (wrongNumBytes) {
      os << "At index i = " << firstBadIndex
         << ", expectedNumBytes > numBytes.";
    }
    else if (outOfBounds) {
      os << "First invalid offset into 'imports' "
         << "unpack buffer at index i = " << firstBadIndex << ".";
    }
    else if (unpackErr) {
      os << "First error in unpackRow() at index i = "
         << firstBadIndex << ".";
    }
    os << "  importLIDs[i]: " << importLIDs[firstBadIndex]
       << ", bufSize: " << bufSize
       << ", offset: " << firstBadOffset
       << ", numBytes: " << firstBadNumBytes
       << ", expectedNumBytes: " << firstBadExpectedNumBytes
       << ", numEnt: " << firstBadNumEnt;
    if (errStr.get() == NULL) {
      errStr = std::unique_ptr<std::string>(new std::string());
    }
    *errStr = os.str();
    return false;
  }
  return true;
}

} // namespace Details
} // namespace Tpetra

#endif // TPETRA_DETAILS_UNPACKCRSMATRIX_HPP
