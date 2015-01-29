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

#ifndef TPETRA_IMPORT_UTIL2_HPP
#define TPETRA_IMPORT_UTIL2_HPP

///
/// \file Tpetra_Import_Util2.hpp
/// \brief Utility functions for packing and unpacking sparse matrix entries.
///

#include "Tpetra_ConfigDefs.hpp" // for map, vector, string, and iostream
#include "Tpetra_Import.hpp"
#include "Tpetra_HashTable.hpp"
#include "Tpetra_Map.hpp"
#include "Tpetra_Util.hpp"
#include "Tpetra_Distributor.hpp"
#include <Teuchos_Array.hpp>
#include <utility>

// Tpetra::CrsMatrix uses the functions below in its implementation.
// To avoid a circular include issue, only include the declarations
// for CrsMatrix.  We will include the definition after the functions
// here have been defined.
#include "Tpetra_CrsMatrix_decl.hpp"

namespace Tpetra {
namespace Import_Util {

/// \brief Special version of Tpetra::CrsMatrix::packAndPrepare that
///   also packs owning process ranks.
///
/// One of several short-cut routines to optimize fill complete for
/// the special case of sparse matrix-matrix multiply.
///
/// Note: The SourcePids vector should contain a list of owning PIDs
/// for each column in the ColMap, as from Epetra_Util::GetPids,
/// without the "-1 for local" option being used.  This routine is
/// basically Tpetra::CrsMatrix::packAndPrepare, but it packs the
/// owning PIDs as well as the GIDs.
///
/// \warning This method is intended for expert developer use only,
///   and should never be called by user code.
template<typename Scalar,
         typename LocalOrdinal,
         typename GlobalOrdinal,
         typename Node>
void
packAndPrepareWithOwningPIDs (const CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>& SourceMatrix,
                              const Teuchos::ArrayView<const LocalOrdinal>& exportLIDs,
                              Teuchos::Array<char>& exports,
                              const Teuchos::ArrayView<size_t>& numPacketsPerLID,
                              size_t& constantNumPackets,
                              Distributor &distor,
                              const Teuchos::ArrayView<const int>& SourcePids);

/// \brief Special version of Tpetra::CrsMatrix::unpackAndCombine
///   that also unpacks owning process ranks.
///
/// It belongs with packAndPrepareWithOwningPIDs() (see above).
///
/// Perform the count for unpacking the imported column indices pids,
/// and values, and combining them into matrix.  Return (a ceiling on)
/// the number of local stored entries ("nonzeros") in the matrix.  If
/// there are no shared rows in the SourceMatrix this count is exact.
///
/// Note: This routine also counts the copyAndPermute nonzeros in
/// addition to those that come in via import.
///
/// \warning This method is intended for expert developer use
///   only, and should never be called by user code.
template<typename Scalar,
         typename LocalOrdinal,
         typename GlobalOrdinal,
         typename Node>
size_t
unpackAndCombineWithOwningPIDsCount (const CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>& SourceMatrix,
                                     const Teuchos::ArrayView<const LocalOrdinal> &importLIDs,
                                     const Teuchos::ArrayView<const char> &imports,
                                     const Teuchos::ArrayView<size_t> &numPacketsPerLID,
                                     size_t constantNumPackets,
                                     Distributor &distor,
                                     CombineMode combineMode,
                                     size_t numSameIDs,
                                     const Teuchos::ArrayView<const LocalOrdinal>& permuteToLIDs,
                                     const Teuchos::ArrayView<const LocalOrdinal>& permuteFromLIDs);

/// \brief unpackAndCombineIntoCrsArrays
///
/// \note You should call unpackAndCombineWithOwningPIDsCount first
///   and allocate all arrays accordingly, before calling this
///   function.
///
/// Note: The SourcePids vector (on input) should contain owning PIDs
/// for each column in the (source) ColMap, as from
/// Tpetra::Import_Util::getPids, with the "-1 for local" option being
/// used.
///
/// Note: The TargetPids vector (on output) will contain owning PIDs
/// for each entry in the matrix, with the "-1 for local" for locally
/// owned entries.
template<typename Scalar,
         typename LocalOrdinal,
         typename GlobalOrdinal,
         typename Node>
void
unpackAndCombineIntoCrsArrays (const CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>& SourceMatrix,
                               const Teuchos::ArrayView<const LocalOrdinal>& importLIDs,
                               const Teuchos::ArrayView<const char>& imports,
                               const Teuchos::ArrayView<size_t>& numPacketsPerLID,
                               size_t constantNumPackets,
                               Distributor &distor,
                               CombineMode combineMode,
                               size_t numSameIDs,
                               const Teuchos::ArrayView<const LocalOrdinal>& permuteToLIDs,
                               const Teuchos::ArrayView<const LocalOrdinal>& permuteFromLIDs,
                               size_t TargetNumRows,
                               size_t TargetNumNonzeros,
                               int MyTargetPID,
                               const Teuchos::ArrayView<size_t>& rowPointers,
                               const Teuchos::ArrayView<GlobalOrdinal>& columnIndices,
                               const Teuchos::ArrayView<Scalar>& values,
                               const Teuchos::ArrayView<const int>& SourcePids,
                               Teuchos::Array<int>& TargetPids);

/// \brief Sort the entries of the (raw CSR) matrix by column index
///   within each row.
template<typename Scalar, typename Ordinal>
void
sortCrsEntries (const Teuchos::ArrayView<size_t>& CRS_rowptr,
                const Teuchos::ArrayView<Ordinal>& CRS_colind,
                const Teuchos::ArrayView<Scalar>&CRS_vals);

/// \brief Sort and merge the entries of the (raw CSR) matrix by
///   column index within each row.
///
/// Entries with the same column index get merged additively.
template<typename Scalar, typename Ordinal>
void
sortAndMergeCrsEntries (const Teuchos::ArrayView<size_t>& CRS_rowptr,
                        const Teuchos::ArrayView<Ordinal>& CRS_colind,
                        const Teuchos::ArrayView<Scalar>& CRS_vals);

/// \brief lowCommunicationMakeColMapAndReindex
///
/// If you know the owning PIDs already, you can make the colmap a lot
/// less expensively.  If LocalOrdinal and GlobalOrdinal are the same,
/// you can (and should) use the same array for both columnIndices_LID
/// and columnIndices_GID.  This routine works just fine "in place."
///
/// Note: The owningPids vector (on input) should contain owning PIDs
/// for each entry in the matrix, like that generated by
/// Tpetra::Import_Util::unpackAndCombineIntoCrsArrays routine.  Note:
/// This method will return a Teuchos::Array of the remotePIDs, used for
/// construction of the importer.
///
/// \warning This method is intended for expert developer use only,
///   and should never be called by user code.
template <typename LocalOrdinal, typename GlobalOrdinal, typename Node>
void
lowCommunicationMakeColMapAndReindex (const Teuchos::ArrayView<const size_t> &rowPointers,
                                      const Teuchos::ArrayView<LocalOrdinal> &columnIndices_LID,
                                      const Teuchos::ArrayView<GlobalOrdinal> &columnIndices_GID,
                                      const Tpetra::RCP<const Tpetra::Map<LocalOrdinal,GlobalOrdinal,Node> > & domainMap,
                                      const Teuchos::ArrayView<const int> &owningPids,
                                      Teuchos::Array<int> &remotePids,
                                      Teuchos::RCP<const Tpetra::Map<LocalOrdinal,GlobalOrdinal,Node> > & colMap);

} // namespace Import_Util
} // namespace Tpetra


//
// Implementations
//

namespace { // (anonymous)

  template<class T>
  size_t
  packArray (const Teuchos::ArrayView<char>& outBuf,
             const Teuchos::ArrayView<const T>& inBuf,
             const size_t numEnt)
  {
    const size_t numBytes = numEnt * sizeof (T);
    memcpy (outBuf.getRawPtr (), inBuf.getRawPtr (), numBytes);
    return numBytes;
  }

  template<class T>
  size_t
  unpackArray (const Teuchos::ArrayView<T>& outBuf,
               const Teuchos::ArrayView<const char>& inBuf,
               const size_t numEnt)
  {
    const size_t numBytes = numEnt * sizeof (T);
    memcpy (outBuf.getRawPtr (), inBuf.getRawPtr (), numBytes);
    return numBytes;
  }

  template<class T>
  size_t
  packValue (const Teuchos::ArrayView<char>& outBuf,
             const T& inVal)
  {
    const size_t numBytes = sizeof (T);
    memcpy (outBuf.getRawPtr (), &inVal, numBytes);
    return numBytes;
  }

  template<class T>
  size_t
  unpackValue (T& outVal,
               const Teuchos::ArrayView<const char>& inBuf)
  {
    const size_t numBytes = sizeof (T);
    memcpy (&outVal, inBuf.getRawPtr (), numBytes);
    return numBytes;
  }

  // Return the number of bytes required to pack that row's entries.
  template<class ST, class LO, class GO>
  size_t
  packRowCount (const size_t numEnt)
  {
    if (numEnt == 0) {
      // Empty rows always take zero bytes, to ensure sparsity.
      return 0;
    }
    else {
      // Store the number of entries as a local index (LO).
      const size_t numEntLen = sizeof (LO);
      const size_t gidsLen = numEnt * sizeof (GO);
      const size_t pidsLen = numEnt * sizeof (int);
      const size_t valsLen = numEnt * sizeof (ST);

      return numEntLen + gidsLen + pidsLen + valsLen;
    }
  }

  template<class LO>
  size_t
  unpackRowCount (const Teuchos::ArrayView<const char>& imports,
                  const size_t offset,
                  const size_t numBytes)
  {
    if (numBytes == 0) {
      // Empty rows always take zero bytes, to ensure sparsity.
      return static_cast<size_t> (0);
    }
    else {
      LO numEntLO = 0;
      memcpy (&numEntLO, imports.getRawPtr () + offset, sizeof (LO));
      return static_cast<size_t> (numEntLO);
    }
  }

  // Return the number of bytes packed.
  template<class ST, class LO, class GO>
  size_t
  packRow (const Teuchos::ArrayView<char>& exports,
           const size_t offset,
           const size_t numEnt,
           const Teuchos::ArrayView<const GO>& gidsIn,
           const Teuchos::ArrayView<const int>& pidsIn,
           const Teuchos::ArrayView<const ST>& valsIn)
  {
    using Teuchos::ArrayView;

    if (numEnt == 0) {
      // Empty rows always take zero bytes, to ensure sparsity.
      return 0;
    }

    const size_t numEntBeg = offset;
    const size_t numEntLen = sizeof (LO);
    const size_t gidsBeg = numEntBeg + numEntLen;
    const size_t gidsLen = numEnt * sizeof (GO);
    const size_t pidsBeg = gidsBeg + gidsLen;
    const size_t pidsLen = numEnt * sizeof (int);
    const size_t valsBeg = pidsBeg + pidsLen;
    const size_t valsLen = numEnt * sizeof (ST);

    ArrayView<char> numEntOut = exports (numEntBeg, numEntLen);
    ArrayView<char> gidsOut = exports (gidsBeg, gidsLen);
    ArrayView<char> pidsOut = exports (pidsBeg, pidsLen);
    ArrayView<char> valsOut = exports (valsBeg, valsLen);

    size_t numBytesOut = 0;
    numBytesOut += packValue<LO> (numEntOut, static_cast<LO> (numEnt));
    numBytesOut += packArray<GO> (gidsOut, gidsIn, numEnt);
    numBytesOut += packArray<int> (pidsOut, pidsIn, numEnt);
    numBytesOut += packArray<ST> (valsOut, valsIn, numEnt);

    const size_t expectedNumBytes = numEntLen + gidsLen + pidsLen + valsLen;
    TEUCHOS_TEST_FOR_EXCEPTION(
      numBytesOut != expectedNumBytes, std::logic_error, "unpackRow: "
      "numBytesOut = " << numBytesOut << " != expectedNumBytes = "
      << expectedNumBytes << ".");

    return numBytesOut;
  }

  // Return the number of bytes actually read / used.
  template<class ST, class LO, class GO>
  size_t
  unpackRow (const Teuchos::ArrayView<GO>& gidsOut,
             const Teuchos::ArrayView<int>& pidsOut,
             const Teuchos::ArrayView<ST>& valsOut,
             const Teuchos::ArrayView<const char>& imports,
             const size_t offset,
             const size_t numBytes,
             const size_t numEnt)
  {
    using Teuchos::ArrayView;
    if (numBytes == 0) {
      // Rows with zero bytes always have zero entries.
      return 0;
    }
    TEUCHOS_TEST_FOR_EXCEPTION(
      static_cast<size_t> (imports.size ()) <= offset, std::logic_error,
      "unpackRow: imports.size() = " << imports.size () << " <= offset = "
      << offset << ".");
    TEUCHOS_TEST_FOR_EXCEPTION(
      static_cast<size_t> (imports.size ()) < offset + numBytes, std::logic_error,
      "unpackRow: imports.size() = " << imports.size () << " < offset + "
      "numBytes = " << (offset + numBytes) << ".");

    const size_t numEntBeg = offset;
    const size_t numEntLen = sizeof (LO);
    const size_t gidsBeg = numEntBeg + numEntLen;
    const size_t gidsLen = numEnt * sizeof (GO);
    const size_t pidsBeg = gidsBeg + gidsLen;
    const size_t pidsLen = numEnt * sizeof (int);
    const size_t valsBeg = pidsBeg + pidsLen;
    const size_t valsLen = numEnt * sizeof (ST);

    ArrayView<const char> numEntIn = imports (numEntBeg, numEntLen);
    ArrayView<const char> gidsIn = imports (gidsBeg, gidsLen);
    ArrayView<const char> pidsIn = imports (pidsBeg, pidsLen);
    ArrayView<const char> valsIn = imports (valsBeg, valsLen);

    size_t numBytesOut = 0;
    LO numEntOut;
    numBytesOut += unpackValue<LO> (numEntOut, numEntIn);
    TEUCHOS_TEST_FOR_EXCEPTION(
      static_cast<size_t> (numEntOut) != numEnt, std::logic_error,
      "unpackRow: Expected number of entries " << numEnt
      << " != actual number of entries " << numEntOut << ".");

    numBytesOut += unpackArray<GO> (gidsOut, gidsIn, numEnt);
    numBytesOut += unpackArray<int> (pidsOut, pidsIn, numEnt);
    numBytesOut += unpackArray<ST> (valsOut, valsIn, numEnt);

    TEUCHOS_TEST_FOR_EXCEPTION(
      numBytesOut != numBytes, std::logic_error, "unpackRow: numBytesOut = "
      << numBytesOut << " != numBytes = " << numBytes << ".");
    const size_t expectedNumBytes = numEntLen + gidsLen + pidsLen + valsLen;
    TEUCHOS_TEST_FOR_EXCEPTION(
      numBytesOut != expectedNumBytes, std::logic_error, "unpackRow: "
      "numBytesOut = " << numBytesOut << " != expectedNumBytes = "
      << expectedNumBytes << ".");
    return numBytesOut;
  }

} // namespace (anonymous)


namespace Tpetra {
namespace Import_Util {

template<typename Scalar,
         typename LocalOrdinal,
         typename GlobalOrdinal,
         typename Node>
void
packAndPrepareWithOwningPIDs (const CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>& SourceMatrix,
                              const Teuchos::ArrayView<const LocalOrdinal>& exportLIDs,
                              Teuchos::Array<char>& exports,
                              const Teuchos::ArrayView<size_t>& numPacketsPerLID,
                              size_t& constantNumPackets,
                              Distributor &distor,
                              const Teuchos::ArrayView<const int>& SourcePids)
{
  using Teuchos::Array;
  using Teuchos::ArrayView;
  using Teuchos::as;
  using Teuchos::RCP;
  typedef LocalOrdinal LO;
  typedef GlobalOrdinal GO;
  typedef Map<LocalOrdinal,GlobalOrdinal,Node> map_type;
  typedef CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> matrix_type;
  typedef typename matrix_type::impl_scalar_type ST;
  typedef typename ArrayView<const LO>::size_type size_type;
  const char prefix[] = "Tpetra::Import_Util::packAndPrepareWithOwningPIDs: ";

  // FIXME (mfh 03 Jan 2015) Currently, it might be the case that if a
  // graph or matrix owns no entries on a process, then it reports
  // that it is neither locally nor globally indexed, even if the
  // graph or matrix has a column Map.  This should change once we get
  // rid of lazy initialization in CrsGraph and CrsMatrix.
  TEUCHOS_TEST_FOR_EXCEPTION(
    ! SourceMatrix.isLocallyIndexed (), std::invalid_argument,
    prefix << "SourceMatrix must be locally indexed.");
  TEUCHOS_TEST_FOR_EXCEPTION(
    SourceMatrix.getColMap ().is_null (), std::logic_error,
    prefix << "The source matrix claims to be locally indexed, but its column "
    "Map is null.  This should never happen.  Please report this bug to the "
    "Tpetra developers.");
  const size_type numExportLIDs = exportLIDs.size ();
  TEUCHOS_TEST_FOR_EXCEPTION(
    numExportLIDs != numPacketsPerLID.size (), std::invalid_argument, prefix
    << "exportLIDs.size() = " << numExportLIDs << "!= numPacketsPerLID.size() "
    << " = " << numPacketsPerLID.size () << ".");
  TEUCHOS_TEST_FOR_EXCEPTION(
    static_cast<size_t> (SourcePids.size ()) != SourceMatrix.getColMap ()->getNodeNumElements (),
    std::invalid_argument, prefix << "SourcePids.size() = "
    << SourcePids.size ()
    << "!= SourceMatrix.getColMap()->getNodeNumElements() = "
    << SourceMatrix.getColMap ()->getNodeNumElements () << ".");

  // This tells the caller that different rows may have different
  // numbers of entries.  That is, the number of packets per LID might
  // not be a constant.
  constantNumPackets = 0;

  // Compute the number of bytes ("packets") per row to pack.  While
  // we're at it, compute the total number of matrix entries to send,
  // and the max number of entries in any of the rows we're sending.
  size_t totalNumBytes = 0;
  size_t totalNumEntries = 0;
  size_t maxRowLength = 0;
  for (size_type i = 0; i < numExportLIDs; ++i) {
    const LO lclRow = exportLIDs[i];
    const size_t numEnt = SourceMatrix.getNumEntriesInLocalRow (lclRow);
    const size_t numBytes = packRowCount<ST, LO, GO> (numEnt);
    numPacketsPerLID[i] = numBytes;
    totalNumBytes += numBytes;
    totalNumEntries += numEnt;
    maxRowLength = std::max (maxRowLength, numEnt);
  }

  // We use a "struct of arrays" approach to packing each row's
  // entries.  All the column indices (as global indices) go first,
  // then all their owning process ranks, and then the values.
  if (totalNumEntries > 0) {
    exports.resize (totalNumBytes);

    // Current position (in bytes) in the 'exports' output array.
    size_t offset = 0;

    // For each row of the matrix owned by the calling process, pack
    // that row's column indices and values into the exports array.

    // Locally indexed matrices always have a column Map.
    const map_type& colMap = * (SourceMatrix.getColMap ());

    // Temporary buffers for a copy of the column gids/pids
    Array<GO>  gids (maxRowLength);
    Array<int> pids (maxRowLength);

    const size_type numExportLIDs = exportLIDs.size();
    for (size_type i = 0; i < numExportLIDs; i++) {
      const LO lclRow = exportLIDs[i];

      // Get a locally indexed view of the current row's data.
      ArrayView<const Scalar> valsView;
      ArrayView<const LO> lidsView;
      SourceMatrix.getLocalRowView (lclRow, lidsView, valsView);
      ArrayView<const ST> valsViewST = av_reinterpret_cast<const ST> (valsView);

      // Convert column indices as LIDs to column indices as GIDs.
      const size_t numEnt = static_cast<size_t> (lidsView.size ());
      ArrayView<GO>  gidsView = gids (0, numEnt);
      ArrayView<int> pidsView = pids (0, numEnt);
      for (size_t k = 0; k < numEnt; ++k) {
        gidsView[k] = colMap.getGlobalElement (lidsView[k]);
        pidsView[k] = SourcePids[lidsView[k]];
      }

      // Copy the row's data into the current spot in the exports array.
      const size_t numBytes =
        packRow<ST, LO, GO> (exports, offset, numEnt,
                             gidsView, pidsView, valsViewST);
      // Keep track of how many bytes we packed.
      offset += numBytes;
    }

#ifdef HAVE_TPETRA_DEBUG
    TEUCHOS_TEST_FOR_EXCEPTION(
      offset != totalNumBytes, std::logic_error, prefix << "At end of method, "
      "the final offset (in bytes) " << offset << " does not equal the total "
      "number of bytes packed " << totalNumBytes << ".  Please report this bug "
      "to the Tpetra developers.");
#endif //  HAVE_TPETRA_DEBUG
  }
}

template<typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Node>
size_t
unpackAndCombineWithOwningPIDsCount (const CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> & SourceMatrix,
                                     const Teuchos::ArrayView<const LocalOrdinal> &importLIDs,
                                     const Teuchos::ArrayView<const char> &imports,
                                     const Teuchos::ArrayView<size_t> &numPacketsPerLID,
                                     size_t constantNumPackets,
                                     Distributor &distor,
                                     CombineMode combineMode,
                                     size_t numSameIDs,
                                     const ArrayView<const LocalOrdinal> &permuteToLIDs,
                                     const ArrayView<const LocalOrdinal> &permuteFromLIDs)
{
  typedef LocalOrdinal LO;
  typedef GlobalOrdinal GO;
  typedef CrsMatrix<Scalar, LO, GO, Node> matrix_type;
  typedef typename matrix_type::impl_scalar_type ST;
  typedef typename Teuchos::ArrayView<const LO>::size_type size_type;
  const char prefix[] = "unpackAndCombineWithOwningPIDsCount: ";

  TEUCHOS_TEST_FOR_EXCEPTION(
    permuteToLIDs.size () != permuteFromLIDs.size (), std::invalid_argument,
    prefix << "permuteToLIDs.size() = " << permuteToLIDs.size () << " != "
    "permuteFromLIDs.size() = " << permuteFromLIDs.size() << ".");
  // FIXME (mfh 26 Jan 2015) If there are no entries on the calling
  // process, then the matrix is neither locally nor globally indexed.
  const bool locallyIndexed = SourceMatrix.isLocallyIndexed ();
  TEUCHOS_TEST_FOR_EXCEPTION(
    ! locallyIndexed, std::invalid_argument, prefix << "The input CrsMatrix "
    "'SourceMatrix' must be locally indexed.");
  TEUCHOS_TEST_FOR_EXCEPTION(
    importLIDs.size () != numPacketsPerLID.size (), std::invalid_argument,
    prefix << "importLIDs.size() = " << importLIDs.size () << " != "
    "numPacketsPerLID.size() = " << numPacketsPerLID.size () << ".");

  // Number of matrix entries to unpack (returned by this function).
  size_t nnz = 0;

  // Count entries copied directly from the source matrix without permuting.
  for (size_t sourceLID = 0; sourceLID < numSameIDs; ++sourceLID) {
    const LO srcLID = static_cast<LO> (sourceLID);
    nnz += SourceMatrix.getNumEntriesInLocalRow (srcLID);
  }

  // Count entries copied directly from the source matrix with permuting.
  const size_type numPermuteToLIDs = permuteToLIDs.size ();
  for (size_type p = 0; p < numPermuteToLIDs; ++p) {
    nnz += SourceMatrix.getNumEntriesInLocalRow (permuteFromLIDs[p]);
  }

  // Count entries received from other MPI processes.
  size_t offset = 0;
  const size_type numImportLIDs = importLIDs.size ();
  for (size_type i = 0; i < numImportLIDs; ++i) {
    const size_t numBytes = numPacketsPerLID[i];
    const size_t numEnt = unpackRowCount<LO> (imports, offset, numBytes);
    nnz += numEnt;
    offset += numBytes;
  }
  return nnz;
}

template<typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Node>
void
unpackAndCombineIntoCrsArrays (const CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> & SourceMatrix,
                               const Teuchos::ArrayView<const LocalOrdinal>& importLIDs,
                               const Teuchos::ArrayView<const char>& imports,
                               const Teuchos::ArrayView<size_t>& numPacketsPerLID,
                               const size_t constantNumPackets,
                               Distributor& distor,
                               const CombineMode combineMode,
                               const size_t numSameIDs,
                               const Teuchos::ArrayView<const LocalOrdinal>& permuteToLIDs,
                               const Teuchos::ArrayView<const LocalOrdinal>& permuteFromLIDs,
                               size_t TargetNumRows,
                               size_t TargetNumNonzeros,
                               const int MyTargetPID,
                               const Teuchos::ArrayView<size_t>& CSR_rowptr,
                               const Teuchos::ArrayView<GlobalOrdinal>& CSR_colind,
                               const Teuchos::ArrayView<Scalar>& CSR_vals,
                               const Teuchos::ArrayView<const int>& SourcePids,
                               Teuchos::Array<int>& TargetPids)
{
  using Teuchos::ArrayView;
  using Teuchos::as;
  using Teuchos::av_reinterpret_cast;
  typedef LocalOrdinal LO;
  typedef GlobalOrdinal GO;
  typedef CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> matrix_type;
  typedef typename matrix_type::impl_scalar_type ST;
  typedef Map<LocalOrdinal,GlobalOrdinal,Node> map_type;
  typedef typename ArrayView<const LO>::size_type size_type;
  const char prefix[] = "Tpetra::Import_Util::unpackAndCombineIntoCrsArrays: ";

  const size_t N = TargetNumRows;
  const size_t mynnz = TargetNumNonzeros;
  // In the case of reduced communicators, the SourceMatrix won't have
  // the right "MyPID", so thus we have to supply it.
  const int MyPID = MyTargetPID;

  TEUCHOS_TEST_FOR_EXCEPTION(
    TargetNumRows + 1 != static_cast<size_t> (CSR_rowptr.size ()),
    std::invalid_argument, prefix << "CSR_rowptr.size() = " <<
    CSR_rowptr.size () << "!= TargetNumRows+1 = " << TargetNumRows+1 << ".");
  TEUCHOS_TEST_FOR_EXCEPTION(
    permuteToLIDs.size () != permuteFromLIDs.size (), std::invalid_argument,
    prefix << "permuteToLIDs.size() = " << permuteToLIDs.size ()
    << "!= permuteFromLIDs.size() = " << permuteFromLIDs.size () << ".");
  const size_type numImportLIDs = importLIDs.size ();
  TEUCHOS_TEST_FOR_EXCEPTION(
    numImportLIDs != numPacketsPerLID.size (), std::invalid_argument,
    prefix << "importLIDs.size() = " << numImportLIDs << " != "
    "numPacketsPerLID.size() = " << numPacketsPerLID.size() << ".");

  // Zero the rowptr
  for (size_t i = 0; i< N+1; ++i) {
    CSR_rowptr[i] = 0;
  }

  // SameIDs: Always first, always in the same place
  for (size_t i = 0; i < numSameIDs; ++i) {
    CSR_rowptr[i] = SourceMatrix.getNumEntriesInLocalRow (static_cast<LO> (i));
  }

  // PermuteIDs: Still local, but reordered
  size_t numPermuteIDs = permuteToLIDs.size ();
  for (size_t i = 0; i < numPermuteIDs; ++i) {
    CSR_rowptr[permuteToLIDs[i]] =
      SourceMatrix.getNumEntriesInLocalRow (permuteFromLIDs[i]);
  }

  // Setup CSR_rowptr for remotes
  const size_t totalNumBytes = imports.size ();
  {
    size_t offset = 0;
    for (size_type k = 0; k < numImportLIDs; ++k) {
      const size_t numBytes = numPacketsPerLID[k];
      const size_t numEnt = unpackRowCount<LO> (imports, offset, numBytes);
      CSR_rowptr[importLIDs[k]] += numEnt;
      offset += numBytes;
    }
  }

  // If multiple processes contribute to the same row, we may need to
  // update row offsets.  This tracks that.
  Teuchos::Array<size_t> NewStartRow (N + 1);

  // Turn row length into a real CSR_rowptr
  size_t last_len = CSR_rowptr[0];
  CSR_rowptr[0] = 0;
  for (size_t i = 1; i < N+1; ++i) {
    size_t new_len = CSR_rowptr[i];
    CSR_rowptr[i]  = last_len + CSR_rowptr[i-1];
    NewStartRow[i] = CSR_rowptr[i];
    last_len       = new_len;
  }

  TEUCHOS_TEST_FOR_EXCEPTION(
    CSR_rowptr[N] != mynnz, std::invalid_argument, prefix << "CSR_rowptr[last]"
    " = " << CSR_rowptr[N] << "!= mynnz = " << mynnz << ".");

  // Preseed TargetPids with -1 for local
  if (static_cast<size_t> (TargetPids.size ()) != mynnz) {
    TargetPids.resize (mynnz);
  }
  TargetPids.assign (mynnz, -1);

  // Grab pointers for SourceMatrix
  ArrayRCP<const size_t> Source_rowptr_RCP;
  ArrayRCP<const LO>     Source_colind_RCP;
  ArrayRCP<const Scalar> Source_vals_RCP;
  SourceMatrix.getAllValues (Source_rowptr_RCP, Source_colind_RCP, Source_vals_RCP);
  ArrayView<const size_t> Source_rowptr = Source_rowptr_RCP ();
  ArrayView<const LO>     Source_colind = Source_colind_RCP ();
  ArrayView<const Scalar> Source_vals = Source_vals_RCP ();

  const map_type& sourceColMap = * (SourceMatrix.getColMap());
  ArrayView<const GO> globalColElts = sourceColMap.getNodeElementList();

  // SameIDs: Copy the data over
  for (size_t i = 0; i < numSameIDs; ++i) {
    size_t FromRow = Source_rowptr[i];
    size_t ToRow   = CSR_rowptr[i];
    NewStartRow[i] += Source_rowptr[i+1] - Source_rowptr[i];

    for (size_t j = Source_rowptr[i]; j < Source_rowptr[i+1]; ++j) {
      CSR_vals[ToRow + j - FromRow]   = Source_vals[j];
      CSR_colind[ToRow + j - FromRow] = globalColElts[Source_colind[j]];
      TargetPids[ToRow + j - FromRow] =
        (SourcePids[Source_colind[j]] != MyPID) ?
        SourcePids[Source_colind[j]] : -1;
    }
  }

  // PermuteIDs: Copy the data over
  for (size_t i = 0; i < numPermuteIDs; ++i) {
    LO FromLID     = permuteFromLIDs[i];
    size_t FromRow = Source_rowptr[FromLID];
    size_t ToRow   = CSR_rowptr[permuteToLIDs[i]];

    NewStartRow[permuteToLIDs[i]] += Source_rowptr[FromLID+1]-Source_rowptr[FromLID];

    for (size_t j = Source_rowptr[FromLID]; j < Source_rowptr[FromLID+1]; ++j) {
      CSR_vals[ToRow + j - FromRow]   = Source_vals[j];
      CSR_colind[ToRow + j - FromRow] = globalColElts[Source_colind[j]];
      TargetPids[ToRow + j - FromRow] =
        (SourcePids[Source_colind[j]] != MyPID) ?
        SourcePids[Source_colind[j]] : -1;
    }
  }

  // RemoteIDs: Loop structure following UnpackAndCombine
  if (imports.size () > 0) {
    size_t offset = 0;
    int lclErr = 0;

    for (size_t i = 0; i < static_cast<size_t> (numImportLIDs); ++i) {
      const size_t numBytes = numPacketsPerLID[i];
      if (numBytes == 0) {
        // Empty buffer for that row means that the row is empty.
        continue;
      }
      const size_t numEnt   = unpackRowCount<LO> (imports, offset, numBytes);
      const LO lclRow       = importLIDs[i];
      const size_t StartRow = NewStartRow[lclRow];
      NewStartRow[lclRow]  += numEnt;

      ArrayView<GO> gidsOut = CSR_colind (StartRow, numEnt);
      ArrayView<int> pidsOut = TargetPids (StartRow, numEnt);
      ArrayView<Scalar> valsOutS = CSR_vals (StartRow, numEnt);
      ArrayView<ST> valsOut = av_reinterpret_cast<ST> (valsOutS);
      const size_t numBytesOut =
        unpackRow<ST, LO, GO> (gidsOut, pidsOut, valsOut, imports,
                               offset, numBytes, numEnt);
      if (numBytesOut != numBytes) {
        lclErr = 1;
        break;
      }
      // Correct target PIDs.
      for (size_t j = 0; j < numEnt; ++j) {
        const int pid = pidsOut[j];
        pidsOut[j] = (pid != MyPID) ? pid : -1;
      }
      offset += numBytes;
    }
#ifdef HAVE_TPETRA_DEBUG
    TEUCHOS_TEST_FOR_EXCEPTION(
      offset != totalNumBytes, std::logic_error, prefix << "After "
      "unpacking and counting all the imports, the final offset in bytes "
      << offset << " != total number of bytes " << totalNumBytes << ".  "
      "Please report this bug to the Tpetra developers.");
    TEUCHOS_TEST_FOR_EXCEPTION(
      lclErr != 0, std::logic_error, prefix << "numBytes != numBytesOut "
      "somewhere in unpack loop.  This should never happen.  "
      "Please report this bug to the Tpetra developers.");
#endif // HAVE_TPETRA_DEBUG
  }
}

// Note: This should get merged with the other Tpetra sort routines eventually.
template<typename Scalar, typename Ordinal>
void
sortCrsEntries (const Teuchos::ArrayView<size_t> &CRS_rowptr,
                const Teuchos::ArrayView<Ordinal> & CRS_colind,
                const Teuchos::ArrayView<Scalar> &CRS_vals)
{
  // For each row, sort column entries from smallest to largest.
  // Use shell sort. Stable sort so it is fast if indices are already sorted.
  // Code copied from  Epetra_CrsMatrix::SortEntries()
  size_t NumRows = CRS_rowptr.size()-1;
  size_t nnz = CRS_colind.size();

  for(size_t i = 0; i < NumRows; i++){
    size_t start=CRS_rowptr[i];
    if(start >= nnz) continue;

    Scalar* locValues   = &CRS_vals[start];
    size_t NumEntries   = CRS_rowptr[i+1] - start;
    Ordinal* locIndices = &CRS_colind[start];

    Ordinal n = NumEntries;
    Ordinal m = n/2;

    while(m > 0) {
      Ordinal max = n - m;
      for(Ordinal j = 0; j < max; j++) {
        for(Ordinal k = j; k >= 0; k-=m) {
          if(locIndices[k+m] >= locIndices[k])
            break;
          Scalar dtemp = locValues[k+m];
          locValues[k+m] = locValues[k];
          locValues[k] = dtemp;
          Ordinal itemp = locIndices[k+m];
          locIndices[k+m] = locIndices[k];
          locIndices[k] = itemp;
        }
      }
      m = m/2;
    }
  }
}

// Note: This should get merged with the other Tpetra sort routines eventually.
template<typename Scalar, typename Ordinal>
void
sortAndMergeCrsEntries (const Teuchos::ArrayView<size_t> &CRS_rowptr,
                        const Teuchos::ArrayView<Ordinal> & CRS_colind,
                        const Teuchos::ArrayView<Scalar> &CRS_vals)
{
  // For each row, sort column entries from smallest to largest,
  // merging column ids that are identify by adding values.  Use shell
  // sort. Stable sort so it is fast if indices are already sorted.
  // Code copied from Epetra_CrsMatrix::SortEntries()

  size_t NumRows = CRS_rowptr.size()-1;
  size_t nnz = CRS_colind.size();
  size_t new_curr=CRS_rowptr[0], old_curr=CRS_rowptr[0];

  for(size_t i = 0; i < NumRows; i++){
    size_t start=CRS_rowptr[i];
    if(start >= nnz) continue;

    Scalar* locValues   = &CRS_vals[start];
    size_t NumEntries   = CRS_rowptr[i+1] - start;
    Ordinal* locIndices = &CRS_colind[start];

    // Sort phase
    Ordinal n = NumEntries;
    Ordinal m = n/2;

    while(m > 0) {
      Ordinal max = n - m;
      for(Ordinal j = 0; j < max; j++) {
        for(Ordinal k = j; k >= 0; k-=m) {
          if(locIndices[k+m] >= locIndices[k])
            break;
          Scalar dtemp = locValues[k+m];
          locValues[k+m] = locValues[k];
          locValues[k] = dtemp;
          Ordinal itemp = locIndices[k+m];
          locIndices[k+m] = locIndices[k];
          locIndices[k] = itemp;
        }
      }
      m = m/2;
    }

    // Merge & shrink
    for(size_t j=CRS_rowptr[i]; j < CRS_rowptr[i+1]; j++) {
      if(j > CRS_rowptr[i] && CRS_colind[j]==CRS_colind[new_curr-1]) {
        CRS_vals[new_curr-1] += CRS_vals[j];
      }
      else if(new_curr==j) {
        new_curr++;
      }
      else {
        CRS_colind[new_curr] = CRS_colind[j];
        CRS_vals[new_curr]   = CRS_vals[j];
        new_curr++;
      }
    }

    CRS_rowptr[i] = old_curr;
    old_curr=new_curr;
  }

  CRS_rowptr[NumRows] = new_curr;
}

template <typename LocalOrdinal, typename GlobalOrdinal, typename Node>
void
lowCommunicationMakeColMapAndReindex (const ArrayView<const size_t> &rowptr,
                                      const ArrayView<LocalOrdinal> &colind_LID,
                                      const ArrayView<GlobalOrdinal> &colind_GID,
                                      const Teuchos::RCP<const Tpetra::Map<LocalOrdinal,GlobalOrdinal,Node> >& domainMapRCP,
                                      const Teuchos::ArrayView<const int> &owningPIDs,
                                      Teuchos::Array<int> &remotePIDs,
                                      Teuchos::RCP<const Tpetra::Map<LocalOrdinal,GlobalOrdinal,Node> > & colMap)
{
  // The domainMap is an RCP because there is a shortcut for a
  // (common) special case to return the columnMap = domainMap.
  const Tpetra::Map<LocalOrdinal,GlobalOrdinal,Node> &domainMap = *domainMapRCP;

  // Scan all column indices and sort into two groups:
  // Local:  those whose GID matches a GID of the domain map on this processor and
  // Remote: All others.
  size_t numDomainElements = domainMap.getNodeNumElements();
  Teuchos::Array<bool> LocalGIDs;
  if (numDomainElements>0) LocalGIDs.resize(numDomainElements,false); // Assume domain GIDs are not local

  // In principle it is good to have RemoteGIDs and RemotGIDList be as
  // long as the number of remote GIDs on this processor, but this
  // would require two passes through the column IDs, so we make it
  // the max of 100 and the number of block rows.
  const size_t numMyRows =  rowptr.size()-1;
  int    hashsize        = Teuchos::as<int>(numMyRows);
  if (hashsize < 100) hashsize = 100;

  Tpetra::Details::HashTable<GlobalOrdinal,LocalOrdinal> RemoteGIDs(hashsize);
  Teuchos::Array<GlobalOrdinal> RemoteGIDList; RemoteGIDList.reserve(hashsize);
  Teuchos::Array<int> PIDList;                 PIDList.reserve(hashsize);

  // Here we start using the *LocalOrdinal* colind_LID array.  This is
  // safe even if both columnIndices arrays are actually the same
  // (because LocalOrdinal==GlobalOrdinal).  For *local* GID's set
  // colind_LID with with their LID in the domainMap.  For *remote*
  // GIDs, we set colind_LID with (numDomainElements+NumRemoteColGIDs)
  // before the increment of the remote count.  These numberings will
  // be separate because no local LID is greater than
  // numDomainElements.

  size_t NumLocalColGIDs = 0;
  LocalOrdinal NumRemoteColGIDs = 0;
  for(size_t i = 0; i < numMyRows; i++) {
    for(size_t j = rowptr[i]; j < rowptr[i+1]; j++) {
      GlobalOrdinal GID = colind_GID[j];
      // Check if GID matches a row GID
      LocalOrdinal LID = domainMap.getLocalElement(GID);
      if(LID != -1) {
        bool alreadyFound = LocalGIDs[LID];
        if (!alreadyFound) {
          LocalGIDs[LID] = true; // There is a column in the graph associated with this domain map GID
          NumLocalColGIDs++;
        }
        colind_LID[j] = LID;
      }
      else {
        LocalOrdinal hash_value=RemoteGIDs.get(GID);
        if(hash_value  == -1) { // This means its a new remote GID
          int PID = owningPIDs[j];
          TEUCHOS_TEST_FOR_EXCEPTION(PID==-1,std::invalid_argument, "lowCommunicationMakeColMapAndReindex: Cannot figure out if PID is owned.");
          colind_LID[j] = Teuchos::as<LocalOrdinal>(numDomainElements + NumRemoteColGIDs);
          RemoteGIDs.add(GID, NumRemoteColGIDs);
          RemoteGIDList.push_back(GID);
          PIDList.push_back(PID);
          NumRemoteColGIDs++;
        }
        else
          colind_LID[j] = Teuchos::as<LocalOrdinal>(numDomainElements + hash_value);
      }
    }
  }

  // Possible short-circuit:  If all domain map GIDs are present as column indices, then set ColMap=domainMap and quit
  if (domainMap.getComm()->getSize()==1) {
    // Sanity check: When there is one processor,there can be no remoteGIDs
    TEUCHOS_TEST_FOR_EXCEPTION(NumRemoteColGIDs!=0,std::runtime_error,"lowCommunicationMakeColMapAndReindex: Some column IDs are not in domainMap.");
    if (Teuchos::as<size_t>(NumLocalColGIDs)==numDomainElements) {
      // In this case, we just use the domainMap's indices, which is, not coincidently, what we clobbered colind with up above anyway.
      // No further reindexing is needed.
      colMap = domainMapRCP;
      return;
    }
  }

  // Now build the array containing column GIDs
  // Build back end, containing remote GIDs, first
  LocalOrdinal numMyCols = NumLocalColGIDs + NumRemoteColGIDs;
  Teuchos::Array<GlobalOrdinal> ColIndices;
  GlobalOrdinal * RemoteColIndices=0;
  if(numMyCols > 0) {
    ColIndices.resize(numMyCols);
    if(NumLocalColGIDs!=Teuchos::as<size_t>(numMyCols)) RemoteColIndices = &ColIndices[NumLocalColGIDs]; // Points to back half of ColIndices
  }

  for(LocalOrdinal i = 0; i < NumRemoteColGIDs; i++)
    RemoteColIndices[i] = RemoteGIDList[i];

  // Build permute array for *remote* reindexing.
  Teuchos::Array<LocalOrdinal> RemotePermuteIDs(NumRemoteColGIDs);
  for(LocalOrdinal i=0; i<NumRemoteColGIDs; i++) RemotePermuteIDs[i]=i;


  // Sort External column indices so that all columns coming from a given remote processor are contiguous
  // This is a sort with two auxillary arrays: RemoteColIndices and RemotePermuteIDs.
  Tpetra::sort3(PIDList.begin(),PIDList.end(),ColIndices.begin()+NumLocalColGIDs,RemotePermuteIDs.begin());

  // Stash the RemotePIDs
  // Note: If Teuchos::Array had a shrink_to_fit like std::vector, we'd call it here.
  remotePIDs = PIDList;

  // Sort external column indices so that columns from a given remote processor are not only contiguous
  // but also in ascending order. NOTE: I don't know if the number of externals associated
  // with a given remote processor is known at this point ... so I count them here.

  // NTS: Only sort the RemoteColIndices this time...
  LocalOrdinal StartCurrent = 0, StartNext = 1;
  while ( StartNext < NumRemoteColGIDs ) {
    if (PIDList[StartNext]==PIDList[StartNext-1]) StartNext++;
    else {
      Tpetra::sort2(ColIndices.begin()+NumLocalColGIDs+StartCurrent,ColIndices.begin()+NumLocalColGIDs+StartNext,RemotePermuteIDs.begin()+StartCurrent);
      StartCurrent = StartNext; StartNext++;
    }
  }
    Tpetra::sort2(ColIndices.begin()+NumLocalColGIDs+StartCurrent,ColIndices.begin()+NumLocalColGIDs+StartNext,RemotePermuteIDs.begin()+StartCurrent);

  // Reverse the permutation to get the information we actually care about
  Teuchos::Array<LocalOrdinal> ReverseRemotePermuteIDs(NumRemoteColGIDs);
  for(LocalOrdinal i=0; i<NumRemoteColGIDs; i++) ReverseRemotePermuteIDs[RemotePermuteIDs[i]]=i;

  // Build permute array for *local* reindexing.
  bool use_local_permute=false;
  Teuchos::Array<LocalOrdinal> LocalPermuteIDs(numDomainElements);

  // Now fill front end. Two cases:
  // (1) If the number of Local column GIDs is the same as the number of Local domain GIDs, we
  //     can simply read the domain GIDs into the front part of ColIndices, otherwise
  // (2) We step through the GIDs of the domainMap, checking to see if each domain GID is a column GID.
  //     we want to do this to maintain a consistent ordering of GIDs between the columns and the domain.
  Teuchos::ArrayView<const GlobalOrdinal> domainGlobalElements = domainMap.getNodeElementList();
  if(Teuchos::as<size_t>(NumLocalColGIDs) == numDomainElements) {
    if(NumLocalColGIDs > 0) {
      // Load Global Indices into first numMyCols elements column GID list
      std::copy(domainGlobalElements.begin(),domainGlobalElements.end(),ColIndices.begin());
    }
  }
  else {
    LocalOrdinal NumLocalAgain = 0;
    use_local_permute = true;
    for(size_t i = 0; i < numDomainElements; i++) {
      if(LocalGIDs[i]) {
        LocalPermuteIDs[i] = NumLocalAgain;
        ColIndices[NumLocalAgain++] = domainGlobalElements[i];
      }
    }
    TEUCHOS_TEST_FOR_EXCEPTION(Teuchos::as<size_t>(NumLocalAgain)!=NumLocalColGIDs,std::runtime_error,"lowCommunicationMakeColMapAndReindex: Local ID count test failed.");
  }

  // Make Column map
  Tpetra::global_size_t minus_one = Teuchos::OrdinalTraits<Tpetra::global_size_t>::invalid();
  colMap = Teuchos::rcp(new Tpetra::Map<LocalOrdinal,GlobalOrdinal,Node>(minus_one,ColIndices,domainMap.getIndexBase(),domainMap.getComm(),domainMap.getNode()));

  // Low-cost reindex of the matrix
  for(size_t i=0; i<numMyRows; i++){
    for(size_t j=rowptr[i]; j<rowptr[i+1]; j++){
      LocalOrdinal ID=colind_LID[j];
      if(Teuchos::as<size_t>(ID) < numDomainElements){
        if(use_local_permute) colind_LID[j] = LocalPermuteIDs[colind_LID[j]];
        // In the case where use_local_permute==false, we just copy the DomainMap's ordering,
        // which it so happens is what we put in colind_LID to begin with.
      }
      else
        colind_LID[j] =  NumLocalColGIDs + ReverseRemotePermuteIDs[colind_LID[j]-numDomainElements];
    }
  }
}

} // namespace Import_Util
} // namespace Tpetra

// We can include the definitions for Tpetra::CrsMatrix now that the above
// functions have been defined.  For ETI, this isn't necessary, so we just
// including the generated hpp
#include "Tpetra_CrsMatrix.hpp"

#endif // TPETRA_IMPORT_UTIL_HPP
