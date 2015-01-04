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

/*!
  \file Tpetra_Import_Util.hpp
  \brief Utility functions and macros designed for use with Tpetra::Import and Tpetra::Export objects.
*/

#include "Tpetra_ConfigDefs.hpp" // for map, vector, string, and iostream
#include "Tpetra_Import.hpp"
#include "Tpetra_HashTable.hpp"
#include "Tpetra_Map.hpp"
#include "Tpetra_Util.hpp"
#include "Tpetra_Distributor.hpp"
#include "Tpetra_CrsMatrix.hpp"
#include <Teuchos_Array.hpp>
#include <utility>

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
                              const Teuchos::ArrayView<int>& SourcePids);

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

template <typename Matrix>
struct MatrixSerializationTraits {
  typedef typename Matrix::impl_scalar_type impl_scalar_type;

  static inline
  size_t scalarSize (const Matrix& mat) { return sizeof (impl_scalar_type); }

  static inline void
  packBuffer (const Matrix& mat,
              const size_t numEntries,
              const Teuchos::ArrayView<const impl_scalar_type>& vals,
              const Teuchos::ArrayView<char> packed_vals)
  {
    Teuchos::ArrayView<impl_scalar_type> packed_vals_scalar =
      Teuchos::av_reinterpret_cast<impl_scalar_type> (packed_vals);
    std::copy (vals.begin (), vals.begin () + numEntries,
               packed_vals_scalar.begin());
  }

  static inline void
  unpackScalar (const Matrix& mat,
                const char* val_char,
                impl_scalar_type& val)
  {
    val = * (reinterpret_cast<const impl_scalar_type*> (val_char));
  }
};

} // namespace Import_Util
} // namespace Tpetra

//
// Implementations
//

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
                              const Teuchos::ArrayView<int>& SourcePids)
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
  typedef MatrixSerializationTraits<matrix_type> serialization_type;
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
    exportLIDs.size () != numPacketsPerLID.size (),
    std::invalid_argument, prefix << "exportLIDs.size() = "
    << exportLIDs.size () << "!= numPacketsPerLID.size() = "
    << numPacketsPerLID.size() << ".");
  TEUCHOS_TEST_FOR_EXCEPTION(
    static_cast<size_t> (SourcePids.size ()) != SourceMatrix.getColMap ()->getNodeNumElements (),
    std::invalid_argument, prefix << "SourcePids.size() = "
    << SourcePids.size ()
    << "!= SourceMatrix.getColMap()->getNodeNumElements() = "
    << SourceMatrix.getColMap ()->getNodeNumElements () << ".");

  // Get a reference to the matrix's row Map.
  const map_type& rowMap = * (SourceMatrix.getRowMap ());
  constantNumPackets = 0;

  // Get the GIDs of the rows we want to pack.
  Array<GO> exportGIDs (exportLIDs.size ());
  const size_type numExportGIDs = exportGIDs.size ();
  for (size_type i = 0; i < numExportGIDs; ++i) {
    exportGIDs[i] = rowMap.getGlobalElement (exportLIDs[i]);
  }

  // Compute the size of a packet.  Each packet contains a matrix
  // entry (global column index and value) and its owning process
  // rank.
  const size_t scalar_size = serialization_type::scalarSize (SourceMatrix);
  const size_t sizeOfPacket = sizeof(GO) + scalar_size + sizeof(int);
  // Record the number of packets in each row.
  size_t totalNumEntries = 0;
  size_t maxRowLength = 0;
  for (size_type i = 0; i < exportGIDs.size(); ++i) {
    const size_t curNumEntries = SourceMatrix.getNumEntriesInGlobalRow(exportGIDs[i]);
    numPacketsPerLID[i] = curNumEntries * sizeOfPacket;
    totalNumEntries += curNumEntries;
    maxRowLength = std::max (curNumEntries, maxRowLength);
  }

  // Each packet contains an entry's global column index, the entry's
  // owning process rank, and the entry's value.  We pack export data
  // by storing packets contiguously, like this:
  //
  // [[inds_row0 pids_row0 vals_row0] [inds_row1 pids_row1 vals_row1] ... ]
  if (totalNumEntries > 0) {
    const size_t totalNumBytes = totalNumEntries * sizeOfPacket;
    exports.resize(totalNumBytes);

    // Current position in the 'exports' output array.
    size_t curOffsetInBytes = 0;

    // For each row of the matrix owned by the calling process, pack
    // that row's column indices and values into the exports array.

    // Locally indexed matrices always have a column Map.
    const map_type& colMap = * (SourceMatrix.getColMap ());
    ArrayView<const LocalOrdinal> lidsView;
    ArrayView<const Scalar> valsView;

    // Temporary buffers for a copy of the column gids/pids
    Array<GO>  gids (static_cast<size_type> (maxRowLength));
    Array<int> pids (static_cast<size_type> (maxRowLength));

    const size_type numExportLIDs = exportLIDs.size();
    for (size_type i = 0; i < numExportLIDs; i++) {
      // Get a (locally indexed) view of the current row's data.
      SourceMatrix.getLocalRowView (exportLIDs[i], lidsView, valsView);

      // Convert column indices as LIDs to column indices as GIDs.
      const size_type curNumEntries = lidsView.size ();
      size_t curNumEntriesST = static_cast<size_t> (curNumEntries);
      ArrayView<GO>  gidsView = gids (0, curNumEntries);
      ArrayView<int> pidsView = pids (0, curNumEntries);
      for (size_type k = 0; k < curNumEntries; ++k) {
        gidsView[k] = colMap.getGlobalElement (lidsView[k]);
        pidsView[k] = SourcePids[lidsView[k]];
      }

      // Views of the right places in each array so everthing looks
      // like the right data type
      ArrayView<char> gidsViewOutChar =
        exports (curOffsetInBytes, curNumEntriesST*sizeof(GO));
      ArrayView<char> pidsViewOutChar =
        exports (curOffsetInBytes+curNumEntriesST*sizeof(GO),
                 curNumEntriesST*sizeof(int));
      ArrayView<char> valsViewOutChar =
        exports (curOffsetInBytes+curNumEntriesST*(sizeof(GO)+sizeof(int)),
                 curNumEntriesST*scalar_size);
      ArrayView<GO> gidsViewOut  = av_reinterpret_cast<GO> (gidsViewOutChar);
      ArrayView<int> pidsViewOut = av_reinterpret_cast<int> (pidsViewOutChar);

      // Copy the row's data into the views of the exports array.
      std::copy (gidsView.begin (), gidsView.begin () + curNumEntriesST,
                 gidsViewOut.begin ());
      std::copy (pidsView.begin (), pidsView.begin () + curNumEntriesST,
                 pidsViewOut.begin ());
      ArrayView<const ST> valsViewST = av_reinterpret_cast<const ST> (valsView);
      serialization_type::packBuffer (SourceMatrix, curNumEntriesST, valsViewST, valsViewOutChar);

      // Keep track of how many bytes we packed.
      curOffsetInBytes += sizeOfPacket * curNumEntries;
    }

#ifdef HAVE_TPETRA_DEBUG
    TEUCHOS_TEST_FOR_EXCEPTION(
      curOffsetInBytes != totalNumBytes, std::logic_error, prefix << "At end "
      "of method, the final offset bytes count curOffsetInBytes=" <<
      curOffsetInBytes << " does not equal the total number of bytes packed "
      "totalNumBytes=" << totalNumBytes <<
      ".  Please report this bug to the Tpetra developers.");
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
  typedef typename ArrayView<const LO>::size_type size_type;
  typedef CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> matrix_type;
  typedef MatrixSerializationTraits<matrix_type> serialization_type;
  size_t nnz = 0;

  // CopyAndPermuteSection
  TEUCHOS_TEST_FOR_EXCEPTION(permuteToLIDs.size() != permuteFromLIDs.size(),
                             std::invalid_argument, "unpackAndCombineWithOwningPIDsCount: permuteToLIDs.size() = " << permuteToLIDs.size()
                             << "!= permuteFromLIDs.size() = " << permuteFromLIDs.size() << ".");
  const bool locallyIndexed = SourceMatrix.isLocallyIndexed();
  TEUCHOS_TEST_FOR_EXCEPTION(!locallyIndexed,std::invalid_argument, "unpackAndCombineWithOwningPIDsCount: SourceMatrix must be locally indexed.");

  // Copy
  const LO numSameIDs_as_LID = Teuchos::as<LO>(numSameIDs);
  for (LO sourceLID = 0; sourceLID < numSameIDs_as_LID; sourceLID++)
    nnz+=SourceMatrix.getNumEntriesInLocalRow(sourceLID);

  // Permute
  const size_t numPermuteToLIDs = Teuchos::as<size_t>(permuteToLIDs.size());
  for (size_t p = 0; p < numPermuteToLIDs; p++)
    nnz+=SourceMatrix.getNumEntriesInLocalRow(permuteFromLIDs[p]);

  // UnpackAndCombine Section
  TEUCHOS_TEST_FOR_EXCEPTION(importLIDs.size() != numPacketsPerLID.size(),
                             std::invalid_argument, "unpackAndCombineWithOwningPIDsCount: importLIDs.size() = " << importLIDs.size()
                             << "!= numPacketsPerLID.size() = " << numPacketsPerLID.size() << ".");

  const size_t scalar_size = serialization_type::scalarSize( SourceMatrix );
  const size_t sizeOfPacket    = sizeof(GlobalOrdinal)  + sizeof(int) + scalar_size;

  size_t curOffsetInBytes = 0;
  for (size_type i = 0; i < importLIDs.size(); ++i) {
    const size_t rowSize = numPacketsPerLID[i] / sizeOfPacket;
    curOffsetInBytes += rowSize * sizeOfPacket;
    nnz +=rowSize;
  }
  return nnz;
}

template<typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Node>
void
unpackAndCombineIntoCrsArrays (const CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> & SourceMatrix,
                               const Teuchos::ArrayView<const LocalOrdinal> &importLIDs,
                               const Teuchos::ArrayView<const char> &imports,
                               const Teuchos::ArrayView<size_t> &numPacketsPerLID,
                               size_t constantNumPackets,
                               Distributor &distor,
                               CombineMode combineMode,
                               size_t numSameIDs,
                               const ArrayView<const LocalOrdinal> &permuteToLIDs,
                               const ArrayView<const LocalOrdinal> &permuteFromLIDs,
                               size_t TargetNumRows,
                               size_t TargetNumNonzeros,
                               int MyTargetPID,
                               const ArrayView<size_t> &CSR_rowptr,
                               const ArrayView<GlobalOrdinal> &CSR_colind,
                               const ArrayView<Scalar> &CSR_vals,
                               const Teuchos::ArrayView<const int> &SourcePids,
                               Teuchos::Array<int> &TargetPids)
{
  using Teuchos::ArrayView;
  using Teuchos::as;
  using Teuchos::av_reinterpret_cast;
  typedef LocalOrdinal LO;
  typedef GlobalOrdinal GO;
  typedef CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> matrix_type;
  typedef typename matrix_type::impl_scalar_type impl_scalar_type;
  typedef Map<LocalOrdinal,GlobalOrdinal,Node> map_type;
  typedef typename ArrayView<const LO>::size_type size_type;
  typedef MatrixSerializationTraits<matrix_type> serialization_type;
  const char prefix[] = "Tpetra::Import_Util::unpackAndCombineIntoCrsArrays: ";

  size_t N = TargetNumRows;
  size_t mynnz = TargetNumNonzeros;
  // In the case of reduced communicators, the SourceMatrix won't have
  // the right "MyPID", so thus we have to supply it.
  int MyPID = MyTargetPID;

  // Zero the rowptr
  TEUCHOS_TEST_FOR_EXCEPTION(
    TargetNumRows + 1 != static_cast<size_t> (CSR_rowptr.size ()),
    std::invalid_argument, prefix << "CSR_rowptr.size() = " <<
    CSR_rowptr.size () << "!= TargetNumRows+1 = " << TargetNumRows+1 << ".");
  for (size_t i = 0; i< N+1; ++i) {
    CSR_rowptr[i] = 0;
  }

  // SameIDs: Always first, always in the same place
  for (size_t i = 0; i < numSameIDs; ++i) {
    CSR_rowptr[i] = SourceMatrix.getNumEntriesInLocalRow (static_cast<LO> (i));
  }

  // PermuteIDs: Still local, but reordered
  TEUCHOS_TEST_FOR_EXCEPTION(
    permuteToLIDs.size () != permuteFromLIDs.size (), std::invalid_argument,
    prefix << "permuteToLIDs.size() = " << permuteToLIDs.size ()
    << "!= permuteFromLIDs.size() = " << permuteFromLIDs.size () << ".");
  size_t numPermuteIDs = permuteToLIDs.size ();
  for (size_t i = 0; i < numPermuteIDs; ++i) {
    CSR_rowptr[permuteToLIDs[i]] =
      SourceMatrix.getNumEntriesInLocalRow (permuteFromLIDs[i]);
  }

  // Setup CSR_rowptr for remotes
  TEUCHOS_TEST_FOR_EXCEPTION(
    importLIDs.size () != numPacketsPerLID.size (), std::invalid_argument,
    prefix << "importLIDs.size() = " << importLIDs.size () << "!= "
    "numPacketsPerLID.size() = " << numPacketsPerLID.size() << ".");

  const size_t scalar_size = serialization_type::scalarSize (SourceMatrix);
  const size_t sizeOfPacket = sizeof(GlobalOrdinal) + sizeof(int) + scalar_size;
  const size_t totalNumBytes = imports.size ();
  const size_t RemoteNumEntries = totalNumBytes / sizeOfPacket;
  for (size_type k = 0; k < importLIDs.size (); ++k) {
    const size_t rowSize = numPacketsPerLID[k] / sizeOfPacket;
    CSR_rowptr[importLIDs[k]] += rowSize;
  }

  // If multiple procs contribute to a row;
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
    TargetPids.resize(mynnz);
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
  if (RemoteNumEntries > 0) {
    // data packed as follows:
    // [inds_row0 pids_row0 vals_row0 inds_row1 pids_row1 vals_row1 ...]
    ArrayView<const char>   avIndsC, avPidsC, avValsC;
    ArrayView<const GO>     avInds;
    ArrayView<const int>    avPids;

    size_t curOffsetInBytes = 0;
    for (size_t i = 0; i < static_cast<size_t> (importLIDs.size ()); ++i) {
      const size_t rowSize = numPacketsPerLID[i] / sizeOfPacket;
      LO ToLID     = importLIDs[i];
      int StartRow = NewStartRow[ToLID];
      NewStartRow[ToLID]+=rowSize;
      if (rowSize == 0) {
        continue;
      }

      // Get views of the import (incoming data) buffers.
      avIndsC = imports (curOffsetInBytes, rowSize*sizeof(GO));
      avPidsC = imports (curOffsetInBytes+rowSize*sizeof(GO),
                         rowSize*sizeof(int));
      avValsC = imports (curOffsetInBytes+rowSize*(sizeof(GO)+sizeof(int)),
                         rowSize*scalar_size);

      avInds = av_reinterpret_cast<const GO> (avIndsC);
      avPids = av_reinterpret_cast<const int> (avPidsC);

      const char * avValsC_ptr = avValsC.getRawPtr ();

      for (size_t j = 0; j < rowSize; ++j) {
        // mfh 03 Jan 2015: In the Kokkos refactor version of Tpetra,
        // Scalar and impl_scalar_type might differ.  impl_scalar_type
        // is what the matrix uses internally.  I've made
        // unpackScalar() (see below) use impl_scalar_type, which is
        // why we need the reinterpret cast here.  The cast is trivial
        // (Scalar == impl_scalar_type) for built-in "plain old data"
        // types like double, float, and int.
        impl_scalar_type& valOut =
          reinterpret_cast<impl_scalar_type&> (CSR_vals[StartRow+j]);
        serialization_type::unpackScalar (SourceMatrix, avValsC_ptr, valOut);
        CSR_colind[StartRow + j] = avInds[j];
        TargetPids[StartRow + j] = (avPids[j] != MyPID) ? avPids[j] : -1;
        avValsC_ptr += scalar_size;
      }
      curOffsetInBytes += rowSize * sizeOfPacket;
    }
#ifdef HAVE_TPETRA_DEBUG
    TEUCHOS_TEST_FOR_EXCEPTION(
      curOffsetInBytes != totalNumBytes, std::logic_error, prefix << "After "
      "unpacking and counting all the imports, the final offset in bytes "
      "curOffsetInBytes=" << curOffsetInBytes << " != total number of bytes "
      "totalNumBytes=" << totalNumBytes << ".  "
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

#endif // TPETRA_IMPORT_UTIL_HPP
