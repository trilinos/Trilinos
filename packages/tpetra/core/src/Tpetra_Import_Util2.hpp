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

#include "Tpetra_ConfigDefs.hpp"
#include "Tpetra_Details_PackTraits.hpp"
#include "Tpetra_Import.hpp"
#include "Tpetra_HashTable.hpp"
#include "Tpetra_Map.hpp"
#include "Tpetra_Util.hpp"
#include "Tpetra_Distributor.hpp"
#include "Tpetra_Details_reallocDualViewIfNeeded.hpp"
#include "Kokkos_DualView.hpp"
#include <Teuchos_Array.hpp>
#include <utility>

// Tpetra::CrsMatrix uses the functions below in its implementation.
// To avoid a circular include issue, only include the declarations
// for CrsMatrix.  We will include the definition after the functions
// here have been defined.
#include "Tpetra_CrsMatrix_decl.hpp"

namespace { // (anonymous)

  template<class T, class D>
  Kokkos::View<T*, D, Kokkos::MemoryUnmanaged>
  getNonconstView (const Teuchos::ArrayView<T>& x)
  {
    typedef Kokkos::View<T*, D, Kokkos::MemoryUnmanaged> view_type;
    typedef typename view_type::size_type size_type;
    const size_type numEnt = static_cast<size_type> (x.size ());
    return view_type (x.getRawPtr (), numEnt);
  }

  template<class T, class D>
    Kokkos::View<const T*, D, Kokkos::MemoryUnmanaged>
  getConstView (const Teuchos::ArrayView<const T>& x)
  {
    typedef Kokkos::View<const T*, D, Kokkos::MemoryUnmanaged> view_type;
    typedef typename view_type::size_type size_type;
    const size_type numEnt = static_cast<size_type> (x.size ());
    return view_type (x.getRawPtr (), numEnt);
  }

  // For a given Kokkos (execution or memory) space, return both its
  // execution space, and the corresponding host execution space.
  template<class Space>
  struct GetHostExecSpace {
    typedef typename Space::execution_space execution_space;
    typedef typename Kokkos::View<int*, execution_space>::HostMirror::execution_space host_execution_space;
  };

} // namespace (anonymous)

namespace Tpetra {
namespace Import_Util {

/// \brief Special version of Tpetra::CrsMatrix::unpackAndCombine
///   that also unpacks owning process ranks.
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
                                     const Teuchos::ArrayView<const size_t>& numPacketsPerLID,
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
                               const Teuchos::ArrayView<const size_t>& numPacketsPerLID,
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


template<typename rowptr_array_type, typename colind_array_type, typename vals_array_type>
void
sortCrsEntries (const rowptr_array_type& CRS_rowptr,
                const colind_array_type& CRS_colind,
                const vals_array_type& CRS_vals);

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
                                      const Teuchos::RCP<const Tpetra::Map<LocalOrdinal,GlobalOrdinal,Node> > & domainMap,
                                      const Teuchos::ArrayView<const int> &owningPids,
                                      Teuchos::Array<int> &remotePids,
                                      Teuchos::RCP<const Tpetra::Map<LocalOrdinal,GlobalOrdinal,Node> > & colMap);


} // namespace Import_Util
} // namespace Tpetra


//
// Implementations
//

namespace { // (anonymous)

  template<class LO, class D>
  size_t
  unpackRowCount (const typename Tpetra::Details::PackTraits<LO, D>::input_buffer_type& imports,
                  const size_t offset,
                  const size_t numBytes,
                  const size_t numBytesPerValue)
  {
    using Kokkos::subview;
    using Tpetra::Details::PackTraits;
    typedef typename PackTraits<LO, D>::input_buffer_type input_buffer_type;
    typedef typename input_buffer_type::size_type size_type;

    if (numBytes == 0) {
      // Empty rows always take zero bytes, to ensure sparsity.
      return static_cast<size_t> (0);
    }
    else {
      LO numEntLO = 0;
      const size_t theNumBytes = PackTraits<LO, D>::packValueCount (numEntLO);
#ifdef HAVE_TPETRA_DEBUG
      TEUCHOS_TEST_FOR_EXCEPTION(
        theNumBytes > numBytes, std::logic_error, "unpackRowCount: "
        "theNumBytes = " << theNumBytes << " < numBytes = " << numBytes
        << ".");
#endif // HAVE_TPETRA_DEBUG
      const std::pair<size_type, size_type> rng (offset, offset + theNumBytes);
      input_buffer_type inBuf = subview (imports, rng); // imports (offset, theNumBytes);
      const size_t actualNumBytes = PackTraits<LO, D>::unpackValue (numEntLO, inBuf);
      (void)actualNumBytes;
#ifdef HAVE_TPETRA_DEBUG
      TEUCHOS_TEST_FOR_EXCEPTION(
        theNumBytes > numBytes, std::logic_error, "unpackRowCount: "
        "actualNumBytes = " << actualNumBytes << " < numBytes = " << numBytes
        << ".");
#endif // HAVE_TPETRA_DEBUG
      return static_cast<size_t> (numEntLO);
    }
  }

  // Return the number of bytes actually read / used.
  template<class ST, class LO, class GO, class D>
  size_t
  unpackRow (const typename Tpetra::Details::PackTraits<GO, D>::output_array_type& gidsOut,
             const typename Tpetra::Details::PackTraits<int, D>::output_array_type& pidsOut,
             const typename Tpetra::Details::PackTraits<ST, D>::output_array_type& valsOut,
             const typename Tpetra::Details::PackTraits<int, D>::input_buffer_type& imports,
             const size_t offset,
             const size_t numBytes,
             const size_t numEnt,
             const size_t numBytesPerValue)
  {
    using Kokkos::subview;
    using Tpetra::Details::PackTraits;
    // NOTE (mfh 02 Feb 2015) This assumes that input_buffer_type is
    // the same, no matter what type we're unpacking.  It's a
    // reasonable assumption, given that we go through the trouble of
    // PackTraits just so that we can pack data of different types in
    // the same buffer.
    typedef typename PackTraits<LO, D>::input_buffer_type input_buffer_type;
    typedef typename input_buffer_type::size_type size_type;
    typedef typename std::pair<size_type, size_type> pair_type;

    if (numBytes == 0) {
      // Rows with zero bytes always have zero entries.
      return 0;
    }
    TEUCHOS_TEST_FOR_EXCEPTION(
      static_cast<size_t> (imports.dimension_0 ()) <= offset, std::logic_error,
      "unpackRow: imports.dimension_0() = " << imports.dimension_0 () <<
      " <= offset = " << offset << ".");
    TEUCHOS_TEST_FOR_EXCEPTION(
      static_cast<size_t> (imports.dimension_0 ()) < offset + numBytes,
      std::logic_error, "unpackRow: imports.dimension_0() = "
      << imports.dimension_0 () << " < offset + numBytes = "
      << (offset + numBytes) << ".");

    const GO gid = 0; // packValueCount wants this
    const LO lid = 0; // packValueCount wants this
    const int pid = 0; // packValueCount wants this

    const size_t numEntBeg = offset;
    const size_t numEntLen = PackTraits<LO, D>::packValueCount (lid);
    const size_t gidsBeg = numEntBeg + numEntLen;
    const size_t gidsLen = numEnt * PackTraits<GO, D>::packValueCount (gid);
    const size_t pidsBeg = gidsBeg + gidsLen;
    const size_t pidsLen = numEnt * PackTraits<int, D>::packValueCount (pid);
    const size_t valsBeg = pidsBeg + pidsLen;
    const size_t valsLen = numEnt * numBytesPerValue;

    input_buffer_type numEntIn = subview (imports, pair_type (numEntBeg, numEntBeg + numEntLen));
    input_buffer_type gidsIn = subview (imports, pair_type (gidsBeg, gidsBeg + gidsLen));
    input_buffer_type pidsIn = subview (imports, pair_type (pidsBeg, pidsBeg + pidsLen));
    input_buffer_type valsIn = subview (imports, pair_type (valsBeg, valsBeg + valsLen));

    size_t numBytesOut = 0;
    int errorCode = 0;
    LO numEntOut;
    numBytesOut += PackTraits<LO, D>::unpackValue (numEntOut, numEntIn);
    TEUCHOS_TEST_FOR_EXCEPTION(
      static_cast<size_t> (numEntOut) != numEnt, std::logic_error,
      "unpackRow: Expected number of entries " << numEnt
      << " != actual number of entries " << numEntOut << ".");

    {
      Kokkos::pair<int, size_t> p;
      p = PackTraits<GO, D>::unpackArray (gidsOut, gidsIn, numEnt);
      errorCode += p.first;
      numBytesOut += p.second;

      p = PackTraits<int, D>::unpackArray (pidsOut, pidsIn, numEnt);
      errorCode += p.first;
      numBytesOut += p.second;

      p = PackTraits<ST, D>::unpackArray (valsOut, valsIn, numEnt);
      errorCode += p.first;
      numBytesOut += p.second;
    }

    TEUCHOS_TEST_FOR_EXCEPTION(
      numBytesOut != numBytes, std::logic_error, "unpackRow: numBytesOut = "
      << numBytesOut << " != numBytes = " << numBytes << ".");

    const size_t expectedNumBytes = numEntLen + gidsLen + pidsLen + valsLen;
    TEUCHOS_TEST_FOR_EXCEPTION(
      numBytesOut != expectedNumBytes, std::logic_error, "unpackRow: "
      "numBytesOut = " << numBytesOut << " != expectedNumBytes = "
      << expectedNumBytes << ".");

    TEUCHOS_TEST_FOR_EXCEPTION(
      errorCode != 0, std::logic_error, "PackTraits<>::unpackArray "
      "returned one or more nonzero error codes");

    return numBytesOut;
  }

} // namespace (anonymous)


namespace Tpetra {
namespace Import_Util {

template<typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Node>
size_t
unpackAndCombineWithOwningPIDsCount (const CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> & SourceMatrix,
                                     const Teuchos::ArrayView<const LocalOrdinal> &importLIDs,
                                     const Teuchos::ArrayView<const char> &imports,
                                     const Teuchos::ArrayView<const size_t>& numPacketsPerLID,
                                     size_t constantNumPackets,
                                     Distributor &distor,
                                     CombineMode combineMode,
                                     size_t numSameIDs,
                                     const Teuchos::ArrayView<const LocalOrdinal>& permuteToLIDs,
                                     const Teuchos::ArrayView<const LocalOrdinal>& permuteFromLIDs)
{
  using Kokkos::MemoryUnmanaged;
  using Kokkos::View;
  typedef LocalOrdinal LO;
  typedef GlobalOrdinal GO;
  typedef CrsMatrix<Scalar, LO, GO, Node> matrix_type;
  typedef typename matrix_type::impl_scalar_type ST;
  typedef typename Teuchos::ArrayView<const LO>::size_type size_type;
  typedef typename Node::execution_space execution_space;
  typedef typename GetHostExecSpace<execution_space>::host_execution_space HES;
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

  View<const char*, HES, MemoryUnmanaged> importsK (imports.getRawPtr (),
                                                    imports.size ());

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
    // FIXME (mfh 07 Feb 2015) Ask the matrix (rather, one of its
    // values, if it has one) for the (possibly run-time-depenendent)
    // number of bytes of one of its entries.
    const size_t numEnt = unpackRowCount<LO, HES> (importsK, offset,
                                                   numBytes, sizeof (ST));
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
                               const Teuchos::ArrayView<const size_t>& numPacketsPerLID,
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
  using Tpetra::Details::PackTraits;
  using Kokkos::MemoryUnmanaged;
  using Kokkos::subview;
  using Kokkos::View;
  using Teuchos::ArrayRCP;
  using Teuchos::ArrayView;
  using Teuchos::as;
  using Teuchos::av_reinterpret_cast;
  using Teuchos::outArg;
  using Teuchos::REDUCE_MAX;
  using Teuchos::reduceAll;
  typedef LocalOrdinal LO;
  typedef GlobalOrdinal GO;
  typedef typename Node::execution_space execution_space;
  typedef typename GetHostExecSpace<execution_space>::host_execution_space HES;
  typedef CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> matrix_type;
  typedef typename matrix_type::impl_scalar_type ST;
  typedef Map<LocalOrdinal,GlobalOrdinal,Node> map_type;
  typedef typename ArrayView<const LO>::size_type size_type;
  //typedef std::pair<typename Kokkos::View<int*, HES>::size_type,
  //                  typename Kokkos::View<int*, HES>::size_type> pair_type;
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

  // Kokkos::View of the input buffer of bytes to unpack.
  View<const char*, HES, MemoryUnmanaged> importsK (imports.getRawPtr (),
                                                    imports.size ());
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
  {
    size_t offset = 0;
    for (size_type k = 0; k < numImportLIDs; ++k) {
      const size_t numBytes = numPacketsPerLID[k];
      // FIXME (mfh 07 Feb 2015) Ask the matrix (rather, one of its
      // values, if it has one) for the (possibly run-time -
      // depenendent) number of bytes of one of its entries.
      const size_t numEnt = unpackRowCount<LO, HES> (importsK, offset,
                                                     numBytes, sizeof (ST));
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

  size_t numBytesPerValue = 0;
  if (PackTraits<ST, HES>::compileTimeSize) {
    ST val; // assume that ST is default constructible
    numBytesPerValue = PackTraits<ST, HES>::packValueCount (val);
  }
  else {
    // Since the packed data come from the source matrix, we can use
    // the source matrix to get the number of bytes per Scalar value
    // stored in the matrix.  This assumes that all Scalar values in
    // the source matrix require the same number of bytes.  If the
    // source matrix has no entries on the calling process, then we
    // have to ask the target matrix (via the output CSR arrays).  If
    // the target matrix has no entries on input on the calling
    // process, then hope that some process does have some idea how
    // big a Scalar value is.  Of course, if no processes have any
    // entries, then no values should be packed (though this does
    // assume that in our packing scheme, rows with zero entries take
    // zero bytes).
    if (Source_rowptr.size () == 0 || Source_rowptr[Source_rowptr.size () - 1] == 0) {
      numBytesPerValue = PackTraits<ST, HES>::packValueCount (CSR_vals[0]);
    }
    else {
      numBytesPerValue = PackTraits<ST, HES>::packValueCount (Source_vals[0]);
    }
    size_t lclNumBytesPerValue = numBytesPerValue;
    Teuchos::RCP<const Teuchos::Comm<int> > comm = SourceMatrix.getComm ();
    reduceAll<int, size_t> (*comm, REDUCE_MAX, lclNumBytesPerValue,
                            outArg (numBytesPerValue));
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
#ifdef HAVE_TPETRA_DEBUG
    int lclErr = 0;
#endif

    for (size_t i = 0; i < static_cast<size_t> (numImportLIDs); ++i) {
      const size_t numBytes = numPacketsPerLID[i];
      if (numBytes == 0) {
        // Empty buffer for that row means that the row is empty.
        continue;
      }
      // FIXME (mfh 07 Feb 2015) Ask the matrix (rather, one of its
      // values, if it has one) for the (possibly run-time -
      // depenendent) number of bytes of one of its entries.
      const size_t numEnt = unpackRowCount<LO, HES> (importsK, offset, numBytes,
                                                     numBytesPerValue);
      const LO lclRow = importLIDs[i];
      const size_t StartRow = NewStartRow[lclRow];
      NewStartRow[lclRow] += numEnt;

      View<GO*, HES, MemoryUnmanaged> gidsOut =
        getNonconstView<GO, HES> (CSR_colind (StartRow, numEnt));
      View<int*, HES, MemoryUnmanaged> pidsOut =
        getNonconstView<int, HES> (TargetPids (StartRow, numEnt));
      ArrayView<Scalar> valsOutS = CSR_vals (StartRow, numEnt);
      View<ST*, HES, MemoryUnmanaged> valsOut =
        getNonconstView<ST, HES> (av_reinterpret_cast<ST> (valsOutS));

      const size_t numBytesOut =
        unpackRow<ST, LO, GO, HES> (gidsOut, pidsOut, valsOut, importsK,
                                    offset, numBytes, numEnt, numBytesPerValue);
      if (numBytesOut != numBytes) {
#ifdef HAVE_TPETRA_DEBUG
        lclErr = 1;
#endif
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
      offset != static_cast<size_t> (imports.size ()), std::logic_error, prefix
      << "After unpacking and counting all the import packets, the final offset"
      " in bytes " << offset << " != imports.size() = " << imports.size () <<
      ".  Please report this bug to the Tpetra developers.");
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
    Ordinal m = 1;
    while (m<n) m = m*3+1;
    m /= 3;

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
      m = m/3;
    }
  }
}

template<typename rowptr_array_type, typename colind_array_type, typename vals_array_type>
void
sortCrsEntries (const rowptr_array_type& CRS_rowptr,
                const colind_array_type& CRS_colind,
                const vals_array_type& CRS_vals)
{
  // For each row, sort column entries from smallest to largest.
  // Use shell sort. Stable sort so it is fast if indices are already sorted.
  // Code copied from  Epetra_CrsMatrix::SortEntries()
  // NOTE: This should not be taken as a particularly efficient way to sort
  // rows of matrices in parallel.  But it is correct, so that's something.
  size_t NumRows = CRS_rowptr.dimension_0()-1;
  size_t nnz = CRS_colind.dimension_0();
  typedef typename colind_array_type::traits::non_const_value_type Ordinal;
  typedef typename vals_array_type::traits::non_const_value_type Scalar;

  typedef size_t index_type; // what this function was using; not my choice
  typedef typename vals_array_type::execution_space execution_space;
  // Specify RangePolicy explicitly, in order to use correct execution
  // space.  See #1345.
  typedef Kokkos::RangePolicy<execution_space, index_type> range_type;

  Kokkos::parallel_for (range_type (0, NumRows), KOKKOS_LAMBDA (const size_t i) {
      size_t start=CRS_rowptr(i);
      if(start < nnz) {
        size_t NumEntries = CRS_rowptr(i+1) - start;

        Ordinal n = (Ordinal) NumEntries;
        Ordinal m = 1;
        while (m<n) m = m*3+1;
        m /= 3;

        while(m > 0) {
          Ordinal max = n - m;
          for(Ordinal j = 0; j < max; j++) {
            for(Ordinal k = j; k >= 0; k-=m) {
              size_t sk = start+k;
              if(CRS_colind(sk+m) >= CRS_colind(sk))
                break;
              Scalar dtemp     = CRS_vals(sk+m);
              CRS_vals(sk+m)   = CRS_vals(sk);
              CRS_vals(sk)     = dtemp;
              Ordinal itemp    = CRS_colind(sk+m);
              CRS_colind(sk+m) = CRS_colind(sk);
              CRS_colind(sk)   = itemp;
            }
          }
          m = m/3;
        }
      }
    });
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
    size_t old_rowptr_i=CRS_rowptr[i];
    CRS_rowptr[i] = old_curr;
    if(old_rowptr_i >= nnz) continue;

    Scalar* locValues   = &CRS_vals[old_rowptr_i];
    size_t NumEntries   = CRS_rowptr[i+1] - old_rowptr_i;
    Ordinal* locIndices = &CRS_colind[old_rowptr_i];

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
    for(size_t j=old_rowptr_i; j < CRS_rowptr[i+1]; j++) {
      if(j > old_rowptr_i && CRS_colind[j]==CRS_colind[new_curr-1]) {
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
    old_curr=new_curr;
  }

  CRS_rowptr[NumRows] = new_curr;
}

template <typename LocalOrdinal, typename GlobalOrdinal, typename Node>
void
lowCommunicationMakeColMapAndReindex (const Teuchos::ArrayView<const size_t> &rowptr,
                                      const Teuchos::ArrayView<LocalOrdinal> &colind_LID,
                                      const Teuchos::ArrayView<GlobalOrdinal> &colind_GID,
                                      const Teuchos::RCP<const Tpetra::Map<LocalOrdinal,GlobalOrdinal,Node> >& domainMapRCP,
                                      const Teuchos::ArrayView<const int> &owningPIDs,
                                      Teuchos::Array<int> &remotePIDs,
                                      Teuchos::RCP<const Tpetra::Map<LocalOrdinal,GlobalOrdinal,Node> > & colMap)
{
  using Teuchos::rcp;
  typedef LocalOrdinal LO;
  typedef GlobalOrdinal GO;
  typedef Tpetra::global_size_t GST;
  typedef Tpetra::Map<LO, GO, Node> map_type;
  const char prefix[] = "lowCommunicationMakeColMapAndReindex: ";

  // The domainMap is an RCP because there is a shortcut for a
  // (common) special case to return the columnMap = domainMap.
  const map_type& domainMap = *domainMapRCP;

  // Scan all column indices and sort into two groups:
  // Local:  those whose GID matches a GID of the domain map on this processor and
  // Remote: All others.
  const size_t numDomainElements = domainMap.getNodeNumElements ();
  Teuchos::Array<bool> LocalGIDs;
  if (numDomainElements > 0) {
    LocalGIDs.resize (numDomainElements, false); // Assume domain GIDs are not local
  }

  // In principle it is good to have RemoteGIDs and RemotGIDList be as
  // long as the number of remote GIDs on this processor, but this
  // would require two passes through the column IDs, so we make it
  // the max of 100 and the number of block rows.
  //
  // FIXME (mfh 11 Feb 2015) Tpetra::Details::HashTable can hold at
  // most INT_MAX entries, but it's possible to have more rows than
  // that (if size_t is 64 bits and int is 32 bits).
  const size_t numMyRows = rowptr.size () - 1;
  const int hashsize = std::max (static_cast<int> (numMyRows), 100);

  Tpetra::Details::HashTable<GO, LO> RemoteGIDs (hashsize);
  Teuchos::Array<GO> RemoteGIDList;
  RemoteGIDList.reserve (hashsize);
  Teuchos::Array<int> PIDList;
  PIDList.reserve (hashsize);

  // Here we start using the *LocalOrdinal* colind_LID array.  This is
  // safe even if both columnIndices arrays are actually the same
  // (because LocalOrdinal==GO).  For *local* GID's set
  // colind_LID with with their LID in the domainMap.  For *remote*
  // GIDs, we set colind_LID with (numDomainElements+NumRemoteColGIDs)
  // before the increment of the remote count.  These numberings will
  // be separate because no local LID is greater than
  // numDomainElements.

  size_t NumLocalColGIDs = 0;
  LO NumRemoteColGIDs = 0;
  for (size_t i = 0; i < numMyRows; ++i) {
    for(size_t j = rowptr[i]; j < rowptr[i+1]; ++j) {
      const GO GID = colind_GID[j];
      // Check if GID matches a row GID
      const LO LID = domainMap.getLocalElement (GID);
      if(LID != -1) {
        const bool alreadyFound = LocalGIDs[LID];
        if (! alreadyFound) {
          LocalGIDs[LID] = true; // There is a column in the graph associated with this domain map GID
          NumLocalColGIDs++;
        }
        colind_LID[j] = LID;
      }
      else {
        const LO hash_value = RemoteGIDs.get (GID);
        if (hash_value == -1) { // This means its a new remote GID
          const int PID = owningPIDs[j];
          TEUCHOS_TEST_FOR_EXCEPTION(
            PID == -1, std::invalid_argument, prefix << "Cannot figure out if "
            "PID is owned.");
          colind_LID[j] = static_cast<LO> (numDomainElements + NumRemoteColGIDs);
          RemoteGIDs.add (GID, NumRemoteColGIDs);
          RemoteGIDList.push_back (GID);
          PIDList.push_back (PID);
          NumRemoteColGIDs++;
        }
        else {
          colind_LID[j] = static_cast<LO> (numDomainElements + hash_value);
        }
      }
    }
  }

  // Possible short-circuit: If all domain map GIDs are present as
  // column indices, then set ColMap=domainMap and quit.
  if (domainMap.getComm ()->getSize () == 1) {
    // Sanity check: When there is only one process, there can be no
    // remoteGIDs.
    TEUCHOS_TEST_FOR_EXCEPTION(
      NumRemoteColGIDs != 0, std::runtime_error, prefix << "There is only one "
      "process in the domain Map's communicator, which means that there are no "
      "\"remote\" indices.  Nevertheless, some column indices are not in the "
      "domain Map.");
    if (static_cast<size_t> (NumLocalColGIDs) == numDomainElements) {
      // In this case, we just use the domainMap's indices, which is,
      // not coincidently, what we clobbered colind with up above
      // anyway.  No further reindexing is needed.
      colMap = domainMapRCP;
      return;
    }
  }

  // Now build the array containing column GIDs
  // Build back end, containing remote GIDs, first
  const LO numMyCols = NumLocalColGIDs + NumRemoteColGIDs;
  Teuchos::Array<GO> ColIndices;
  GO* RemoteColIndices = NULL;
  if (numMyCols > 0) {
    ColIndices.resize (numMyCols);
    if (NumLocalColGIDs != static_cast<size_t> (numMyCols)) {
      RemoteColIndices = &ColIndices[NumLocalColGIDs]; // Points to back half of ColIndices
    }
  }

  for (LO i = 0; i < NumRemoteColGIDs; ++i) {
    RemoteColIndices[i] = RemoteGIDList[i];
  }

  // Build permute array for *remote* reindexing.
  Teuchos::Array<LO> RemotePermuteIDs (NumRemoteColGIDs);
  for (LO i = 0; i < NumRemoteColGIDs; ++i) {
    RemotePermuteIDs[i]=i;
  }

  // Sort External column indices so that all columns coming from a
  // given remote processor are contiguous.  This is a sort with two
  // auxillary arrays: RemoteColIndices and RemotePermuteIDs.
  Tpetra::sort3 (PIDList.begin (), PIDList.end (),
                 ColIndices.begin () + NumLocalColGIDs,
                 RemotePermuteIDs.begin ());

  // Stash the RemotePIDs.
  //
  // Note: If Teuchos::Array had a shrink_to_fit like std::vector,
  // we'd call it here.
  remotePIDs = PIDList;

  // Sort external column indices so that columns from a given remote
  // processor are not only contiguous but also in ascending
  // order. NOTE: I don't know if the number of externals associated
  // with a given remote processor is known at this point ... so I
  // count them here.

  // NTS: Only sort the RemoteColIndices this time...
  LO StartCurrent = 0, StartNext = 1;
  while (StartNext < NumRemoteColGIDs) {
    if (PIDList[StartNext]==PIDList[StartNext-1]) {
      StartNext++;
    }
    else {
      Tpetra::sort2 (ColIndices.begin () + NumLocalColGIDs + StartCurrent,
                     ColIndices.begin () + NumLocalColGIDs + StartNext,
                     RemotePermuteIDs.begin () + StartCurrent);
      StartCurrent = StartNext;
      StartNext++;
    }
  }
  Tpetra::sort2 (ColIndices.begin () + NumLocalColGIDs + StartCurrent,
                 ColIndices.begin () + NumLocalColGIDs + StartNext,
                 RemotePermuteIDs.begin () + StartCurrent);

  // Reverse the permutation to get the information we actually care about
  Teuchos::Array<LO> ReverseRemotePermuteIDs (NumRemoteColGIDs);
  for (LO i = 0; i < NumRemoteColGIDs; ++i) {
    ReverseRemotePermuteIDs[RemotePermuteIDs[i]] = i;
  }

  // Build permute array for *local* reindexing.
  bool use_local_permute = false;
  Teuchos::Array<LO> LocalPermuteIDs (numDomainElements);

  // Now fill front end. Two cases:
  //
  // (1) If the number of Local column GIDs is the same as the number
  //     of Local domain GIDs, we can simply read the domain GIDs into
  //     the front part of ColIndices, otherwise
  //
  // (2) We step through the GIDs of the domainMap, checking to see if
  //     each domain GID is a column GID.  we want to do this to
  //     maintain a consistent ordering of GIDs between the columns
  //     and the domain.
  Teuchos::ArrayView<const GO> domainGlobalElements = domainMap.getNodeElementList();
  if (static_cast<size_t> (NumLocalColGIDs) == numDomainElements) {
    if (NumLocalColGIDs > 0) {
      // Load Global Indices into first numMyCols elements column GID list
      std::copy (domainGlobalElements.begin (), domainGlobalElements.end (),
                 ColIndices.begin ());
    }
  }
  else {
    LO NumLocalAgain = 0;
    use_local_permute = true;
    for (size_t i = 0; i < numDomainElements; ++i) {
      if (LocalGIDs[i]) {
        LocalPermuteIDs[i] = NumLocalAgain;
        ColIndices[NumLocalAgain++] = domainGlobalElements[i];
      }
    }
    TEUCHOS_TEST_FOR_EXCEPTION(
      static_cast<size_t> (NumLocalAgain) != NumLocalColGIDs,
      std::runtime_error, prefix << "Local ID count test failed.");
  }

  // Make column Map
  const GST minus_one = Teuchos::OrdinalTraits<GST>::invalid ();
  colMap = rcp (new map_type (minus_one, ColIndices, domainMap.getIndexBase (),
                              domainMap.getComm (), domainMap.getNode ()));

  // Low-cost reindex of the matrix
  for (size_t i = 0; i < numMyRows; ++i) {
    for (size_t j = rowptr[i]; j < rowptr[i+1]; ++j) {
      const LO ID = colind_LID[j];
      if (static_cast<size_t> (ID) < numDomainElements) {
        if (use_local_permute) {
          colind_LID[j] = LocalPermuteIDs[colind_LID[j]];
        }
        // In the case where use_local_permute==false, we just copy
        // the DomainMap's ordering, which it so happens is what we
        // put in colind_LID to begin with.
      }
      else {
        colind_LID[j] =  NumLocalColGIDs + ReverseRemotePermuteIDs[colind_LID[j]-numDomainElements];
      }
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
