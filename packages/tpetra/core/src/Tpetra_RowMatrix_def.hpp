// @HEADER
// *****************************************************************************
//          Tpetra: Templated Linear Algebra Services Package
//
// Copyright 2008 NTESS and the Tpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef TPETRA_ROWMATRIX_DEF_HPP
#define TPETRA_ROWMATRIX_DEF_HPP

#include "Tpetra_CrsMatrix.hpp"
#include "Tpetra_Map.hpp"
#include "Tpetra_RowGraph.hpp"

namespace Tpetra {

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  RowMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node>::~RowMatrix() {}

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  Teuchos::RCP<RowMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> >
  RowMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
  add (const Scalar& alpha,
       const RowMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>& A,
       const Scalar& beta,
       const Teuchos::RCP<const Map<LocalOrdinal, GlobalOrdinal, Node> >& domainMap,
       const Teuchos::RCP<const Map<LocalOrdinal, GlobalOrdinal, Node> >& rangeMap,
       const Teuchos::RCP<Teuchos::ParameterList>& params) const
  {
    using Teuchos::Array;
    using Teuchos::ArrayView;
    using Teuchos::ParameterList;
    using Teuchos::RCP;
    using Teuchos::rcp;
    using Teuchos::rcp_implicit_cast;
    using Teuchos::sublist;
    typedef LocalOrdinal LO;
    typedef GlobalOrdinal GO;
    typedef Teuchos::ScalarTraits<Scalar> STS;
    typedef Map<LO, GO, Node> map_type;
    typedef RowMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> this_type;
    typedef CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> crs_matrix_type;

    const this_type& B = *this; // a convenient abbreviation

    // If the user didn't supply a domain or range Map, then try to
    // get one from B first (if it has them), then from A (if it has
    // them).  If we don't have any domain or range Maps, scold the
    // user.
    RCP<const map_type> A_domainMap = A.getDomainMap ();
    RCP<const map_type> A_rangeMap = A.getRangeMap ();
    RCP<const map_type> B_domainMap = B.getDomainMap ();
    RCP<const map_type> B_rangeMap = B.getRangeMap ();

    RCP<const map_type> theDomainMap = domainMap;
    RCP<const map_type> theRangeMap = rangeMap;

    if (domainMap.is_null ()) {
      if (B_domainMap.is_null ()) {
        TEUCHOS_TEST_FOR_EXCEPTION(
          A_domainMap.is_null (), std::invalid_argument,
          "Tpetra::RowMatrix::add: If neither A nor B have a domain Map, "
          "then you must supply a nonnull domain Map to this method.");
        theDomainMap = A_domainMap;
      } else {
        theDomainMap = B_domainMap;
      }
    }
    if (rangeMap.is_null ()) {
      if (B_rangeMap.is_null ()) {
        TEUCHOS_TEST_FOR_EXCEPTION(
          A_rangeMap.is_null (), std::invalid_argument,
          "Tpetra::RowMatrix::add: If neither A nor B have a range Map, "
          "then you must supply a nonnull range Map to this method.");
        theRangeMap = A_rangeMap;
      } else {
        theRangeMap = B_rangeMap;
      }
    }

#ifdef HAVE_TPETRA_DEBUG
    // In a debug build, check that A and B have matching domain and
    // range Maps, if they have domain and range Maps at all.  (If
    // they aren't fill complete, then they may not yet have them.)
    if (! A_domainMap.is_null () && ! A_rangeMap.is_null ()) {
      if (! B_domainMap.is_null () && ! B_rangeMap.is_null ()) {
        TEUCHOS_TEST_FOR_EXCEPTION(
          ! B_domainMap->isSameAs (*A_domainMap), std::invalid_argument,
          "Tpetra::RowMatrix::add: The input RowMatrix A must have a domain Map "
          "which is the same as (isSameAs) this RowMatrix's domain Map.");
        TEUCHOS_TEST_FOR_EXCEPTION(
          ! B_rangeMap->isSameAs (*A_rangeMap), std::invalid_argument,
          "Tpetra::RowMatrix::add: The input RowMatrix A must have a range Map "
          "which is the same as (isSameAs) this RowMatrix's range Map.");
        TEUCHOS_TEST_FOR_EXCEPTION(
          ! domainMap.is_null () && ! domainMap->isSameAs (*B_domainMap),
          std::invalid_argument,
          "Tpetra::RowMatrix::add: The input domain Map must be the same as "
          "(isSameAs) this RowMatrix's domain Map.");
        TEUCHOS_TEST_FOR_EXCEPTION(
          ! rangeMap.is_null () && ! rangeMap->isSameAs (*B_rangeMap),
          std::invalid_argument,
          "Tpetra::RowMatrix::add: The input range Map must be the same as "
          "(isSameAs) this RowMatrix's range Map.");
      }
    }
    else if (! B_domainMap.is_null () && ! B_rangeMap.is_null ()) {
      TEUCHOS_TEST_FOR_EXCEPTION(
        ! domainMap.is_null () && ! domainMap->isSameAs (*B_domainMap),
        std::invalid_argument,
        "Tpetra::RowMatrix::add: The input domain Map must be the same as "
        "(isSameAs) this RowMatrix's domain Map.");
      TEUCHOS_TEST_FOR_EXCEPTION(
        ! rangeMap.is_null () && ! rangeMap->isSameAs (*B_rangeMap),
        std::invalid_argument,
        "Tpetra::RowMatrix::add: The input range Map must be the same as "
        "(isSameAs) this RowMatrix's range Map.");
    }
    else {
      TEUCHOS_TEST_FOR_EXCEPTION(
        domainMap.is_null () || rangeMap.is_null (), std::invalid_argument,
        "Tpetra::RowMatrix::add: If neither A nor B have a domain and range "
        "Map, then you must supply a nonnull domain and range Map to this "
        "method.");
    }
#endif // HAVE_TPETRA_DEBUG

    // What parameters do we pass to C's constructor?  Do we call
    // fillComplete on C after filling it?  And if so, what parameters
    // do we pass to C's fillComplete call?
    bool callFillComplete = true;
    RCP<ParameterList> constructorSublist;
    RCP<ParameterList> fillCompleteSublist;
    if (! params.is_null ()) {
      callFillComplete = params->get ("Call fillComplete", callFillComplete);
      constructorSublist = sublist (params, "Constructor parameters");
      fillCompleteSublist = sublist (params, "fillComplete parameters");
    }

    RCP<const map_type> A_rowMap = A.getRowMap ();
    RCP<const map_type> B_rowMap = B.getRowMap ();
    RCP<const map_type> C_rowMap = B_rowMap; // see discussion in documentation
    RCP<crs_matrix_type> C; // The result matrix.

    // If A and B's row Maps are the same, we can compute an upper
    // bound on the number of entries in each row of C, before
    // actually computing the sum.  A reasonable upper bound is the
    // sum of the two entry counts in each row.  If we choose this as
    // the actual per-row upper bound, we can use static profile.
    if (A_rowMap->isSameAs (*B_rowMap)) {
      const LO localNumRows = static_cast<LO> (A_rowMap->getLocalNumElements ());
      Array<size_t> C_maxNumEntriesPerRow (localNumRows, 0);

      // Get the number of entries in each row of A.
      if (alpha != STS::zero ()) {
        for (LO localRow = 0; localRow < localNumRows; ++localRow) {
          const size_t A_numEntries = A.getNumEntriesInLocalRow (localRow);
          C_maxNumEntriesPerRow[localRow] += A_numEntries;
        }
      }
      // Get the number of entries in each row of B.
      if (beta != STS::zero ()) {
        for (LO localRow = 0; localRow < localNumRows; ++localRow) {
          const size_t B_numEntries = B.getNumEntriesInLocalRow (localRow);
          C_maxNumEntriesPerRow[localRow] += B_numEntries;
        }
      }
      // Construct the result matrix C.
      if (constructorSublist.is_null ()) {
        C = rcp (new crs_matrix_type (C_rowMap, C_maxNumEntriesPerRow ()));
      } else {
        C = rcp (new crs_matrix_type (C_rowMap, C_maxNumEntriesPerRow (),
                                      constructorSublist));
      }
      // Since A and B have the same row Maps, we could add them
      // together all at once and merge values before we call
      // insertGlobalValues.  However, we don't really need to, since
      // we've already allocated enough space in each row of C for C
      // to do the merge itself.
    }
    else { // the row Maps of A and B are not the same
      // Construct the result matrix C.
      // true: !A_rowMap->isSameAs (*B_rowMap)
      TEUCHOS_TEST_FOR_EXCEPTION(true,
				 std::invalid_argument,
				 "Tpetra::RowMatrix::add: The row maps must be the same for statically "
				 "allocated matrices in order to be sure that there is sufficient space "
				 "to do the addition");
    }

#ifdef HAVE_TPETRA_DEBUG
    TEUCHOS_TEST_FOR_EXCEPTION(C.is_null (), std::logic_error,
      "Tpetra::RowMatrix::add: C should not be null at this point.  "
      "Please report this bug to the Tpetra developers.");
#endif // HAVE_TPETRA_DEBUG
    //
    // Compute C = alpha*A + beta*B.
    //
    using gids_type = nonconst_global_inds_host_view_type;
    using vals_type = nonconst_values_host_view_type;
    gids_type ind;
    vals_type val;

    if (alpha != STS::zero ()) {
      const LO A_localNumRows = static_cast<LO> (A_rowMap->getLocalNumElements ());
      for (LO localRow = 0; localRow < A_localNumRows; ++localRow) {
        size_t A_numEntries = A.getNumEntriesInLocalRow (localRow);
        const GO globalRow = A_rowMap->getGlobalElement (localRow);
        if (A_numEntries > static_cast<size_t> (ind.size ())) {
          Kokkos::resize(ind,A_numEntries);
          Kokkos::resize(val,A_numEntries);
        }
        gids_type indView = Kokkos::subview(ind, std::make_pair((size_t)0, A_numEntries));
        vals_type valView = Kokkos::subview(val, std::make_pair((size_t)0, A_numEntries));
        A.getGlobalRowCopy (globalRow, indView, valView, A_numEntries);

        if (alpha != STS::one ()) {
          for (size_t k = 0; k < A_numEntries; ++k) {
            valView[k] *= alpha;
          }
        }
        C->insertGlobalValues (globalRow, A_numEntries, 
                               reinterpret_cast<const Scalar*>(valView.data()),
                               indView.data());
      }
    }

    if (beta != STS::zero ()) {
      const LO B_localNumRows = static_cast<LO> (B_rowMap->getLocalNumElements ());
      for (LO localRow = 0; localRow < B_localNumRows; ++localRow) {
        size_t B_numEntries = B.getNumEntriesInLocalRow (localRow);
        const GO globalRow = B_rowMap->getGlobalElement (localRow);
        if (B_numEntries > static_cast<size_t> (ind.size ())) {
          Kokkos::resize(ind,B_numEntries);
          Kokkos::resize(val,B_numEntries);
        }
        gids_type indView = Kokkos::subview(ind, std::make_pair((size_t)0, B_numEntries));
        vals_type valView = Kokkos::subview(val, std::make_pair((size_t)0, B_numEntries));
        B.getGlobalRowCopy (globalRow, indView, valView, B_numEntries);

        if (beta != STS::one ()) {
          for (size_t k = 0; k < B_numEntries; ++k) {
            valView[k] *= beta;
          }
        }
        C->insertGlobalValues (globalRow, B_numEntries, 
                               reinterpret_cast<const Scalar*>(valView.data()),
                               indView.data());
      }
    }

    if (callFillComplete) {
      if (fillCompleteSublist.is_null ()) {
        C->fillComplete (theDomainMap, theRangeMap);
      } else {
        C->fillComplete (theDomainMap, theRangeMap, fillCompleteSublist);
      }
    }

    return rcp_implicit_cast<this_type> (C);
  }


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void
  RowMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
  pack (const Teuchos::ArrayView<const LocalOrdinal>& exportLIDs,
        Teuchos::Array<char>& exports,
        const Teuchos::ArrayView<size_t>& numPacketsPerLID,
        size_t& constantNumPackets) const
  {
#ifdef HAVE_TPETRA_DEBUG
    const char tfecfFuncName[] = "pack: ";
    {
      using Teuchos::reduceAll;
      std::ostringstream msg;
      int lclBad = 0;
      try {
        this->packImpl (exportLIDs, exports, numPacketsPerLID,
                        constantNumPackets);
      } catch (std::exception& e) {
        lclBad = 1;
        msg << e.what ();
      }
      int gblBad = 0;
      const Teuchos::Comm<int>& comm = * (this->getComm ());
      reduceAll<int, int> (comm, Teuchos::REDUCE_MAX,
                           lclBad, Teuchos::outArg (gblBad));
      if (gblBad != 0) {
        const int myRank = comm.getRank ();
        const int numProcs = comm.getSize ();
        for (int r = 0; r < numProcs; ++r) {
          if (r == myRank && lclBad != 0) {
            std::ostringstream os;
            os << "Proc " << myRank << ": " << msg.str () << std::endl;
            std::cerr << os.str ();
          }
          comm.barrier ();
          comm.barrier ();
          comm.barrier ();
        }
        TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
          true, std::logic_error, "packImpl() threw an exception on one or "
          "more participating processes.");
      }
    }
#else
    this->packImpl (exportLIDs, exports, numPacketsPerLID,
                    constantNumPackets);
#endif // HAVE_TPETRA_DEBUG
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void
  RowMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
  allocatePackSpace (Teuchos::Array<char>& exports,
                     size_t& totalNumEntries,
                     const Teuchos::ArrayView<const LocalOrdinal>& exportLIDs) const
  {
    typedef LocalOrdinal LO;
    typedef GlobalOrdinal GO;
    typedef typename Teuchos::ArrayView<const LO>::size_type size_type;
    //const char tfecfFuncName[] = "allocatePackSpace: ";
    const size_type numExportLIDs = exportLIDs.size ();

    // Count the total number of entries to send.
    totalNumEntries = 0;
    for (size_type i = 0; i < numExportLIDs; ++i) {
      const LO lclRow = exportLIDs[i];
      size_t curNumEntries = this->getNumEntriesInLocalRow (lclRow);
      // FIXME (mfh 25 Jan 2015) We should actually report invalid row
      // indices as an error.  Just consider them nonowned for now.
      if (curNumEntries == Teuchos::OrdinalTraits<size_t>::invalid ()) {
        curNumEntries = 0;
      }
      totalNumEntries += curNumEntries;
    }

    // FIXME (mfh 24 Feb 2013) This code is only correct if
    // sizeof(Scalar) is a meaningful representation of the amount of
    // data in a Scalar instance.  (LO and GO are always built-in
    // integer types.)
    //
    // Allocate the exports array.  It does NOT need padding for
    // alignment, since we use memcpy to write to / read from send /
    // receive buffers.
    const size_t allocSize =
      static_cast<size_t> (numExportLIDs) * sizeof (LO) +
      totalNumEntries * (sizeof (Scalar) + sizeof (GO));
    if (static_cast<size_t> (exports.size ()) < allocSize) {
      exports.resize (allocSize);
    }
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  bool
  RowMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
  packRow (char* const numEntOut,
           char* const valOut,
           char* const indOut,
           const size_t numEnt,
           const LocalOrdinal lclRow) const
  {
    using Teuchos::Array;
    using Teuchos::ArrayView;
    typedef LocalOrdinal LO;
    typedef GlobalOrdinal GO;
    typedef Tpetra::Map<LO, GO, Node> map_type;

    const LO numEntLO = static_cast<LO> (numEnt);
    memcpy (numEntOut, &numEntLO, sizeof (LO));

    if (this->supportsRowViews ()) {
      if (this->isLocallyIndexed ()) {
        // If the matrix is locally indexed on the calling process, we
        // have to use its column Map (which it _must_ have in this
        // case) to convert to global indices.
        local_inds_host_view_type indIn;
        values_host_view_type valIn;
        this->getLocalRowView (lclRow, indIn, valIn);
        const map_type& colMap = * (this->getColMap ());
        // Copy column indices one at a time, so that we don't need
        // temporary storage.
        for (size_t k = 0; k < numEnt; ++k) {
          const GO gblIndIn = colMap.getGlobalElement (indIn[k]);
          memcpy (indOut + k * sizeof (GO), &gblIndIn, sizeof (GO));
        }
        memcpy (valOut, valIn.data (), numEnt * sizeof (Scalar));
      }
      else if (this->isGloballyIndexed ()) {
        // If the matrix is globally indexed on the calling process,
        // then we can use the column indices directly.  However, we
        // have to get the global row index.  The calling process must
        // have a row Map, since otherwise it shouldn't be participating
        // in packing operations.
        global_inds_host_view_type indIn;
        values_host_view_type valIn;
        const map_type& rowMap = * (this->getRowMap ());
        const GO gblRow = rowMap.getGlobalElement (lclRow);
        this->getGlobalRowView (gblRow, indIn, valIn);
        memcpy (indOut, indIn.data (), numEnt * sizeof (GO));
        memcpy (valOut, valIn.data (), numEnt * sizeof (Scalar));
      }
      else {
        if (numEnt != 0) {
          return false;
        }
      }
    }
    else {
      // FIXME (mfh 25 Jan 2015) Pass in valIn and indIn as scratch
      // space, instead of allocating them on each call.
      if (this->isLocallyIndexed ()) {
        nonconst_local_inds_host_view_type indIn("indIn",numEnt);
        nonconst_values_host_view_type valIn("valIn",numEnt);
        size_t theNumEnt = 0;
        this->getLocalRowCopy (lclRow, indIn, valIn, theNumEnt);
        if (theNumEnt != numEnt) {
          return false;
        }
        const map_type& colMap = * (this->getColMap ());
        // Copy column indices one at a time, so that we don't need
        // temporary storage.
        for (size_t k = 0; k < numEnt; ++k) {
          const GO gblIndIn = colMap.getGlobalElement (indIn[k]);
          memcpy (indOut + k * sizeof (GO), &gblIndIn, sizeof (GO));
        }
        memcpy (valOut, valIn.data(), numEnt * sizeof (Scalar));
      }
      else if (this->isGloballyIndexed ()) {
        nonconst_global_inds_host_view_type indIn("indIn",numEnt);
        nonconst_values_host_view_type valIn("valIn",numEnt);
        const map_type& rowMap = * (this->getRowMap ());
        const GO gblRow = rowMap.getGlobalElement (lclRow);
        size_t theNumEnt = 0;
        this->getGlobalRowCopy (gblRow, indIn, valIn, theNumEnt);
        if (theNumEnt != numEnt) {
          return false;
        }
        memcpy (indOut, indIn.data(), numEnt * sizeof (GO));
        memcpy (valOut, valIn.data(), numEnt * sizeof (Scalar));
      }
      else {
        if (numEnt != 0) {
          return false;
        }
      }
    }
    return true;
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void
  RowMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
  packImpl (const Teuchos::ArrayView<const LocalOrdinal>& exportLIDs,
            Teuchos::Array<char>& exports,
            const Teuchos::ArrayView<size_t>& numPacketsPerLID,
            size_t& constantNumPackets) const
  {
    using Teuchos::Array;
    using Teuchos::ArrayView;
    using Teuchos::as;
    using Teuchos::av_reinterpret_cast;
    using Teuchos::RCP;
    typedef LocalOrdinal LO;
    typedef GlobalOrdinal GO;
    typedef typename ArrayView<const LO>::size_type size_type;
    const char tfecfFuncName[] = "packImpl: ";

    const size_type numExportLIDs = exportLIDs.size ();
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
      numExportLIDs != numPacketsPerLID.size (), std::invalid_argument,
      "exportLIDs.size() = " << numExportLIDs << " != numPacketsPerLID.size()"
      " = " << numPacketsPerLID.size () << ".");

    // Setting this to zero tells the caller to expect a possibly
    // different ("nonconstant") number of packets per local index
    // (i.e., a possibly different number of entries per row).
    constantNumPackets = 0;

    // The pack buffer 'exports' enters this method possibly
    // unallocated.  Do the first two parts of "Count, allocate, fill,
    // compute."
    size_t totalNumEntries = 0;
    allocatePackSpace (exports, totalNumEntries, exportLIDs);
    const size_t bufSize = static_cast<size_t> (exports.size ());

    // Compute the number of "packets" (in this case, bytes) per
    // export LID (in this case, local index of the row to send), and
    // actually pack the data.
    //
    // FIXME (mfh 24 Feb 2013, 25 Jan 2015) This code is only correct
    // if sizeof(Scalar) is a meaningful representation of the amount
    // of data in a Scalar instance.  (LO and GO are always built-in
    // integer types.)

    // Variables for error reporting in the loop.
    size_type firstBadIndex = 0; // only valid if outOfBounds == true.
    size_t firstBadOffset = 0;   // only valid if outOfBounds == true.
    size_t firstBadNumBytes = 0; // only valid if outOfBounds == true.
    bool outOfBounds = false;
    bool packErr = false;

    char* const exportsRawPtr = exports.getRawPtr ();
    size_t offset = 0; // current index into 'exports' array.
    for (size_type i = 0; i < numExportLIDs; ++i) {
      const LO lclRow = exportLIDs[i];
      const size_t numEnt = this->getNumEntriesInLocalRow (lclRow);

      // Only pad this row if it has a nonzero number of entries.
      if (numEnt == 0) {
        numPacketsPerLID[i] = 0;
      }
      else {
        char* const numEntBeg = exportsRawPtr + offset;
        char* const numEntEnd = numEntBeg + sizeof (LO);
        char* const valBeg = numEntEnd;
        char* const valEnd = valBeg + numEnt * sizeof (Scalar);
        char* const indBeg = valEnd;
        const size_t numBytes = sizeof (LO) +
          numEnt * (sizeof (Scalar) + sizeof (GO));
        if (offset > bufSize || offset + numBytes > bufSize) {
          firstBadIndex = i;
          firstBadOffset = offset;
          firstBadNumBytes = numBytes;
          outOfBounds = true;
          break;
        }
        packErr = ! packRow (numEntBeg, valBeg, indBeg, numEnt, lclRow);
        if (packErr) {
          firstBadIndex = i;
          firstBadOffset = offset;
          firstBadNumBytes = numBytes;
          break;
        }
        // numPacketsPerLID[i] is the number of "packets" in the
        // current local row i.  Packet=char (really "byte") so use
        // the number of bytes of the packed data for that row.
        numPacketsPerLID[i] = numBytes;
        offset += numBytes;
      }
    }

    // The point of these exceptions is to stop computation if the
    // above checks found a bug.  If HAVE_TPETRA_DEBUG is defined,
    // Tpetra will do extra all-reduces to check for global
    // consistency of the error state.  Otherwise, throwing an
    // exception here might cause deadlock, but we accept that as
    // better than just continuing.
    TEUCHOS_TEST_FOR_EXCEPTION(
      outOfBounds, std::logic_error, "First invalid offset into 'exports' "
      "pack buffer at index i = " << firstBadIndex << ".  exportLIDs[i]: "
      << exportLIDs[firstBadIndex] << ", bufSize: " << bufSize << ", offset: "
      << firstBadOffset << ", numBytes: " << firstBadNumBytes << ".");
    TEUCHOS_TEST_FOR_EXCEPTION(
      packErr, std::logic_error, "First error in packRow() at index i = "
      << firstBadIndex << ".  exportLIDs[i]: " << exportLIDs[firstBadIndex]
      << ", bufSize: " << bufSize << ", offset: " << firstBadOffset
      << ", numBytes: " << firstBadNumBytes << ".");
  }


} // namespace Tpetra

//
// Explicit instantiation macro
//
// Must be expanded from within the Tpetra namespace!
//

#define TPETRA_ROWMATRIX_INSTANT(SCALAR,LO,GO,NODE) \
  template class RowMatrix< SCALAR , LO , GO , NODE >;


#endif // TPETRA_ROWMATRIX_DEF_HPP

