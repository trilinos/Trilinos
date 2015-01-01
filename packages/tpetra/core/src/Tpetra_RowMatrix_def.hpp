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

#ifndef TPETRA_ROWMATRIX_DEF_HPP
#define TPETRA_ROWMATRIX_DEF_HPP

#include <Tpetra_ConfigDefs.hpp>
#include <Tpetra_CrsMatrix.hpp>
#include <Tpetra_Map.hpp>

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
    using Teuchos::ArrayRCP;
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
      const LO localNumRows = static_cast<LO> (A_rowMap->getNodeNumElements ());
      ArrayRCP<size_t> C_maxNumEntriesPerRow (localNumRows, 0);

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
        C = rcp (new crs_matrix_type (C_rowMap, C_maxNumEntriesPerRow,
                                      StaticProfile));
      } else {
        C = rcp (new crs_matrix_type (C_rowMap, C_maxNumEntriesPerRow,
                                      StaticProfile, constructorSublist));
      }
      // Since A and B have the same row Maps, we could add them
      // together all at once and merge values before we call
      // insertGlobalValues.  However, we don't really need to, since
      // we've already allocated enough space in each row of C for C
      // to do the merge itself.
    }
    else { // the row Maps of A and B are not the same
      // Construct the result matrix C.
      if (constructorSublist.is_null ()) {
        C = rcp (new crs_matrix_type (C_rowMap, 0, DynamicProfile));
      } else {
        C = rcp (new crs_matrix_type (C_rowMap, 0, DynamicProfile,
                                      constructorSublist));
      }
    }

#ifdef HAVE_TPETRA_DEBUG
    TEUCHOS_TEST_FOR_EXCEPTION(C.is_null (), std::logic_error,
      "Tpetra::RowMatrix::add: C should not be null at this point.  "
      "Please report this bug to the Tpetra developers.");
#endif // HAVE_TPETRA_DEBUG
    //
    // Compute C = alpha*A + beta*B.
    //
    Array<GO> ind;
    Array<Scalar> val;

    if (alpha != STS::zero ()) {
      const LO A_localNumRows = static_cast<LO> (A_rowMap->getNodeNumElements ());
      for (LO localRow = 0; localRow < A_localNumRows; ++localRow) {
        size_t A_numEntries = A.getNumEntriesInLocalRow (localRow);
        const GO globalRow = A_rowMap->getGlobalElement (localRow);
        if (A_numEntries > static_cast<size_t> (ind.size ())) {
          ind.resize (A_numEntries);
          val.resize (A_numEntries);
        }
        ArrayView<GO> indView = ind (0, A_numEntries);
        ArrayView<Scalar> valView = val (0, A_numEntries);
        A.getGlobalRowCopy (globalRow, indView, valView, A_numEntries);

        if (alpha != STS::one ()) {
          for (size_t k = 0; k < A_numEntries; ++k) {
            valView[k] *= alpha;
          }
        }
        C->insertGlobalValues (globalRow, indView, valView);
      }
    }

    if (beta != STS::zero ()) {
      const LO B_localNumRows = static_cast<LO> (B_rowMap->getNodeNumElements ());
      for (LO localRow = 0; localRow < B_localNumRows; ++localRow) {
        size_t B_numEntries = B.getNumEntriesInLocalRow (localRow);
        const GO globalRow = B_rowMap->getGlobalElement (localRow);
        if (B_numEntries > static_cast<size_t> (ind.size ())) {
          ind.resize (B_numEntries);
          val.resize (B_numEntries);
        }
        ArrayView<GO> indView = ind (0, B_numEntries);
        ArrayView<Scalar> valView = val (0, B_numEntries);
        B.getGlobalRowCopy (globalRow, indView, valView, B_numEntries);

        if (beta != STS::one ()) {
          for (size_t k = 0; k < B_numEntries; ++k) {
            valView[k] *= beta;
          }
        }
        C->insertGlobalValues (globalRow, indView, valView);
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
        size_t& constantNumPackets,
        Distributor &distor) const
  {
    using Teuchos::Array;
    using Teuchos::ArrayView;
    using Teuchos::as;
    using Teuchos::av_reinterpret_cast;
    using Teuchos::RCP;
    typedef LocalOrdinal LO;
    typedef GlobalOrdinal GO;
    typedef typename ArrayView<const LO>::size_type size_type;
    typedef Map<LO, GO, Node> map_type;
    const char tfecfFuncName[] = "pack";

    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
      exportLIDs.size() != numPacketsPerLID.size(),
      std::invalid_argument, "exportLIDs.size() = " << exportLIDs.size()
      << "!= numPacketsPerLID.size() = " << numPacketsPerLID.size() << ".");

    // Row Map of the source matrix.
    RCP<const map_type> rowMapPtr = this->getRowMap ();
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
      rowMapPtr.is_null (), std::runtime_error,
      "The source object's row Map is null.");
    const map_type& rowMap = *rowMapPtr;

    // Each input LID corresponds to a different row of the sparse
    // matrix.  In general, a RowMatrix doesn't have the same number
    // of entries in each row.
    constantNumPackets = 0;

    // Get the GIDs of the rows we want to pack.
    Array<GO> exportGIDs (exportLIDs.size ());
    const size_type numExportGIDs = exportGIDs.size ();
    for (size_type i = 0; i < numExportGIDs; ++i) {
      exportGIDs[i] = rowMap.getGlobalElement (exportLIDs[i]);
    }

    // We say "Packet" is char (really a "byte"), but the actual unit
    // of packing is a (GID, value) pair.  The GID is the column index
    // in that row of the sparse matrix, and the value is the value at
    // that entry of the sparse matrix.  Thus, we have to scale
    // numPacketsPerLID by the number of bytes in a _packed_ (GID,
    // value) pair.  (We pack the GID and value in each pair
    // separately, so the number of bytes in a packed pair is actually
    // sizeof(GO) + sizeof(Scalar).)
    //
    // FIXME (mfh 24 Feb 2013) This code is only correct if
    // sizeof(Scalar) is a meaningful representation of the amount of
    // data in a Scalar instance.  (GO is always a built-in integer
    // type.)
    //
    // Compute the number of packets per export LID, and accumulate
    // the total number of packages.  While doing so, find the max
    // number of entries in each row owned by this process; we will
    // use that to size temporary arrays below.
    const size_t sizeOfOrdValPair = sizeof (GO) + sizeof (Scalar);
    size_t totalNumEntries = 0;
    size_t maxRowLength = 0;
    for (size_type i = 0; i < exportGIDs.size(); ++i) {
      const size_t curNumEntries =
        this->getNumEntriesInGlobalRow (exportGIDs[i]);
      numPacketsPerLID[i] = curNumEntries * sizeOfOrdValPair;
      totalNumEntries += curNumEntries;
      maxRowLength = std::max (curNumEntries, maxRowLength);
    }

    // Pack export data by interleaving rows' indices and values in
    // the following way:
    //
    // [inds_row0 vals_row0 inds_row1 vals_row1 ... ]
    if (totalNumEntries > 0) {
      // exports is an array of char (bytes), so scale the total
      // number of entries by the number of bytes per entry (where
      // "entry" includes both the column index and the value).
      const size_t totalNumBytes = totalNumEntries * sizeOfOrdValPair;
      exports.resize (totalNumBytes);

      // Temporary buffers for a copy of the entries in each row.
      Array<GO> inds (as<size_type> (maxRowLength));
      Array<Scalar> vals (as<size_type> (maxRowLength));
      // Current position in the 'exports' output array.
      size_t curOffsetInBytes = 0;

      // For each row of the matrix owned by the calling process, pack
      // that row's column indices and values into the exports array.
      // We'll use getGlobalRowCopy, since it will always work, no
      // matter how the subclass stores its data.  Subclasses (like
      // CrsMatrix) should reimplement this method to use (e.g.,)
      // getGlobalRowView if appropriate, since that will likely be
      // more efficient.
      //
      // FIXME (mfh 28 Jun 2013) This could be made a (shared-memory)
      // parallel kernel, by using the CSR data layout to calculate
      // positions in the output buffer.
      for (size_type i = 0; i < exportGIDs.size(); ++i) {
        // Get a copy of the current row's data.
        size_t curNumEntries = 0;
        this->getGlobalRowCopy (exportGIDs[i], inds (), vals (), curNumEntries);
        // inds and vals arrays might have more entries than the row.
        // curNumEntries is the number of valid entries of these arrays.
        ArrayView<const GO> curInds = inds (0, as<size_type> (curNumEntries));
        ArrayView<const Scalar> curVals = vals (0, as<size_type> (curNumEntries));

        // Get views of the spots in the exports array in which to
        // put the indices resp. values.  See notes and FIXME above.
        ArrayView<char> outIndsChar =
          exports (curOffsetInBytes, curNumEntries * sizeof (GO));
        ArrayView<char> outValsChar =
          exports (curOffsetInBytes + curNumEntries * sizeof (GO),
                   curNumEntries * sizeof (Scalar));

        // Cast the above views of char as views of GO resp. Scalar.
        ArrayView<GO> outInds = av_reinterpret_cast<GO> (outIndsChar);
        ArrayView<Scalar> outVals = av_reinterpret_cast<Scalar> (outValsChar);

        // Copy the source matrix's row data into the views of the
        // exports array for indices resp. values.
        std::copy (curInds.begin(), curInds.end(), outInds.begin());
        std::copy (curVals.begin(), curVals.end(), outVals.begin());
        curOffsetInBytes += sizeOfOrdValPair * curNumEntries;
      }

#ifdef HAVE_TPETRA_DEBUG
      TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(curOffsetInBytes != totalNumBytes,
        std::logic_error, ": At end of method, the final offset bytes count "
        "curOffsetInBytes=" << curOffsetInBytes << " does not equal the total "
        "number of bytes packed totalNumBytes=" << totalNumBytes << ".  Please "
        "report this bug to the Tpetra developers.");
#endif //  HAVE_TPETRA_DEBUG
    }
  }

} // namespace Tpetra

//
// Explicit instantiation macro
//
// Must be expanded from within the Tpetra namespace!
//

#define TPETRA_ROWMATRIX_INSTANT(SCALAR,LO,GO,NODE) \
  \
  template class RowMatrix< SCALAR , LO , GO , NODE >;


#endif // TPETRA_ROWMATRIX_DEF_HPP

