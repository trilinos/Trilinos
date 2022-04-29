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

#ifndef TPETRA_ROWGRAPH_DEF_HPP
#define TPETRA_ROWGRAPH_DEF_HPP

namespace Tpetra {
  template<class LocalOrdinal, class GlobalOrdinal, class Node>
  void
  RowGraph<LocalOrdinal,GlobalOrdinal,Node>::
  pack (const Teuchos::ArrayView<const LocalOrdinal>& exportLIDs,
        Teuchos::Array<GlobalOrdinal>& exports,
        const Teuchos::ArrayView<size_t>& numPacketsPerLID,
        size_t& constantNumPackets) const
  {
    using Teuchos::Array;
    typedef LocalOrdinal LO;
    typedef GlobalOrdinal GO;
    typedef Map<LO, GO, Node> map_type;
    const char tfecfFuncName[] = "packAndPrepare";

    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
      exportLIDs.size() != numPacketsPerLID.size(), std::runtime_error,
      ": exportLIDs and numPacketsPerLID must have the same size.");

    const map_type& srcMap = * (this->getRowMap ());
    constantNumPackets = 0;

    // Set numPacketsPerLID[i] to the number of entries owned by the
    // calling process in (local) row exportLIDs[i] of the graph, that
    // the caller wants us to send out.  Compute the total number of
    // packets (that is, entries) owned by this process in all the
    // rows that the caller wants us to send out.
    size_t totalNumPackets = 0;
    nonconst_global_inds_host_view_type row;
    for (LO i = 0; i < exportLIDs.size (); ++i) {
      const GO GID = srcMap.getGlobalElement (exportLIDs[i]);
      size_t row_length = this->getNumEntriesInGlobalRow (GID);
      numPacketsPerLID[i] = row_length;
      totalNumPackets += row_length;
    }

    exports.resize (totalNumPackets);

    // Loop again over the rows to export, and pack rows of indices
    // into the output buffer.
    size_t exportsOffset = 0;
    for (LO i = 0; i < exportLIDs.size (); ++i) {
      const GO GID = srcMap.getGlobalElement (exportLIDs[i]);
      size_t row_length = this->getNumEntriesInGlobalRow (GID);
      Kokkos::resize(row,row_length);
      size_t check_row_length = 0;
      this->getGlobalRowCopy (GID, row, check_row_length);

      for (size_t j=0; j<row_length; ++j) {
        exports[exportsOffset+j] = row[j];
      }
      exportsOffset += row.extent(0);
    }
  }

#ifdef TPETRA_ENABLE_DEPRECATED_CODE
  template<class LocalOrdinal, class GlobalOrdinal, class Node>
  void
  RowGraph<LocalOrdinal,GlobalOrdinal,Node>::
  getLocalRowView (const LocalOrdinal /* lclRow */,
                   Teuchos::ArrayView<const LocalOrdinal>& /* lclColInds */) const
  {
    const char prefix[] = "Tpetra::RowGraph::getLocalRowView: ";

    // Subclasses are expected to implement this method.  We would
    // have made this method pure virtual, but that would have broken
    // backwards compatibility, since we added the method at least one
    // major release after introducing this class.
    TEUCHOS_TEST_FOR_EXCEPTION
      (! this->supportsRowViews (), std::runtime_error,
       prefix << "This object does not support row views.");
    TEUCHOS_TEST_FOR_EXCEPTION
      (this->supportsRowViews (), std::runtime_error,
       prefix << "This object claims to support row views, "
       "but this method is not implemented.");
  }

  template<class LocalOrdinal, class GlobalOrdinal, class Node>
  void
  RowGraph<LocalOrdinal,GlobalOrdinal,Node>::
  getGlobalRowView (const GlobalOrdinal /* gblRow */,
                    Teuchos::ArrayView<const GlobalOrdinal>& /* gblColInds */) const
  {
    const char prefix[] = "Tpetra::RowGraph::getGlobalRowView: ";

    // Subclasses are expected to implement this method.  We would
    // have made this method pure virtual, but that would have broken
    // backwards compatibility, since we added the method at least one
    // major release after introducing this class.
    TEUCHOS_TEST_FOR_EXCEPTION
      (! this->supportsRowViews (), std::runtime_error,
       prefix << "This object does not support row views.");
    TEUCHOS_TEST_FOR_EXCEPTION
      (this->supportsRowViews (), std::runtime_error,
       prefix << "This object claims to support row views, "
       "but this method is not implemented.");
  }
#endif // TPETRA_ENABLE_DEPRECATED_CODE

} // namespace Tpetra

//
// Explicit instantiation macro
//
// Must be expanded from within the Tpetra namespace!
//

#define TPETRA_ROWGRAPH_INSTANT(LO,GO,NODE) \
  template class RowGraph< LO , GO , NODE >;

#endif // TPETRA_ROWGRAPH_DEF_HPP
