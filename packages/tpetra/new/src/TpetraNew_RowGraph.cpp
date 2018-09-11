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

#include "TpetraNew_RowGraph.hpp"
#include "TpetraNew_Map.hpp"
#include "Tpetra_Distributor.hpp"
#include "Teuchos_Array.hpp"
#include "Teuchos_ArrayView.hpp"
#include "Teuchos_StandardCatchMacros.hpp"
#include <stdexcept>

namespace TpetraNew {
  
  void
  RowGraph::
  pack (const Teuchos::ArrayView<const local_ordinal_type>& exportLIDs,
        Teuchos::Array<global_ordinal_type>& exports,
        const Teuchos::ArrayView<size_t>& numPacketsPerLID,
        size_t& constantNumPackets,
        ::Tpetra::Distributor& /* distor */) const
  {
    using LO = local_ordinal_type;
    using GO = global_ordinal_type;    
    const char tfecfFuncName[] = "packAndPrepare: ";

    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC
      (exportLIDs.size() != numPacketsPerLID.size(), std::runtime_error,
       "exportLIDs.size() = " << exportLIDs.size ()
       << " != numPacketsPerLID.size() = " << numPacketsPerLID.size ()
       << ".");
    const Map& srcMap = * (this->getRowMap ());
    constantNumPackets = 0;

    // Set numPacketsPerLID[i] to the number of entries owned by the
    // calling process in (local) row exportLIDs[i] of the graph, that
    // the caller wants us to send out.  Compute the total number of
    // packets (that is, entries) owned by this process in all the
    // rows that the caller wants us to send out.
    size_t totalNumPackets = 0;

    LO max_row_length = 0;
    for (LO i = 0; i < exportLIDs.size (); ++i) {
      const LO LID = exportLIDs[i];
      const GO GID = srcMap.getGlobalIndex (LID);
      const LO row_length = this->getNumEntriesInGlobalRow (GID);
      max_row_length = row_length > max_row_length ? row_length : max_row_length;
      numPacketsPerLID[i] = row_length;
      totalNumPackets += size_t (row_length);
    }
    Teuchos::Array<GO> row (max_row_length);
    exports.resize (totalNumPackets);

    // Loop again over the rows to export, and pack rows of indices
    // into the output buffer.
    size_t exportsOffset = 0;
    LO lastBadLID = 0;
    bool badness = false;
    for (LO LID : exportLIDs) {
      const GO GID = srcMap.getGlobalIndex (LID);
      const LO row_length = this->getGlobalRowCopy (GID, row ());
      if (row_length > LO (row.size ())) {
	badness = true;
	lastBadLID = LID;
      }

      Teuchos::ArrayView<const GO> rowView = row.view (0, row_length);
      size_t j = 0;
      for (auto&& rowEnt : rowView) {
        exports[exportsOffset+j] = rowEnt;
	++j;
      }
      exportsOffset += size_t (row_length);
    }

    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC
      (badness, std::logic_error, "getGlobalRowCopy gives inconsistent "
       "return values.  Last bad LID on this process: " << lastBadLID << ".");
  }

  void
  RowGraph::
  getLocalRowView (const local_ordinal_type /* lclRow */,
                   Teuchos::ArrayView<const local_ordinal_type>& /* lclColInds */) const
  {
    const char tfecfFuncName[] = "getLocalRowView: ";

    // Subclasses are expected to implement this method.  We would
    // have made this method pure virtual, but that would have broken
    // backwards compatibility, since we added the method at least one
    // major release after introducing this class.
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC
      (! this->supportsRowViews (), std::runtime_error,
       "This object does not support row views.");
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC
      (this->supportsRowViews (), std::runtime_error,
       "This object claims to support row views, "
       "but this method is not implemented.");
  }

  void
  RowGraph::
  getGlobalRowView (const global_ordinal_type /* gblRow */,
                    Teuchos::ArrayView<const global_ordinal_type>& /* gblColInds */) const
  {
    const char tfecfFuncName[] = "getGlobalRowView: ";

    // Subclasses are expected to implement this method.  We would
    // have made this method pure virtual, but that would have broken
    // backwards compatibility, since we added the method at least one
    // major release after introducing this class.
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC
      (! this->supportsRowViews (), std::runtime_error,
       "This object does not support row views.");
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC
      (this->supportsRowViews (), std::runtime_error,
       "This object claims to support row views, "
       "but this method is not implemented.");
  }

} // namespace TpetraNew
