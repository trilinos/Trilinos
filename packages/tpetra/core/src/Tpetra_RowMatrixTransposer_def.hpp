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

#ifndef TPETRA_ROWMATRIXTRANSPOSER_DEF_HPP
#define TPETRA_ROWMATRIXTRANSPOSER_DEF_HPP

#include "Tpetra_CrsMatrix.hpp"
#include "Tpetra_Export.hpp"
#include "Tpetra_Import.hpp"
#include "Tpetra_Details_computeOffsets.hpp"
#include "Tpetra_Details_shortSort.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_TimeMonitor.hpp"

namespace Tpetra {

namespace { // (anonymous)

auto
view_no_init (const std::string& label) ->
  decltype (Kokkos::view_alloc (label, Kokkos::WithoutInitializing))
{
  return Kokkos::view_alloc (label, Kokkos::WithoutInitializing);
}

} // namespace (anonymous)

template<class Scalar,
         class LocalOrdinal,
         class GlobalOrdinal,
         class Node>
RowMatrixTransposer<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
RowMatrixTransposer (const Teuchos::RCP<const crs_matrix_type>& origMatrix,
                     const std::string& label)
  : origMatrix_ (origMatrix), label_ (label)
{}

template<class Scalar,
         class LocalOrdinal,
         class GlobalOrdinal,
         class Node>
Teuchos::RCP<CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> >
RowMatrixTransposer<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
createTranspose (const Teuchos::RCP<Teuchos::ParameterList> &params)
{
  using Teuchos::RCP;
  // Do the local transpose
  RCP<crs_matrix_type> transMatrixWithSharedRows = createTransposeLocal (params);

#ifdef HAVE_TPETRA_MMM_TIMINGS
  const std::string prefix = std::string ("Tpetra ") + label_ + ": ";
  using Teuchos::TimeMonitor;
  TimeMonitor MM (*TimeMonitor::getNewTimer (prefix + "Transpose TAFC"));
#endif

  // If transMatrixWithSharedRows has an exporter, that's what we
  // want.  If it doesn't, the rows aren't actually shared, and we're
  // done!
  using export_type = Export<LocalOrdinal, GlobalOrdinal, Node>;
  RCP<const export_type> exporter =
    transMatrixWithSharedRows->getGraph ()->getExporter ();
  if (exporter.is_null ()) {
    return transMatrixWithSharedRows;
  }
  else {
    Teuchos::ParameterList labelList;
#ifdef HAVE_TPETRA_MMM_TIMINGS
    labelList.set("Timer Label", label_);
#endif
    if(! params.is_null ()) {
      const char paramName[] = "compute global constants";
      labelList.set (paramName, params->get (paramName, true));
    }
    // Use the Export object to do a fused Export and fillComplete.
    return exportAndFillCompleteCrsMatrix<crs_matrix_type>
      (transMatrixWithSharedRows, *exporter, Teuchos::null,
       Teuchos::null, Teuchos::rcpFromRef (labelList));
  }
}


// mfh 03 Feb 2013: In a definition outside the class like this, the
// return value is considered outside the class scope (for things like
// resolving typedefs), but the arguments are considered inside the
// class scope.
template<class Scalar,
         class LocalOrdinal,
         class GlobalOrdinal,
         class Node>
Teuchos::RCP<CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> >
RowMatrixTransposer<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
createTransposeLocal (const Teuchos::RCP<Teuchos::ParameterList>& params)
{
  using Teuchos::Array;
  using Teuchos::ArrayRCP;
  using Teuchos::ArrayView;
  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::rcp_dynamic_cast;
  using LO = LocalOrdinal;
  using GO = GlobalOrdinal;
  using IST = typename CrsMatrix<Scalar, LO, GO, Node>::impl_scalar_type;
  using import_type = Tpetra::Import<LO, GO, Node>;
  using export_type = Tpetra::Export<LO, GO, Node>;

#ifdef HAVE_TPETRA_MMM_TIMINGS
  std::string prefix = std::string("Tpetra ") + label_ + ": ";
  using Teuchos::TimeMonitor;
  TimeMonitor MM (*TimeMonitor::getNewTimer (prefix + "Transpose Local"));
#endif

  const bool sortMerge = [&] () {
    constexpr bool sortMergeDefault = false; // see #4607 discussion
    const char sortParamName[] = "sort and merge";
    return params.get () == nullptr ? sortMergeDefault :
      params->get (sortParamName, sortMergeDefault);
  } ();

  // Prebuild the importers and exporters the no-communication way,
  // flipping the importers and exporters around.
  RCP<const import_type> myImport;
  RCP<const export_type> myExport;
  if (! origMatrix_->getGraph ()->getImporter ().is_null ()) {
    myExport = rcp (new export_type (*origMatrix_->getGraph ()->getImporter ()));
  }
  if (! origMatrix_->getGraph ()->getExporter ().is_null ()) {
    myImport = rcp (new import_type (*origMatrix_->getGraph ()->getExporter ()));
  }

  //
  // This transpose is based upon the approach in EpetraExt.
  //

  const LO numLocalCols (origMatrix_->getNodeNumCols ());
  const LO numLocalRows (origMatrix_->getNodeNumRows ());
  size_t numLocalNnz  = origMatrix_->getNodeNumEntries();

  RCP<const crs_matrix_type> crsMatrix =
    rcp_dynamic_cast<const crs_matrix_type> (origMatrix_);

  RCP<crs_matrix_type> transMatrixWithSharedRows;
  if (crsMatrix != Teuchos::null) {
    // The matrix is a CrsMatrix
    using local_matrix_type = typename crs_matrix_type::local_matrix_type;
    using local_graph_type = typename crs_matrix_type::local_graph_type;
    using offset_type = typename local_graph_type::size_type;
    using row_map_type = typename local_matrix_type::row_map_type::non_const_type;
    using index_type = typename local_matrix_type::index_type::non_const_type;
    using values_type = typename local_matrix_type::values_type::non_const_type;
    using execution_space = typename local_matrix_type::execution_space;

    local_matrix_type lclMatrix = crsMatrix->getLocalMatrix ();
    local_graph_type lclGraph = lclMatrix.graph;

    // Determine how many nonzeros there are per row in the transpose.
    using DT = typename crs_matrix_type::device_type;
    Kokkos::View<LO*, DT> t_counts ("transpose row counts", numLocalCols);
    using range_type = Kokkos::RangePolicy<LO, execution_space>;
    Kokkos::parallel_for
      ("Compute row counts of local transpose",
       range_type (0, numLocalRows),
       KOKKOS_LAMBDA (const LO row) {
        auto rowView = lclGraph.rowConst(row);
        const LO length  = rowView.length;

        for (LO colID = 0; colID < length; ++colID) {
          const LO col = rowView(colID);
          Kokkos::atomic_fetch_add (&t_counts[col], LO (1));
        }
    });

    row_map_type t_offsets
      (view_no_init ("transpose ptr"), numLocalCols + 1);

    // TODO (mfh 10 Jul 2019): This returns the sum of all counts,
    // which could be useful for checking numLocalNnz.
    (void) Details::computeOffsetsFromCounts (t_offsets, t_counts);

    index_type t_cols
      (view_no_init ("transpose lcl ind"), numLocalNnz);
    values_type t_vals
      (view_no_init ("transpose val"), numLocalNnz);
    Kokkos::parallel_for
      ("Compute local transpose",
       range_type (0, numLocalRows),
       KOKKOS_LAMBDA (const LO row) {
        auto rowView = lclMatrix.rowConst(row);
        const LO length  = rowView.length;

        for (LO colID = 0; colID < length; colID++) {
          const LO col = rowView.colidx(colID);
          const offset_type beg = t_offsets[col];
          const LO old_count =
            Kokkos::atomic_fetch_sub (&t_counts[col], LO (1));
          const LO len (t_offsets[col+1] - beg);
          const offset_type insert_pos = beg + (len - old_count);
          t_cols[insert_pos] = row;
          t_vals[insert_pos] = rowView.value(colID);
        }
    });

    // Invariant: At this point, all entries of t_counts are zero.
    // This means we can use it to store new post-merge counts.

    if (sortMerge) {
      using Details::shellSortKeysAndValues;
      offset_type new_nnz = 0;
      Kokkos::parallel_reduce
        ("Sort & merge rows of local transpose",
         range_type (0, numLocalCols),
         KOKKOS_LAMBDA (const LO lclCol, offset_type& nnz) {
          const offset_type beg = t_offsets[lclCol];
          const LO len (t_offsets[lclCol+1] - t_offsets[lclCol]);

          LO* cols_beg = t_cols.data () + beg;
          LO* cols_end = cols_beg + len;
          IST* vals_beg = t_vals.data () + beg;
          IST* vals_end = vals_beg + len;
          shellSortKeysAndValues (cols_beg, vals_beg, len);
          merge2 (cols_end, vals_end, cols_beg, cols_end, vals_beg, vals_end);
          const offset_type new_len = (cols_end - cols_beg);
          t_counts[lclCol] = new_len;
          nnz += new_len;
        }, new_nnz);

      TEUCHOS_ASSERT( new_nnz <= numLocalNnz );
      if (new_nnz != numLocalNnz) {

        // FIXME (mfh 11 Jul 2019) We don't have tests for this case
        // in Tpetra.
        TEUCHOS_ASSERT( false );

        index_type t_cols_packed
          (view_no_init ("transpose lcl ind packed"), new_nnz);
        values_type t_vals_packed
          (view_no_init ("transpose val packed"), new_nnz);
        row_map_type t_offsets_packed
          (view_no_init ("transpose ptr packed"), numLocalCols + 1);
        Kokkos::parallel_scan
          ("Recompute offsets & pack rows of local transpose",
           range_type (0, numLocalCols+1),
           KOKKOS_LAMBDA (const LO lclCol,
                          offset_type& running_count,
                          const bool final) {
            const offset_type new_count = lclCol < numLocalCols ?
              t_counts[lclCol] : offset_type (0);
            if (final) {
              t_offsets_packed[lclCol] = running_count;

              if (lclCol < numLocalCols) {
                const offset_type unpacked_beg = t_offsets[lclCol];
                const offset_type packed_beg = running_count;
                for (offset_type k = 0; k < new_count; ++k) {
                  t_cols_packed[packed_beg + k] = t_cols[unpacked_beg + k];
                  t_vals_packed[packed_beg + k] = t_vals[unpacked_beg + k];
                }
              }
            }
            running_count += new_count;
          });

        t_vals = t_vals_packed;
        t_cols = t_cols_packed;
        t_offsets = t_offsets_packed;
      }
    }

    local_matrix_type lclTransposeMatrix ("transpose", numLocalCols,
                                          numLocalRows, numLocalNnz,
                                          t_vals, t_offsets, t_cols);
    transMatrixWithSharedRows =
        rcp (new crs_matrix_type (lclTransposeMatrix,
                                  origMatrix_->getColMap (),
                                  origMatrix_->getRowMap (),
                                  origMatrix_->getRangeMap (),
                                  origMatrix_->getDomainMap ()));
  }
  else {
    // FIXME: are we ever here? There is no RowMatrix constructor, we seem to
    // be working only with CrsMatrices

    // Determine how many nonzeros there are per row in the transpose.
    Array<size_t> CurrentStart(numLocalCols,0);
    ArrayView<const LO> localIndices;
    ArrayView<const Scalar> localValues;
    // RowMatrix path
    for (LO i = 0; i < numLocalRows; ++i) {
      const size_t numEntriesInRow = origMatrix_->getNumEntriesInLocalRow(i);
      origMatrix_->getLocalRowView(i, localIndices, localValues);
      for (size_t j=0; j<numEntriesInRow; ++j) {
        ++CurrentStart[ localIndices[j] ];
      }
    }

    // create temporary row-major storage for the transposed matrix

    ArrayRCP<size_t> rowptr_rcp(numLocalCols+1);
    ArrayRCP<LO>     colind_rcp(numLocalNnz);
    ArrayRCP<Scalar> values_rcp(numLocalNnz);

    // Since ArrayRCP's are slow...
    ArrayView<size_t> TransRowptr = rowptr_rcp();
    ArrayView<LO>     TransColind = colind_rcp();
    ArrayView<Scalar> TransValues = values_rcp();

    // Scansum the TransRowptr; reset CurrentStart
    TransRowptr[0] = 0;
    for (LO i = 1; i < numLocalCols+1; ++i) {
      TransRowptr[i]  = CurrentStart[i-1] + TransRowptr[i-1];
    }
    for (LO i = 0; i < numLocalCols; ++i) {
      CurrentStart[i] = TransRowptr[i];
    }

    // populate the row-major storage so that the data for the transposed
    // matrix is easy to access
    for (LO i = 0; i < numLocalRows; ++i) {
      const size_t numEntriesInRow = origMatrix_->getNumEntriesInLocalRow (i);
      origMatrix_->getLocalRowView(i, localIndices, localValues);

      for (size_t j=0; j<numEntriesInRow; ++j) {
        size_t idx = CurrentStart[localIndices[j]];
        TransColind[idx] = Teuchos::as<LO>(i);
        TransValues[idx] = localValues[j];
        ++CurrentStart[localIndices[j]];
      }
    } //for (size_t i=0; i<numLocalRows; ++i)

    if (sortMerge) {
      // FIXME (mfh 11 Jul 2019) Just convert the input RowMatrix to
      // local_matrix_type and use the above routine.
      TEUCHOS_ASSERT( false ); // not implemented
    }

    // Allocate and populate temporary matrix with rows not uniquely owned
    transMatrixWithSharedRows =
        rcp (new crs_matrix_type (origMatrix_->getColMap (),
                                  origMatrix_->getRowMap (), 0));
    transMatrixWithSharedRows->setAllValues (rowptr_rcp, colind_rcp, values_rcp);

    // Call ESFC & return
    Teuchos::ParameterList eParams;
#ifdef HAVE_TPETRA_MMM_TIMINGS
    eParams.set("Timer Label",label_);
#endif
    if (! params.is_null ()) {
      const char paramName[] = "compute global constants";
      eParams.set(paramName, params->get(paramName, true));
    }

    using Teuchos::rcpFromRef;
    transMatrixWithSharedRows->expertStaticFillComplete (origMatrix_->getRangeMap (),
                                                         origMatrix_->getDomainMap (),
                                                         myImport, myExport,
                                                         rcpFromRef (eParams));
  }

  return transMatrixWithSharedRows;
}
//
// Explicit instantiation macro
//
// Must be expanded from within the Tpetra namespace!
//

#define TPETRA_ROWMATRIXTRANSPOSER_INSTANT(SCALAR,LO,GO,NODE) \
  template class RowMatrixTransposer< SCALAR, LO , GO , NODE >;

} // namespace Tpetra

#endif
