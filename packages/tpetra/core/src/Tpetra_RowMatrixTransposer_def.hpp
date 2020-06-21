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
#include "Tpetra_Details_computeOffsets.hpp"
#include "Tpetra_Details_shortSort.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_TimeMonitor.hpp"

namespace Tpetra {

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
    // This always sorts the local matrix after communication, so
    //   no need to set "sorted = false" in parameters.
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

  const bool sort = [&] () {
    constexpr bool sortDefault = true; // see #4607 discussion
    const char sortParamName[] = "sort";
    return params.get () == nullptr ? sortDefault :
      params->get (sortParamName, sortDefault);
  } ();

  const LO lclNumCols (origMatrix_->getNodeNumCols ());
  const LO lclNumRows (origMatrix_->getNodeNumRows ());
  const size_t nnz (origMatrix_->getNodeNumEntries ());

  RCP<const crs_matrix_type> crsMatrix =
    rcp_dynamic_cast<const crs_matrix_type> (origMatrix_);
  if (crsMatrix.is_null ()) {
    auto rowMap = origMatrix_->getRowMap ();
    if (rowMap->isOneToOne ()) {
      Teuchos::Array<size_t> numEntPerRow (lclNumRows);
      for (LO lclRow = 0; lclRow < lclNumRows; ++lclRow) {
        numEntPerRow[lclRow] = origMatrix_->getNumEntriesInLocalRow (lclRow);
      }
      auto colMap = origMatrix_->getColMap ();

      RCP<crs_matrix_type> crsMatrix_nc =
        rcp (new crs_matrix_type (rowMap, colMap, numEntPerRow ()));

      // When source & target Maps are same, Import just copies.
      import_type imp (rowMap, rowMap);
      crsMatrix_nc->doImport (*origMatrix_, imp, Tpetra::REPLACE);
      crsMatrix_nc->fillComplete (origMatrix_->getDomainMap (),
                                  origMatrix_->getRangeMap ());
      crsMatrix = crsMatrix_nc;
    }
    else {
      TEUCHOS_ASSERT( false ); // not implemented (it wasn't before)
    }
  }

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
  Kokkos::View<LO*, DT> t_counts ("transpose row counts", lclNumCols);
  using range_type = Kokkos::RangePolicy<LO, execution_space>;
  Kokkos::parallel_for
    ("Compute row counts of local transpose",
     range_type (0, lclNumRows),
     KOKKOS_LAMBDA (const LO row) {
      auto rowView = lclGraph.rowConst(row);
      const LO length  = rowView.length;

      for (LO colID = 0; colID < length; ++colID) {
        const LO col = rowView(colID);
        Kokkos::atomic_fetch_add (&t_counts[col], LO (1));
      }
    });

  using Kokkos::view_alloc;
  using Kokkos::WithoutInitializing;
  row_map_type t_offsets
    (view_alloc ("transpose ptr", WithoutInitializing), lclNumCols + 1);

  // TODO (mfh 10 Jul 2019): This returns the sum of all counts,
  // which could be useful for checking nnz.
  (void) Details::computeOffsetsFromCounts (t_offsets, t_counts);

  index_type t_cols
    (view_alloc ("transpose lcl ind", WithoutInitializing), nnz);
  values_type t_vals
    (view_alloc ("transpose val", WithoutInitializing), nnz);
  Kokkos::parallel_for
    ("Compute local transpose",
     range_type (0, lclNumRows),
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

  if (sort) {
    // NOTE (mfh 11 Jul 2019) Merging is unnecessary: above
    // parallel_for visits each row of the original matrix once, so
    // there can be no duplicate column indices in the transpose.
    using Details::shellSortKeysAndValues;
    Kokkos::parallel_for
      ("Sort rows of local transpose",
       range_type (0, lclNumCols),
       KOKKOS_LAMBDA (const LO lclCol) {
        const offset_type beg = t_offsets[lclCol];
        const LO len (t_offsets[lclCol+1] - t_offsets[lclCol]);

        LO* cols_beg = t_cols.data () + beg;
        IST* vals_beg = t_vals.data () + beg;
        shellSortKeysAndValues (cols_beg, vals_beg, len);
      });
  }

  local_matrix_type lclTransposeMatrix ("transpose", lclNumCols,
                                        lclNumRows, nnz,
                                        t_vals, t_offsets, t_cols);

  // Prebuild the importers and exporters the no-communication way,
  // flipping the importers and exporters around.
  const auto origExport = origMatrix_->getGraph ()->getExporter ();
  RCP<const import_type> myImport = origExport.is_null () ?
    Teuchos::null : rcp (new import_type (*origExport));
  const auto origImport = origMatrix_->getGraph ()->getImporter ();
  RCP<const export_type> myExport = origImport.is_null () ?
    Teuchos::null : rcp (new export_type (*origImport));

  RCP<Teuchos::ParameterList> graphParams = Teuchos::null;
  if(!sort) {
    graphParams = rcp(new Teuchos::ParameterList);
    graphParams->set("sorted", false);
  }

  return rcp (new crs_matrix_type (lclTransposeMatrix,
                                   origMatrix_->getColMap (),
                                   origMatrix_->getRowMap (),
                                   origMatrix_->getRangeMap (),
                                   origMatrix_->getDomainMap (),
                                   myImport, myExport, graphParams));
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
