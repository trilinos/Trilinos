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
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_TimeMonitor.hpp"

namespace Tpetra {

template<class Scalar,
     class LocalOrdinal,
     class GlobalOrdinal,
     class Node>
RowMatrixTransposer<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
RowMatrixTransposer (const Teuchos::RCP<const crs_matrix_type>& origMatrix,const std::string & label)
  : origMatrix_(origMatrix), label_(label) {}

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
  std::string prefix = std::string("Tpetra ")+ label_ + std::string(": ");
  using Teuchos::TimeMonitor;
  Teuchos::RCP<Teuchos::TimeMonitor> MM = Teuchos::rcp(new TimeMonitor(*TimeMonitor::getNewTimer(prefix + std::string("Transpose TAFC"))));
#endif

  // If transMatrixWithSharedRows has an exporter, that's what we
  // want.  If it doesn't, the rows aren't actually shared, and we're
  // done!
  RCP<const Export<LocalOrdinal,GlobalOrdinal,Node> > exporter =
    transMatrixWithSharedRows->getGraph ()->getExporter ();
  if (exporter.is_null ()) {
    return transMatrixWithSharedRows;
  }
  else {
    Teuchos::ParameterList labelList;
#ifdef HAVE_TPETRA_MMM_TIMINGS
    labelList.set("Timer Label",label_);
#endif
    if(!params.is_null()) labelList.set("compute global constants",params->get("compute global constants",true));
    // Use the Export object to do a fused Export and fillComplete.
    return exportAndFillCompleteCrsMatrix<crs_matrix_type> (transMatrixWithSharedRows, *exporter,Teuchos::null,Teuchos::null,Teuchos::rcp(&labelList,false));
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
createTransposeLocal (const Teuchos::RCP<Teuchos::ParameterList> &params)
{
  using Teuchos::Array;
  using Teuchos::ArrayRCP;
  using Teuchos::ArrayView;
  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::rcp_dynamic_cast;
  typedef LocalOrdinal LO;
  typedef GlobalOrdinal GO;
  typedef Tpetra::Import<LO, GO, Node> import_type;
  typedef Tpetra::Export<LO, GO, Node> export_type;

#ifdef HAVE_TPETRA_MMM_TIMINGS
  std::string prefix = std::string("Tpetra ")+ label_ + std::string(": ");
  using Teuchos::TimeMonitor;
  Teuchos::RCP<Teuchos::TimeMonitor> MM = Teuchos::rcp(new TimeMonitor(*TimeMonitor::getNewTimer(prefix + std::string("Transpose Local"))));
#endif

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

  size_t numLocalCols = origMatrix_->getNodeNumCols();
  size_t numLocalRows = origMatrix_->getNodeNumRows();
  size_t numLocalNnz  = origMatrix_->getNodeNumEntries();

  RCP<const crs_matrix_type> crsMatrix =
    rcp_dynamic_cast<const crs_matrix_type> (origMatrix_);

  RCP<crs_matrix_type> transMatrixWithSharedRows;
  if (crsMatrix != Teuchos::null) {
    // The matrix is a CrsMatrix
    using local_matrix_type = typename crs_matrix_type::local_matrix_type;
    using row_map_type = typename local_matrix_type::row_map_type::non_const_type;
    using index_type = typename local_matrix_type::index_type::non_const_type;
    using values_type = typename local_matrix_type::values_type::non_const_type;
    using execution_space = typename local_matrix_type::execution_space;

    // Determine how many nonzeros there are per row in the transpose.
    auto lclMatrix = crsMatrix->getLocalMatrix();
    auto lclGraph  = lclMatrix.graph;

    using range_type = Kokkos::RangePolicy<LO,execution_space>;

    // Determine how many nonzeros there are per row in the transpose.
    row_map_type t_rows("transpose_rows", numLocalCols+1);
    Kokkos::parallel_for("compute_number_of_indices_per_column", range_type(0, numLocalRows),
      KOKKOS_LAMBDA(const LO row) {
        auto rowView = lclGraph.rowConst(row);
        auto length  = rowView.length;

        for (decltype(length) colID = 0; colID < length; colID++) {
          auto col = rowView(colID);
          Kokkos::atomic_fetch_add(&t_rows[col], 1);
        }
    });

    // Compute offsets
    Kokkos::parallel_scan("compute_transpose_row_offsets", range_type(0, numLocalCols+1),
      KOKKOS_LAMBDA(const LO i, LO& update, const bool& final_pass) {
        const LO val = t_rows(i);
        if (final_pass)
          t_rows(i) = update;
        update += val;
      });

    row_map_type offsets("transpose_row_offsets_aux", numLocalCols+1);
    Kokkos::deep_copy(offsets, t_rows);

    index_type  t_cols("transpose_cols", numLocalNnz);
    values_type t_vals("transpose_vals", numLocalNnz);
    Kokkos::parallel_for("compute_transposed_rows", range_type(0, numLocalRows),
      KOKKOS_LAMBDA(const LO row) {
        auto rowView = lclMatrix.rowConst(row);
        auto length  = rowView.length;

        for (decltype(length) colID = 0; colID < length; colID++) {
          auto col = rowView.colidx(colID);

          LO insert_pos = Kokkos::atomic_fetch_add(&offsets[col], 1);

          t_cols[insert_pos] = row;
          t_vals[insert_pos] = rowView.value(colID);
        }
    });

    local_matrix_type lclTransposeMatrix("transpose", numLocalCols, numLocalRows, numLocalNnz, t_vals, t_rows, t_cols);

    transMatrixWithSharedRows =
        rcp (new crs_matrix_type (lclTransposeMatrix,
                                  origMatrix_->getColMap (), origMatrix_->getRowMap (),
                                  origMatrix_->getRangeMap (), origMatrix_->getDomainMap ()));
  } else {
    // FIXME: are we ever here? There is no RowMatrix constructor, we seem to
    // be working only with CrsMatrices

    // Determine how many nonzeros there are per row in the transpose.
    Array<size_t> CurrentStart(numLocalCols,0);
    ArrayView<const LO> localIndices;
    ArrayView<const Scalar> localValues;
    // RowMatrix path
    for (size_t i=0; i<numLocalRows; ++i) {
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
    TransRowptr[0]=0;
    for (size_t i=1; i<numLocalCols+1; ++i) TransRowptr[i]  = CurrentStart[i-1] + TransRowptr[i-1];
    for (size_t i=0; i<numLocalCols;   ++i) CurrentStart[i] = TransRowptr[i];

    // populate the row-major storage so that the data for the transposed
    // matrix is easy to access
    for (size_t i=0; i<numLocalRows; ++i) {
      const size_t numEntriesInRow = origMatrix_->getNumEntriesInLocalRow (i);
      origMatrix_->getLocalRowView(i, localIndices, localValues);

      for (size_t j=0; j<numEntriesInRow; ++j) {
        size_t idx = CurrentStart[localIndices[j]];
        TransColind[idx] = Teuchos::as<LO>(i);
        TransValues[idx] = localValues[j];
        ++CurrentStart[localIndices[j]];
      }
    } //for (size_t i=0; i<numLocalRows; ++i)

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
    if (!params.is_null())
      eParams.set("compute global constants", params->get("compute global constants", true));

    transMatrixWithSharedRows->expertStaticFillComplete (origMatrix_->getRangeMap (),
                                                         origMatrix_->getDomainMap (),
                                                         myImport, myExport,rcp(&eParams,false));

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
