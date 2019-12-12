#ifndef TPETRA_CREATEDEEPCOPY_CRSMATRIX_DEF_HPP
#define TPETRA_CREATEDEEPCOPY_CRSMATRIX_DEF_HPP

#include "Tpetra_CrsMatrix.hpp"
#include "Tpetra_Map.hpp"
#include "Tpetra_RowMatrix.hpp"
#include "Tpetra_Details_localDeepCopyRowMatrix.hpp"
#include "Teuchos_Array.hpp"
#include "Teuchos_ArrayView.hpp"
#include <memory>

namespace Tpetra {

namespace { // (anonymous)

template<class SC, class LO, class GO, class NT>
typename CrsMatrix<SC, LO, GO, NT>::local_matrix_type
localDeepCopyFillCompleteCrsMatrix (const CrsMatrix<SC, LO, GO, NT>& A)
{
  using Kokkos::view_alloc;
  using Kokkos::WithoutInitializing;
  using crs_matrix_type = CrsMatrix<SC, LO, GO, NT>;
  using local_matrix_type =
    typename crs_matrix_type::local_matrix_type;
  local_matrix_type A_lcl = A.getLocalMatrix ();

  using local_graph_type = typename crs_matrix_type::local_graph_type;
  using inds_type = typename local_graph_type::entries_type;
  inds_type ind (view_alloc ("ind", WithoutInitializing),
                 A_lcl.graph.entries.extent (0));
  Kokkos::deep_copy (ind, A_lcl.graph.entries);

  using offsets_type =
    typename local_graph_type::row_map_type::non_const_type;
  offsets_type ptr (view_alloc ("ptr", WithoutInitializing),
                    A_lcl.graph.row_map.extent (0));
  Kokkos::deep_copy (ptr, A_lcl.graph.row_map);

  using values_type = typename local_matrix_type::values_type;
  values_type val (view_alloc ("val", WithoutInitializing),
                   A_lcl.values.extent (0));
  Kokkos::deep_copy (val, A_lcl.values);

  local_graph_type lclGraph (ind, ptr);
  const size_t numCols = A.getColMap ()->getNodeNumElements ();
  return local_matrix_type (A.getObjectLabel (), numCols, val, lclGraph);
}

} // namespace // (anonymous)

template<class SC, class LO, class GO, class NT>
CrsMatrix<SC, LO, GO, NT>
createDeepCopy (const RowMatrix<SC, LO, GO, NT>& A)
{
  using crs_matrix_type = CrsMatrix<SC, LO, GO, NT>;
  const crs_matrix_type* A_crs =
    dynamic_cast<const crs_matrix_type*> (&A);

  if (A_crs != nullptr && A_crs->isFillComplete ()) {
    auto A_lcl = localDeepCopyFillCompleteCrsMatrix (*A_crs);
    auto G = A_crs->getCrsGraph ();
    return crs_matrix_type (A_lcl, A_crs->getRowMap (),
                            A_crs->getColMap (),
                            A_crs->getDomainMap (),
                            A_crs->getRangeMap (),
                            G->getImporter (),
                            G->getExporter ());
  }
  else if (A.isGloballyIndexed ()) {
    const LO lclNumRows (A.getNodeNumRows ());

    std::unique_ptr<size_t[]> entPerRow (new size_t [lclNumRows]);
    size_t maxNumEnt = 0;
    for (LO i = 0; i < lclNumRows; ++i) {
      const size_t lclNumEnt = A.getNumEntriesInLocalRow (i);
      entPerRow[i] = lclNumEnt;
      maxNumEnt = maxNumEnt < lclNumEnt ? lclNumEnt : maxNumEnt;
    }

    Teuchos::ArrayView<const size_t> entPerRow_av
      (entPerRow.get (), lclNumRows);

    const bool hasColMap =
      A.hasColMap () && ! A.getColMap ().is_null ();

    crs_matrix_type A_copy = hasColMap ?
#ifdef TPETRA_ENABLE_DEPRECATED_CODE
      crs_matrix_type (A.getRowMap (), A.getColMap (),
                       entPerRow_av,
                       Tpetra::StaticProfile) :
      crs_matrix_type (A.getRowMap (), entPerRow_av,
                       Tpetra::StaticProfile);
#else // TPETRA_ENABLE_DEPRECATED_CODE
      crs_matrix_type (A.getRowMap (), A.getColMap (),
                       entPerRow_av) :
      crs_matrix_type (A.getRowMap (), entPerRow_av);
#endif // TPETRA_ENABLE_DEPRECATED_CODE

    const bool hasViews = A.supportsRowViews ();

    Teuchos::Array<GO> inputIndsBuf;
    Teuchos::Array<SC> inputValsBuf;
    if (! hasViews) {
      inputIndsBuf.resize (maxNumEnt);
      inputValsBuf.resize (maxNumEnt);
    }

    const auto& rowMap = * (A.getRowMap ());
    for (LO lclRow = 0; lclRow < lclNumRows; ++lclRow) {
      const GO gblRow = rowMap.getGlobalElement (lclRow);
      Teuchos::ArrayView<const GO> inputInds_av;
      Teuchos::ArrayView<const SC> inputVals_av;
      size_t numEnt = 0;
      if (hasViews) {
        A.getGlobalRowView (gblRow, inputInds_av, inputVals_av);
        numEnt = static_cast<size_t> (inputInds_av.size ());
      }
      else {
        const size_t lclNumEnt = A.getNumEntriesInLocalRow (lclRow);
        TEUCHOS_ASSERT(lclNumEnt <= maxNumEnt);
        A.getGlobalRowCopy (gblRow, inputIndsBuf (),
                            inputValsBuf (), numEnt);
        inputInds_av = inputIndsBuf.view (0, numEnt);
        inputVals_av = inputValsBuf.view (0, numEnt);
      }
      A_copy.insertGlobalValues (gblRow, inputInds_av, inputVals_av);
    }

    if (A.isFillComplete ()) {
      A_copy.fillComplete (A.getDomainMap (), A.getRangeMap ());
    }
    return A_copy;
  }
  else { // locally indexed or empty
    using Details::localDeepCopyLocallyIndexedRowMatrix;
    auto A_lcl = localDeepCopyLocallyIndexedRowMatrix (A, "A");

    Teuchos::RCP<const Export<LO, GO, NT>> exp;
    Teuchos::RCP<const Import<LO, GO, NT>> imp;
    auto G = A.getGraph ();
    if (! G.is_null ()) {
      imp = G->getImporter ();
      exp = G->getExporter ();
    }
    if (! imp.is_null () || ! exp.is_null ()) {
      return crs_matrix_type (A_lcl, A.getRowMap (),
                              A.getColMap (),
                              A.getDomainMap (),
                              A.getRangeMap (),
                              imp, exp);
    }
    else {
      return crs_matrix_type (A_lcl, A.getRowMap (),
                              A.getColMap (),
                              A.getDomainMap (),
                              A.getRangeMap ());
    }
  }
}

} // namespace Tpetra

//
// Explicit instantiation macro
//
// Must be expanded from within the Tpetra namespace!
//

#define TPETRA_CREATEDEEPCOPY_CRSMATRIX_INSTANT(SC, LO, GO, NT) \
  template CrsMatrix< SC , LO , GO , NT > \
  createDeepCopy (const RowMatrix<SC, LO, GO, NT>& );

#endif // TPETRA_CREATEDEEPCOPY_CRSMATRIX_DEF_HPP
