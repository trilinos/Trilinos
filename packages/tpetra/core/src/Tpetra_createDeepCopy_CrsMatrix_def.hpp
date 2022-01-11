#ifndef TPETRA_CREATEDEEPCOPY_CRSMATRIX_DEF_HPP
#define TPETRA_CREATEDEEPCOPY_CRSMATRIX_DEF_HPP

#ifdef TPETRA_ENABLE_DEPRECATED_CODE 

#include "Tpetra_CrsMatrix.hpp"
#include "Tpetra_Map.hpp"
#include "Tpetra_RowMatrix.hpp"
#include "Tpetra_Details_localDeepCopyRowMatrix.hpp"
#include "Teuchos_Array.hpp"
#include "Teuchos_ArrayView.hpp"
#include <memory>


namespace Tpetra {

namespace { // (anonymous)

#ifdef TPETRA_ENABLE_DEPRECATED_CODE
template<class SC, class LO, class GO, class NT>
// This function is deprecated, but users don't call it directly; it is a 
// helper from createDeepCopy.  createDeepCopy is also deprecated.
// We silence TPETRA_DEPRECATED warnings here to prevent noise from
// compilation of createDeepCopy.
// TPETRA_DEPRECATED
typename CrsMatrix<SC, LO, GO, NT>::local_matrix_type
localDeepCopyFillCompleteCrsMatrix (const CrsMatrix<SC, LO, GO, NT>& A)
{
  using Kokkos::view_alloc;
  using Kokkos::WithoutInitializing;
  using crs_matrix_type = CrsMatrix<SC, LO, GO, NT>;
  using local_matrix_type =
    typename crs_matrix_type::local_matrix_type;
  local_matrix_type A_lcl = A.getLocalMatrixDevice ();

  using local_graph_device_type = typename crs_matrix_type::local_graph_device_type;
  using inds_type = typename local_graph_device_type::entries_type;
  inds_type ind (view_alloc ("ind", WithoutInitializing),
                 A_lcl.graph.entries.extent (0));
  Kokkos::deep_copy (ind, A_lcl.graph.entries);

  using offsets_type =
    typename local_graph_device_type::row_map_type::non_const_type;
  offsets_type ptr (view_alloc ("ptr", WithoutInitializing),
                    A_lcl.graph.row_map.extent (0));
  Kokkos::deep_copy (ptr, A_lcl.graph.row_map);

  using values_type = typename local_matrix_type::values_type;
  values_type val (view_alloc ("val", WithoutInitializing),
                   A_lcl.values.extent (0));
  Kokkos::deep_copy (val, A_lcl.values);

  local_graph_device_type lclGraph (ind, ptr);
  const size_t numCols = A.getColMap ()->getNodeNumElements ();
  return local_matrix_type (A.getObjectLabel (), numCols, val, lclGraph);
}

} // namespace // (anonymous)

template<class SC, class LO, class GO, class NT>
TPETRA_DEPRECATED
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
      crs_matrix_type (A.getRowMap (), A.getColMap (),
                       entPerRow_av) :
      crs_matrix_type (A.getRowMap (), entPerRow_av);

    const bool hasViews = A.supportsRowViews ();
    
    typename crs_matrix_type::nonconst_global_inds_host_view_type inputIndsBuf;
    typename crs_matrix_type::nonconst_values_host_view_type inputValsBuf;
    if (! hasViews) {
      Kokkos::resize(inputIndsBuf,maxNumEnt);
      Kokkos::resize(inputValsBuf,maxNumEnt);
    }

    const auto& rowMap = * (A.getRowMap ());
    for (LO lclRow = 0; lclRow < lclNumRows; ++lclRow) {
      const GO gblRow = rowMap.getGlobalElement (lclRow);
      if (hasViews) {
        typename crs_matrix_type::global_inds_host_view_type inputInds;
        typename crs_matrix_type::values_host_view_type inputVals;
        A.getGlobalRowView (gblRow, inputInds, inputVals);
        // BAD BAD BAD
        // we want a better way than reinterpret casting back and forth between scalar_type and
        // impl_scalar_type everywhere
        A_copy.insertGlobalValues (gblRow, inputInds.extent(0),
                                   reinterpret_cast<const typename crs_matrix_type::scalar_type*>(inputVals.data()),
                                   inputInds.data());
      }
      else {
        const size_t lclNumEnt = A.getNumEntriesInLocalRow (lclRow);
        TEUCHOS_ASSERT(lclNumEnt <= maxNumEnt);
        size_t numEnt = 0;
        A.getGlobalRowCopy (gblRow, inputIndsBuf, inputValsBuf, numEnt);
        A_copy.insertGlobalValues (gblRow, numEnt, 
                                   reinterpret_cast<const typename crs_matrix_type::scalar_type*>(inputValsBuf.data()),
                                   inputIndsBuf.data());

      }
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
#endif

} // namespace Tpetra

//
// Explicit instantiation macro
//
// Must be expanded from within the Tpetra namespace!
//

#define TPETRA_CREATEDEEPCOPY_CRSMATRIX_INSTANT(SC, LO, GO, NT) \
  template CrsMatrix< SC , LO , GO , NT > \
  createDeepCopy (const RowMatrix<SC, LO, GO, NT>& );

#endif // TPETRA_ENABLE_DEPRECATED_CODE

#endif // TPETRA_CREATEDEEPCOPY_CRSMATRIX_DEF_HPP

