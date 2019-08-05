#ifndef TPETRA_DETAILS_LOCALDEEPCOPYROWMATRIX_DEF_HPP
#define TPETRA_DETAILS_LOCALDEEPCOPYROWMATRIX_DEF_HPP

#include "Tpetra_Map.hpp"
#include "Tpetra_RowMatrix.hpp"
#include "Tpetra_Details_computeOffsets.hpp"
#include "Teuchos_TestForException.hpp"
#include "Teuchos_Array.hpp"
#include "Kokkos_Core.hpp"
#include <algorithm>
#include <stdexcept>

namespace Tpetra {
namespace Details {

template <class SC, class LO, class GO, class NT>
Impl::LocalRowOffsetsResult<NT>
localRowOffsets (const RowMatrix<SC, LO, GO, NT>& A)
{
  using Kokkos::view_alloc;
  using Kokkos::WithoutInitializing;
  using offsets_type =
    typename KokkosSparse::CrsMatrix<
      double, LO, typename NT::execution_space, void>::
        staticcrsgraph_type::row_map_type::non_const_type;
  using offset_type = typename offsets_type::non_const_value_type;

  const LO lclNumRows (A.getNodeNumRows ());
  offsets_type entPerRow;
  if (lclNumRows != 0) {
    entPerRow =
      offsets_type (view_alloc ("entPerRow", WithoutInitializing),
                    lclNumRows);
  }
  using host = Kokkos::DefaultHostExecutionSpace;
  auto entPerRow_h = Kokkos::create_mirror_view (host (), entPerRow);
  size_t maxNumEnt = 0;
  for (LO i = 0; i < lclNumRows; ++i) {
    const size_t lclNumEnt = A.getNumEntriesInLocalRow (i);
    entPerRow_h[i] = offset_type (lclNumEnt);
    maxNumEnt = maxNumEnt > lclNumEnt ? lclNumEnt : maxNumEnt;
  }
  Kokkos::deep_copy (entPerRow, entPerRow_h);

  offsets_type ptr;
  offset_type nnz = 0;
  if (lclNumRows != 0) {
    ptr = offsets_type (view_alloc ("ptr", WithoutInitializing),
                        lclNumRows + 1);
    using ::Tpetra::Details::computeOffsetsFromCounts;
    nnz = computeOffsetsFromCounts (ptr, entPerRow);
  }
  return {ptr, nnz, maxNumEnt};
}

template <class SC, class LO, class GO, class NT>
KokkosSparse::CrsMatrix<
  typename Kokkos::ArithTraits<SC>::val_type,
    LO,
    typename NT::execution_space,
    void>
localDeepCopyLocallyIndexedRowMatrix
(const RowMatrix<SC, LO, GO, NT>& A,
 const char label[])
{
  TEUCHOS_TEST_FOR_EXCEPTION
    (A.isGloballyIndexed (), std::logic_error,
     "Tpetra::Details::localDeepCopyLocallyIndexedRowMatrix: "
     "Input RowMatrix is globally indexed.");
  const LO lclNumRows (A.getNodeNumRows ());

  using lro_result_type = Impl::LocalRowOffsetsResult<NT>;
  using offsets_type = typename lro_result_type::offsets_type;
  using offset_type = typename lro_result_type::offset_type;

  offsets_type ptr;
  offset_type nnz = 0;
  size_t maxNumEnt = 0;
  {
    const lro_result_type result = localRowOffsets (A);
    ptr = result.ptr;
    nnz = result.nnz;
    maxNumEnt = result.maxNumEnt;
  }

  using Kokkos::view_alloc;
  using Kokkos::WithoutInitializing;
  using IST = typename Kokkos::ArithTraits<SC>::val_type;
  using local_matrix_type = KokkosSparse::CrsMatrix<
    IST, LO, typename NT::execution_space, void>;
  using local_graph_type =
    typename local_matrix_type::staticcrsgraph_type;
  using inds_type = typename local_graph_type::entries_type;
  inds_type ind (view_alloc ("ind", WithoutInitializing), nnz);
  auto ind_h = Kokkos::create_mirror_view (ind);

  using values_type = typename local_matrix_type::values_type;
  values_type val (view_alloc ("val", WithoutInitializing), nnz);
  auto val_h = Kokkos::create_mirror_view (val);

  const bool hasViews = A.supportsRowViews ();

  Teuchos::Array<LO> inputIndsBuf;
  Teuchos::Array<SC> inputValsBuf;
  if (! hasViews) {
    inputIndsBuf.resize (maxNumEnt);
    inputValsBuf.resize (maxNumEnt);
  }

  offset_type curPos = 0;
  for (LO lclRow = 0; lclRow < lclNumRows; ++lclRow) {
    Teuchos::ArrayView<const LO> inputInds_av;
    Teuchos::ArrayView<const SC> inputVals_av;
    size_t numEnt = 0;
    if (hasViews) {
      A.getLocalRowView (lclRow, inputInds_av,
                         inputVals_av);
      numEnt = static_cast<size_t> (inputInds_av.size ());
    }
    else {
      A.getLocalRowCopy (lclRow, inputIndsBuf (),
                         inputValsBuf (), numEnt);
      inputInds_av = inputIndsBuf.view (0, numEnt);
      inputVals_av = inputValsBuf.view (0, numEnt);
    }
    const IST* inVals =
      reinterpret_cast<const IST*> (inputVals_av.getRawPtr ());
    const LO* inInds = inputInds_av.getRawPtr ();
    std::copy (inInds, inInds + numEnt, ind_h.data () + curPos);
    std::copy (inVals, inVals + numEnt, val_h.data () + curPos);
    curPos += offset_type (numEnt);
  }
  Kokkos::deep_copy (ind, ind_h);
  Kokkos::deep_copy (val, val_h);

  local_graph_type lclGraph (ind, ptr);
  const size_t numCols = A.getColMap ()->getNodeNumElements ();
  return local_matrix_type (label, numCols, val, lclGraph);
}

} // namespace Details
} // namespace Tpetra

//
// Explicit instantiation macros
//
// Must be expanded from within the Tpetra namespace!
//
#define TPETRA_DETAILS_LOCALDEEPCOPYROWMATRIX_INSTANT(SC, LO, GO, NT) \
namespace Details { \
  template Impl::LocalRowOffsetsResult<NT> \
  localRowOffsets (const RowMatrix<SC, LO, GO, NT>& A); \
  \
  template KokkosSparse::CrsMatrix< \
    Kokkos::ArithTraits<SC>::val_type, \
    LO, NT::execution_space, void> \
  localDeepCopyLocallyIndexedRowMatrix<SC, LO, GO, NT> \
    (const RowMatrix<SC, LO, GO, NT>& A, \
     const char label[]); \
}

#endif // TPETRA_DETAILS_LOCALDEEPCOPYROWMATRIX_DEF_HPP
