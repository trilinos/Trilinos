// @HEADER
// *****************************************************************************
//          Tpetra: Templated Linear Algebra Services Package
//
// Copyright 2008 NTESS and the Tpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef TPETRA_MATRIX_APPLY_HELPER_HPP
#define TPETRA_MATRIX_APPLY_HELPER_HPP

#include "Kokkos_Core.hpp"
#include "KokkosSparse_spmv_handle.hpp"
#include "Tpetra_Details_IntRowPtrHelper.hpp"

namespace Tpetra {
namespace Details {

/// Helper for CrsMatrix::apply and BlockCrsMatrix::apply. Converts rowptrs
/// of the local matrix to int if it's not already int-valued, and overflow won't occur.
/// This enables TPLs to be used in KokkosSparse::spmv.
///
/// This functionality is temporary. When KokkosKernels 4.4 is released, the default offset (rowptrs element type)
/// in KK's ETI will switch to int, and Tpetra's matrix types will also switch to int. This will make the int rowptrs handling
/// here unnecessary, and CrsMatrix/BlockCrsMatrix can simply store the SPMVHandle themselves.
template<typename LocalMatrix, typename IntLocalMatrix, typename MultiVectorLocalView>
struct MatrixApplyHelper : public IntRowPtrHelper<LocalMatrix, IntLocalMatrix>
{
  // using Rowptrs = typename LocalMatrix::row_map_type;
  // using IntRowptrs = typename IntLocalMatrix::row_map_type;
  using XVectorType = typename MultiVectorLocalView::const_type;
  using YVectorType = MultiVectorLocalView;
  using SPMVHandle = KokkosSparse::SPMVHandle<typename LocalMatrix::device_type, LocalMatrix, XVectorType, YVectorType>;
  using SPMVHandleInt = KokkosSparse::SPMVHandle<typename LocalMatrix::device_type, IntLocalMatrix, XVectorType, YVectorType>;

  MatrixApplyHelper(size_t nnz, const typename LocalMatrix::row_map_type& rowptrs, KokkosSparse::SPMVAlgorithm algo = KokkosSparse::SPMV_DEFAULT)
    : IntRowPtrHelper<LocalMatrix, IntLocalMatrix>(nnz, rowptrs), handle_int(algo) {}
  

  // SPMVHandles are lazily initialized by actual spmv calls.
  // We declare both here because we don't know until runtime (based on nnz) which should be used.
  SPMVHandle handle;
  SPMVHandleInt handle_int;
};
}
}

#endif // TPETRA_MATRIX_APPLY_HELPER_HPP

