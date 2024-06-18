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
// ************************************************************************
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

