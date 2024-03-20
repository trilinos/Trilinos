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

#if KOKKOSKERNELS_VERSION < 40299
#error "Tpetra::Details::MatrixApplyHelper can only be included with KokkosKernels 4.3+"
#endif

namespace Tpetra {
namespace Details {

// Class that works pretty much like std::unique_ptr, except it has
// operator= and copy-constructor that leave the new object null.
//
// unique_ptr would force CrsMatrix and BlockCrsMatrix to have deleted operator=,
// but this version doesn't.
template<typename T>
struct LazyUniquePtr
{
  LazyUniquePtr() : ptr(nullptr) {}
  LazyUniquePtr(T* ptr_) : ptr(ptr_) {}
  LazyUniquePtr(const LazyUniquePtr<T>&) {ptr = nullptr;}
  void operator=(const LazyUniquePtr<T>&) {ptr = nullptr;}
  ~LazyUniquePtr() {if(ptr) delete ptr;}
  T* operator->(){return ptr;}
  T* get() {return ptr;}
  void assign(T* ptr_)
  {
    reset();
    ptr = ptr_;
  }
  void reset()
  {
    if(ptr)
      delete ptr;
    ptr = nullptr;
  }
  T* ptr;
};

/// Helper for CrsMatrix::apply and BlockCrsMatrix::apply. Converts rowptrs
/// of the local matrix to int if it's not already int-valued, and overflow won't occur.
/// This enables TPLs to be used in KokkosSparse::spmv.
///
/// This functionality is temporary. When KokkosKernels 4.4 is released, the default offset (rowptrs element type)
/// in KK's ETI will switch to int, and Tpetra's matrix types will also switch to int. This will make the int rowptrs handling
/// here unnecessary, and CrsMatrix/BlockCrsMatrix can simply store the SPMVHandle themselves.
template<typename LocalMatrix, typename IntLocalMatrix, typename MultiVectorLocalView>
struct MatrixApplyHelper
{
  using Rowptrs = typename LocalMatrix::row_map_type;
  using IntRowptrs = typename IntLocalMatrix::row_map_type;
  using XVectorType = typename MultiVectorLocalView::const_type;
  using YVectorType = MultiVectorLocalView;
  using SPMVHandle = KokkosSparse::SPMVHandle<typename LocalMatrix::device_type, LocalMatrix, XVectorType, YVectorType>;
  using SPMVHandleInt = KokkosSparse::SPMVHandle<typename LocalMatrix::device_type, IntLocalMatrix, XVectorType, YVectorType>;

  MatrixApplyHelper(size_t nnz_, const typename LocalMatrix::row_map_type& rowptrs, KokkosSparse::SPMVAlgorithm algo = KokkosSparse::SPMV_DEFAULT)
    : nnz(nnz_), handle(algo), handle_int(algo)
  {
    if(shouldUseIntRowptrs())
      fillRowptrsInt(rowptrs);
    // otherwise, leave rowptrs_int empty
  }

  bool shouldUseIntRowptrs() const
  {
    // Use int-valued rowptrs if:
    // - number of entries doesn't overflow int
    // - LocalMatrix's rowptrs is not already int-typed
    return nnz <= size_t(INT_MAX) && !std::is_same_v<int, typename Rowptrs::non_const_value_type>;
  }

  // Given the local matrix, return a version with int-typed rowptrs
  IntLocalMatrix getIntRowptrMatrix(const LocalMatrix& A)
  {
    if constexpr(KokkosSparse::is_crs_matrix_v<LocalMatrix>)
    {
      return IntLocalMatrix("", A.numRows(), A.numCols(), A.nnz(),
          A.values, rowptrs_int, A.graph.entries);
    }
    else
    {
      return IntLocalMatrix("", A.numRows(), A.numCols(), A.nnz(),
          A.values, rowptrs_int, A.graph.entries, A.blockDim());
    }
  }

  void fillRowptrsInt(const Rowptrs& rowptrs)
  {
    // Populate rowptrs_int as a deep copy of rowptrs, converted to int.
    // This needs to be in its own public function because KOKKOS_LAMBDA cannot be used inside constructor or protected/private
    typename IntRowptrs::non_const_type rowptrs_temp(
        Kokkos::view_alloc(Kokkos::WithoutInitializing, "rowptrs_int"), rowptrs.extent(0));
    Kokkos::parallel_for(Kokkos::RangePolicy<typename LocalMatrix::execution_space>(0, rowptrs.extent(0)),
      KOKKOS_LAMBDA(int i)
      {
        rowptrs_temp(i) = int(rowptrs(i));
      });
    rowptrs_int = rowptrs_temp;
  }

  size_t nnz;
  // SPMVHandles are lazily initialized by actual spmv calls.
  // We declare both here because we don't know until runtime (based on nnz) which should be used.
  SPMVHandle handle;
  SPMVHandleInt handle_int;
  IntRowptrs rowptrs_int;
};
}
}

#endif // TPETRA_MATRIX_APPLY_HELPER_HPP

