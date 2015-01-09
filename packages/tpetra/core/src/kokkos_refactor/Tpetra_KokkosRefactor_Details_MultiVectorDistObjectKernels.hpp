/*
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
*/

// mfh 13/14 Sep 2013 The "should use as<size_t>" comments are both
// incorrect (as() is not a device function) and usually irrelevant
// (it would only matter if LocalOrdinal were bigger than size_t on a
// particular platform, which is unlikely).

#ifndef TPETRA_KOKKOS_REFACTOR_DETAILS_MULTI_VECTOR_DIST_OBJECT_KERNELS_HPP
#define TPETRA_KOKKOS_REFACTOR_DETAILS_MULTI_VECTOR_DIST_OBJECT_KERNELS_HPP

#include "Kokkos_Core.hpp"
#include "Kokkos_ArithTraits.hpp"

namespace Tpetra {
namespace KokkosRefactor {
namespace Details {

  // Functors for implementing packAndPrepare and unpackAndCombine
  // through parallel_for

  template <typename DstView, typename SrcView, typename IdxView>
  struct PackArraySingleColumn {
    typedef typename DstView::device_type device_type;
    typedef typename device_type::size_type size_type;

    DstView dst;
    SrcView src;
    IdxView idx;
    size_t col;

    PackArraySingleColumn(const DstView& dst_,
                          const SrcView& src_,
                          const IdxView& idx_,
                          size_t col_) :
      dst(dst_), src(src_), idx(idx_), col(col_) {}

    KOKKOS_INLINE_FUNCTION
    void operator()( const size_type k ) const {
      dst(k) = src(idx(k), col);
    }

    static void pack(const DstView& dst,
                     const SrcView& src,
                     const IdxView& idx,
                     size_t col) {
      Kokkos::parallel_for( idx.size(),
                            PackArraySingleColumn(dst,src,idx,col) );
    }
  };

  // To do:  Add enable_if<> restrictions on DstView::Rank == 1,
  // SrcView::Rank == 2
  template <typename DstView, typename SrcView, typename IdxView>
  void pack_array_single_column(const DstView& dst,
                                const SrcView& src,
                                const IdxView& idx,
                                size_t col) {
    PackArraySingleColumn<DstView,SrcView,IdxView>::pack(
      dst, src, idx, col);
  }

  template <typename DstView, typename SrcView, typename IdxView>
  struct PackArrayMultiColumn {
    typedef typename DstView::device_type device_type;
    typedef typename device_type::size_type size_type;

    DstView dst;
    SrcView src;
    IdxView idx;
    size_t numCols;

    PackArrayMultiColumn(const DstView& dst_,
                         const SrcView& src_,
                         const IdxView& idx_,
                         size_t numCols_) :
      dst(dst_), src(src_), idx(idx_), numCols(numCols_) {}

    KOKKOS_INLINE_FUNCTION
    void operator()( const size_type k ) const {
      const typename IdxView::value_type localRow = idx(k);
      const size_t offset = k*numCols;
      for (size_t j = 0; j < numCols; ++j)
        dst(offset + j) = src(localRow, j);
    }

    static void pack(const DstView& dst,
                     const SrcView& src,
                     const IdxView& idx,
                     size_t numCols) {
      Kokkos::parallel_for( idx.size(),
                            PackArrayMultiColumn(dst,src,idx,numCols) );
    }
  };

  // To do:  Add enable_if<> restrictions on DstView::Rank == 1,
  // SrcView::Rank == 2
  template <typename DstView, typename SrcView, typename IdxView>
  void pack_array_multi_column(const DstView& dst,
                               const SrcView& src,
                               const IdxView& idx,
                               size_t numCols) {
    PackArrayMultiColumn<DstView,SrcView,IdxView>::pack(
      dst, src, idx, numCols);
  }

  template <typename DstView, typename SrcView, typename IdxView,
            typename ColView>
  struct PackArrayMultiColumnVariableStride {
    typedef typename DstView::device_type device_type;
    typedef typename device_type::size_type size_type;

    DstView dst;
    SrcView src;
    IdxView idx;
    ColView col;
    size_t numCols;

    PackArrayMultiColumnVariableStride(const DstView& dst_,
                                       const SrcView& src_,
                                       const IdxView& idx_,
                                       const ColView& col_,
                                       size_t numCols_) :
      dst(dst_), src(src_), idx(idx_), col(col_), numCols(numCols_) {}

    KOKKOS_INLINE_FUNCTION
    void operator()( const size_type k ) const {
      const typename IdxView::value_type localRow = idx(k);
      const size_t offset = k*numCols;
      for (size_t j = 0; j < numCols; ++j)
        dst(offset + j) = src(localRow, col(j));
    }

    static void pack(const DstView& dst,
                     const SrcView& src,
                     const IdxView& idx,
                     const ColView& col,
                     size_t numCols) {
      Kokkos::parallel_for( idx.size(),
                            PackArrayMultiColumnVariableStride(
                              dst,src,idx,col,numCols) );
    }
  };

  // To do:  Add enable_if<> restrictions on DstView::Rank == 1,
  // SrcView::Rank == 2
  template <typename DstView, typename SrcView, typename IdxView,
            typename ColView>
  void pack_array_multi_column_variable_stride(const DstView& dst,
                                               const SrcView& src,
                                               const IdxView& idx,
                                               const ColView& col,
                                               size_t numCols) {
    PackArrayMultiColumnVariableStride<DstView,SrcView,IdxView,ColView>::pack(
      dst, src, idx, col, numCols);
  }

  struct InsertOp {
    template <typename Scalar>
    KOKKOS_INLINE_FUNCTION
    void operator() (Scalar& dest, const Scalar& src) const {
      Kokkos::atomic_assign(&dest, src);
    }
  };
  struct AddOp {
    template <typename Scalar>
    KOKKOS_INLINE_FUNCTION
    void operator() (Scalar& dest, const Scalar& src) const {
      Kokkos::atomic_add(&dest, src);
    }
  };
  struct AbsMaxOp {
    // ETP:  Is this really what we want?  This seems very odd if
    // Scalar != SCT::mag_type (e.g., Scalar == std::complex<T>)
    template <typename T>
    KOKKOS_INLINE_FUNCTION
    T max(const T& a, const T& b) const { return a > b ? a : b; }

    template <typename Scalar>
    KOKKOS_INLINE_FUNCTION
    void operator() (Scalar& dest, const Scalar& src) const {
      typedef Kokkos::Details::ArithTraits<Scalar> SCT;
      Kokkos::atomic_assign(&dest, Scalar(max(SCT::abs(dest),SCT::abs(src))));
    }
  };

  template <typename DstView, typename SrcView, typename IdxView, typename Op>
  struct UnpackArrayMultiColumn {
    typedef typename DstView::device_type device_type;
    typedef typename device_type::size_type size_type;

    DstView dst;
    SrcView src;
    IdxView idx;
    Op op;
    size_t numCols;

    UnpackArrayMultiColumn(const DstView& dst_,
                           const SrcView& src_,
                           const IdxView& idx_,
                           const Op& op_,
                           size_t numCols_) :
      dst(dst_), src(src_), idx(idx_), op(op_), numCols(numCols_) {}

    KOKKOS_INLINE_FUNCTION
    void operator()( const size_type k ) const {
      const typename IdxView::value_type localRow = idx(k);
      const size_t offset = k*numCols;
      for (size_t j = 0; j < numCols; ++j)
        op( dst(localRow,j), src(offset+j) );
    }

    static void unpack(const DstView& dst,
                       const SrcView& src,
                       const IdxView& idx,
                       const Op& op,
                       size_t numCols) {
      Kokkos::parallel_for( idx.size(),
                            UnpackArrayMultiColumn(dst,src,idx,op,numCols) );
    }
  };

  // To do:  Add enable_if<> restrictions on DstView::Rank == 2,
  // SrcView::Rank == 1
  template <typename DstView, typename SrcView, typename IdxView, typename Op>
  void unpack_array_multi_column(const DstView& dst,
                                 const SrcView& src,
                                 const IdxView& idx,
                                 const Op& op,
                                 size_t numCols) {
    UnpackArrayMultiColumn<DstView,SrcView,IdxView,Op>::unpack(
      dst, src, idx, op, numCols);
  }

  template <typename DstView, typename SrcView, typename IdxView,
            typename ColView, typename Op>
  struct UnpackArrayMultiColumnVariableStride {
    typedef typename DstView::device_type device_type;
    typedef typename device_type::size_type size_type;

    DstView dst;
    SrcView src;
    IdxView idx;
    ColView col;
    Op op;
    size_t numCols;

    UnpackArrayMultiColumnVariableStride(const DstView& dst_,
                                         const SrcView& src_,
                                         const IdxView& idx_,
                                         const ColView& col_,
                                         const Op& op_,
                                         size_t numCols_) :
      dst(dst_), src(src_), idx(idx_), col(col_), op(op_), numCols(numCols_) {}

    KOKKOS_INLINE_FUNCTION
    void operator()( const size_type k ) const {
      const typename IdxView::value_type localRow = idx(k);
      const size_t offset = k*numCols;
      for (size_t j = 0; j < numCols; ++j)
        op( dst(localRow,col(j)), src(offset+j) );
    }

    static void unpack(const DstView& dst,
                       const SrcView& src,
                       const IdxView& idx,
                       const ColView& col,
                       const Op& op,
                       size_t numCols) {
      Kokkos::parallel_for( idx.size(),
                            UnpackArrayMultiColumnVariableStride(
                              dst,src,idx,col,op,numCols) );
    }
  };

  // To do:  Add enable_if<> restrictions on DstView::Rank == 2,
  // SrcView::Rank == 1
  template <typename DstView, typename SrcView,typename IdxView,
            typename ColView, typename Op>
  void unpack_array_multi_column_variable_stride(const DstView& dst,
                                                 const SrcView& src,
                                                 const IdxView& idx,
                                                 const ColView& col,
                                                 const Op& op,
                                                 size_t numCols) {
    UnpackArrayMultiColumnVariableStride<DstView,SrcView,IdxView,ColView,Op>::unpack(
      dst, src, idx, col, op, numCols);
  }

  template <typename DstView, typename SrcView,
            typename DstIdxView, typename SrcIdxView>
  struct PermuteArrayMultiColumn {
    typedef typename DstView::device_type device_type;
    typedef typename device_type::size_type size_type;

    DstView dst;
    SrcView src;
    DstIdxView dst_idx;
    SrcIdxView src_idx;
    size_t numCols;

    PermuteArrayMultiColumn(const DstView& dst_,
                            const SrcView& src_,
                            const DstIdxView& dst_idx_,
                            const SrcIdxView& src_idx_,
                            size_t numCols_) :
      dst(dst_), src(src_), dst_idx(dst_idx_), src_idx(src_idx_),
      numCols(numCols_) {}

    KOKKOS_INLINE_FUNCTION
    void operator()( const size_type k ) const {
      const typename DstIdxView::value_type toRow = dst_idx(k);
      const typename SrcIdxView::value_type fromRow = src_idx(k);
      for (size_t j = 0; j < numCols; ++j)
        dst(toRow, j) = src(fromRow, j);
    }

    static void permute(const DstView& dst,
                        const SrcView& src,
                        const DstIdxView& dst_idx,
                        const SrcIdxView& src_idx,
                        size_t numCols) {
      const size_type n = std::min( dst_idx.size(), src_idx.size() );
      Kokkos::parallel_for(
        n, PermuteArrayMultiColumn(dst,src,dst_idx,src_idx,numCols) );
    }
  };

  // To do:  Add enable_if<> restrictions on DstView::Rank == 1,
  // SrcView::Rank == 2
  template <typename DstView, typename SrcView,
            typename DstIdxView, typename SrcIdxView>
  void permute_array_multi_column(const DstView& dst,
                                  const SrcView& src,
                                  const DstIdxView& dst_idx,
                                  const SrcIdxView& src_idx,
                                  size_t numCols) {
    PermuteArrayMultiColumn<DstView,SrcView,DstIdxView,SrcIdxView>::permute(
      dst, src, dst_idx, src_idx, numCols);
  }

  template <typename DstView, typename SrcView,
            typename DstIdxView, typename SrcIdxView,
            typename DstColView, typename SrcColView>
  struct PermuteArrayMultiColumnVariableStride {
    typedef typename DstView::device_type device_type;
    typedef typename device_type::size_type size_type;

    DstView dst;
    SrcView src;
    DstIdxView dst_idx;
    SrcIdxView src_idx;
    DstColView dst_col;
    SrcColView src_col;
    size_t numCols;

    PermuteArrayMultiColumnVariableStride(const DstView& dst_,
                                          const SrcView& src_,
                                          const DstIdxView& dst_idx_,
                                          const SrcIdxView& src_idx_,
                                          const DstColView& dst_col_,
                                          const SrcColView& src_col_,
                                          size_t numCols_) :
      dst(dst_), src(src_), dst_idx(dst_idx_), src_idx(src_idx_),
      dst_col(dst_col_), src_col(src_col_),
      numCols(numCols_) {}

    KOKKOS_INLINE_FUNCTION
    void operator()( const size_type k ) const {
      const typename DstIdxView::value_type toRow = dst_idx(k);
      const typename SrcIdxView::value_type fromRow = src_idx(k);
      for (size_t j = 0; j < numCols; ++j)
        dst(toRow, dst_col(j)) = src(fromRow, src_col(j));
    }

    static void permute(const DstView& dst,
                        const SrcView& src,
                        const DstIdxView& dst_idx,
                        const SrcIdxView& src_idx,
                        const DstColView& dst_col,
                        const SrcColView& src_col,
                        size_t numCols) {
      const size_type n = std::min( dst_idx.size(), src_idx.size() );
      Kokkos::parallel_for(
        n, PermuteArrayMultiColumnVariableStride(
          dst,src,dst_idx,src_idx,dst_col,src_col,numCols) );
    }
  };

  // To do:  Add enable_if<> restrictions on DstView::Rank == 1,
  // SrcView::Rank == 2
  template <typename DstView, typename SrcView,
            typename DstIdxView, typename SrcIdxView,
            typename DstColView, typename SrcColView>
  void permute_array_multi_column_variable_stride(const DstView& dst,
                                                  const SrcView& src,
                                                  const DstIdxView& dst_idx,
                                                  const SrcIdxView& src_idx,
                                                  const DstColView& dst_col,
                                                  const SrcColView& src_col,
                                                  size_t numCols) {
    PermuteArrayMultiColumnVariableStride<DstView,SrcView,
      DstIdxView,SrcIdxView,DstColView,SrcColView>::permute(
      dst, src, dst_idx, src_idx, dst_col, src_col, numCols);
  }

} // Details namespace
} // KokkosRefactor namespace
} // Tpetra namespace

#endif // TPETRA_KOKKOS_REFACTOR_DETAILS_MULTI_VECTOR_DIST_OBJECT_KERNELS_HPP
