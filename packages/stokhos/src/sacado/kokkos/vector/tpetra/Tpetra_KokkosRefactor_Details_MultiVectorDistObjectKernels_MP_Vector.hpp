// @HEADER
// *****************************************************************************
//          Tpetra: Templated Linear Algebra Services Package
//
// Copyright 2008 NTESS and the Tpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef TPETRA_KOKKOS_REFACTOR_DETAILS_MULTI_VECTOR_DIST_OBJECT_KERNELS_MP_VECTOR_HPP
#define TPETRA_KOKKOS_REFACTOR_DETAILS_MULTI_VECTOR_DIST_OBJECT_KERNELS_MP_VECTOR_HPP

//----------------------------------------------------------------------------
// Specializations of Tpetra::MultiVector pack/unpack kernels for MP::Vector
//----------------------------------------------------------------------------

#include "Tpetra_KokkosRefactor_Details_MultiVectorDistObjectKernels.hpp"
#include "Stokhos_Sacado_Kokkos_MP_Vector.hpp"
#include "Kokkos_Parallel_MP_Vector.hpp"

namespace Tpetra {
namespace KokkosRefactor {
namespace Details {

  // Functors for implementing packAndPrepare and unpackAndCombine
  // through parallel_for

  template <typename DstView, typename SrcView, typename IdxView>
  struct PackArraySingleColumn<
    DstView, SrcView, IdxView,
    typename std::enable_if< Kokkos::is_view_mp_vector<DstView>::value &&
                             Kokkos::is_view_mp_vector<SrcView>::value >::type >
  {
    typedef typename DstView::execution_space execution_space;
    typedef typename execution_space::size_type size_type;

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

    KOKKOS_INLINE_FUNCTION
    void operator()( const size_type k, const size_type tidx ) const {
      dst(k).fastAccessCoeff(tidx) = src(idx(k), col).fastAccessCoeff(tidx);
    }

    static void pack(const DstView& dst,
                     const SrcView& src,
                     const IdxView& idx,
                     size_t col,
                     const execution_space &space) {
      Kokkos::parallel_for(
        Kokkos::MPVectorWorkConfig<execution_space>(
          space, idx.size(), Kokkos::dimension_scalar(dst) ),
        PackArraySingleColumn(dst,src,idx,col) );
    }
  };

  template <typename DstView, typename SrcView, typename IdxView>
  struct PackArrayMultiColumn<
    DstView, SrcView, IdxView,
    typename std::enable_if< Kokkos::is_view_mp_vector<DstView>::value &&
                             Kokkos::is_view_mp_vector<SrcView>::value >::type >
  {
    typedef typename DstView::execution_space execution_space;
    typedef typename execution_space::size_type size_type;

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

    KOKKOS_INLINE_FUNCTION
    void operator()( const size_type k, const size_type tidx ) const {
      const typename IdxView::value_type localRow = idx(k);
      const size_t offset = k*numCols;
      for (size_t j = 0; j < numCols; ++j)
        dst(offset + j).fastAccessCoeff(tidx) =
          src(localRow, j).fastAccessCoeff(tidx);
    }

    static void pack(const DstView& dst,
                     const SrcView& src,
                     const IdxView& idx,
                     size_t numCols,
                     const execution_space &space) {
      Kokkos::parallel_for(
        Kokkos::MPVectorWorkConfig<execution_space>(space, idx.size(), Kokkos::dimension_scalar(dst) ),
        PackArrayMultiColumn(dst,src,idx,numCols) );
    }
  };

  template <typename DstView, typename SrcView, typename IdxView,
            typename ColView>
  struct PackArrayMultiColumnVariableStride<
    DstView, SrcView, IdxView, ColView,
    typename std::enable_if< Kokkos::is_view_mp_vector<DstView>::value &&
                             Kokkos::is_view_mp_vector<SrcView>::value >::type >
  {
    typedef typename DstView::execution_space execution_space;
    typedef typename execution_space::size_type size_type;

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

    KOKKOS_INLINE_FUNCTION
    void operator()( const size_type k, const size_type tidx ) const {
      const typename IdxView::value_type localRow = idx(k);
      const size_t offset = k*numCols;
      for (size_t j = 0; j < numCols; ++j)
        dst(offset + j).fastAccessCoeff(tidx) =
          src(localRow, col(j)).fastAccessCoeff(tidx);
    }

    static void pack(const DstView& dst,
                     const SrcView& src,
                     const IdxView& idx,
                     const ColView& col,
                     size_t numCols,
                     const execution_space &space) {
      Kokkos::parallel_for(
        Kokkos::MPVectorWorkConfig<execution_space>(
          space, idx.size(), Kokkos::dimension_scalar(dst) ),
        PackArrayMultiColumnVariableStride(dst,src,idx,col,numCols) );
    }
  };

  template <typename ExecutionSpace,
            typename DstView,
            typename SrcView,
            typename IdxView,
            typename Op>
  struct UnpackArrayMultiColumn<
    ExecutionSpace, DstView, SrcView, IdxView, Op,
    typename std::enable_if< Kokkos::is_view_mp_vector<DstView>::value &&
                             Kokkos::is_view_mp_vector<SrcView>::value >::type >
  {
    typedef typename ExecutionSpace::execution_space execution_space;
    typedef typename execution_space::size_type size_type;

    DstView dst;
    SrcView src;
    IdxView idx;
    Op op;
    size_t numCols;

    UnpackArrayMultiColumn (const ExecutionSpace& /* execSpace */,
                            const DstView& dst_,
                            const SrcView& src_,
                            const IdxView& idx_,
                            const Op& op_,
                            const size_t numCols_) :
      dst (dst_),
      src (src_),
      idx (idx_),
      op (op_),
      numCols (numCols_)
    {}

    template <class TagType>
    KOKKOS_INLINE_FUNCTION void
    operator() (TagType tag, const size_type k) const
    {
      const typename IdxView::value_type localRow = idx(k);
      const size_t offset = k*numCols;
      for (size_t j = 0; j < numCols; ++j) {
        op (tag, dst(localRow, j), src(offset+j));
      }
    }

    template <class TagType>
    KOKKOS_INLINE_FUNCTION void
    operator() (TagType tag, const size_type k, const size_type tidx) const
    {
      const typename IdxView::value_type localRow = idx(k);
      const size_t offset = k*numCols;
      for (size_t j = 0; j < numCols; ++j) {
        op (tag, dst(localRow, j).fastAccessCoeff(tidx),
            src(offset+j).fastAccessCoeff(tidx));
      }
    }

    static void
    unpack (const ExecutionSpace& execSpace,
            const DstView& dst,
            const SrcView& src,
            const IdxView& idx,
            const Op& op,
            const size_t numCols,
            const bool use_atomic_updates)
    {
      if (use_atomic_updates) {
        Kokkos::parallel_for
          ("Tpetra::MultiVector (Sacado::MP::Vector) unpack (constant stride)",
           Kokkos::MPVectorWorkConfig<execution_space, atomic_tag> (
             idx.size (), Kokkos::dimension_scalar (dst)),
           UnpackArrayMultiColumn (execSpace, dst, src, idx, op, numCols));
      }
      else {
        Kokkos::parallel_for
          ("Tpetra::MultiVector (Sacado::MP::Vector) unpack (constant stride)",
           Kokkos::MPVectorWorkConfig<execution_space, nonatomic_tag> (
             idx.size (), Kokkos::dimension_scalar (dst)),
           UnpackArrayMultiColumn (execSpace, dst, src, idx, op, numCols));
      }
    }
  };

  template <typename ExecutionSpace,
            typename DstView,
            typename SrcView,
            typename IdxView,
            typename ColView,
            typename Op>
  struct UnpackArrayMultiColumnVariableStride<
    ExecutionSpace, DstView, SrcView, IdxView, ColView, Op,
    typename std::enable_if< Kokkos::is_view_mp_vector<DstView>::value &&
                             Kokkos::is_view_mp_vector<SrcView>::value >::type >
  {
    typedef typename ExecutionSpace::execution_space execution_space;
    typedef typename execution_space::size_type size_type;

    DstView dst;
    SrcView src;
    IdxView idx;
    ColView col;
    Op op;
    size_t numCols;

    UnpackArrayMultiColumnVariableStride (const ExecutionSpace& /* execSpace */,
                                          const DstView& dst_,
                                          const SrcView& src_,
                                          const IdxView& idx_,
                                          const ColView& col_,
                                          const Op& op_,
                                          const size_t numCols_) :
      dst (dst_),
      src (src_),
      idx (idx_),
      col (col_),
      op (op_),
      numCols (numCols_)
    {}

    template <class TagType>
    KOKKOS_INLINE_FUNCTION void
    operator() (TagType tag, const size_type k) const
    {
      const typename IdxView::value_type localRow = idx(k);
      const size_t offset = k*numCols;
      for (size_t j = 0; j < numCols; ++j) {
        op (tag, dst(localRow, col(j)), src(offset+j));
      }
    }

    template <class TagType>
    KOKKOS_INLINE_FUNCTION void
    operator() (TagType tag, const size_type k, const size_type tidx) const
    {
      const typename IdxView::value_type localRow = idx(k);
      const size_t offset = k*numCols;
      for (size_t j = 0; j < numCols; ++j) {
        op (tag, dst(localRow, col(j)).fastAccessCoeff(tidx),
            src(offset+j).fastAccessCoeff(tidx));
      }
    }

    static void
    unpack (const ExecutionSpace& execSpace,
            const DstView& dst,
            const SrcView& src,
            const IdxView& idx,
            const ColView& col,
            const Op& op,
            const size_t numCols,
            const bool use_atomic_updates)
    {
      if (use_atomic_updates) {
        Kokkos::parallel_for
          ("Tpetra::MultiVector unpack (Sacado::MP::Vector) (nonconstant stride)",
           Kokkos::MPVectorWorkConfig<execution_space, atomic_tag> (
             idx.size (), Kokkos::dimension_scalar (dst)),
           UnpackArrayMultiColumnVariableStride (execSpace, dst, src, idx, col, op, numCols));
      }
      else {
        Kokkos::parallel_for
          ("Tpetra::MultiVector unpack (Sacado::MP::Vector) (nonconstant stride)",
           Kokkos::MPVectorWorkConfig<execution_space, nonatomic_tag> (
             idx.size (), Kokkos::dimension_scalar (dst)),
           UnpackArrayMultiColumnVariableStride (execSpace, dst, src, idx, col, op, numCols));
      }
    }
  };

  template <typename DstView, typename SrcView,
            typename DstIdxView, typename SrcIdxView, typename Op>
  struct PermuteArrayMultiColumn<
    DstView, SrcView, DstIdxView, SrcIdxView, Op,
    typename std::enable_if< Kokkos::is_view_mp_vector<DstView>::value &&
                             Kokkos::is_view_mp_vector<SrcView>::value >::type >
  {
    typedef typename DstView::execution_space execution_space;
    typedef typename execution_space::size_type size_type;

    DstView dst;
    SrcView src;
    DstIdxView dst_idx;
    SrcIdxView src_idx;
    size_t numCols;
    Op op;

    PermuteArrayMultiColumn(const DstView& dst_,
                            const SrcView& src_,
                            const DstIdxView& dst_idx_,
                            const SrcIdxView& src_idx_,
                            size_t numCols_, const Op& op_) :
      dst(dst_), src(src_), dst_idx(dst_idx_), src_idx(src_idx_),
      numCols(numCols_), op(op_) {}

    KOKKOS_INLINE_FUNCTION
    void operator()( const size_type k ) const {
      const typename DstIdxView::value_type toRow = dst_idx(k);
      const typename SrcIdxView::value_type fromRow = src_idx(k);
      nonatomic_tag tag;  // permute does not need atomics
      for (size_t j = 0; j < numCols; ++j)
        op(tag, dst(toRow, j),src(fromRow, j));
    }

    KOKKOS_INLINE_FUNCTION
    void operator()( const size_type k, const size_type tidx ) const {
      const typename DstIdxView::value_type toRow = dst_idx(k);
      const typename SrcIdxView::value_type fromRow = src_idx(k);
      nonatomic_tag tag;  // permute does not need atomics
      for (size_t j = 0; j < numCols; ++j)
        op(tag, dst(toRow, j).fastAccessCoeff(tidx),
           src(fromRow, j).fastAccessCoeff(tidx));
    }

    template <typename ExecutionSpace>
    static void permute(const ExecutionSpace &space,
                        const DstView& dst,
                        const SrcView& src,
                        const DstIdxView& dst_idx,
                        const SrcIdxView& src_idx,
                        size_t numCols,
                        const Op& op) {
      const size_type n = std::min( dst_idx.size(), src_idx.size() );
      Kokkos::parallel_for(
        Kokkos::MPVectorWorkConfig<execution_space>(
          space, n,  Kokkos::dimension_scalar(dst) ),
        PermuteArrayMultiColumn(dst,src,dst_idx,src_idx,numCols,op) );
    }
  };

  template <typename DstView, typename SrcView,
            typename DstIdxView, typename SrcIdxView,
            typename DstColView, typename SrcColView, typename Op>
  struct PermuteArrayMultiColumnVariableStride<
    DstView, SrcView, DstIdxView, SrcIdxView, DstColView, SrcColView, Op,
    typename std::enable_if< Kokkos::is_view_mp_vector<DstView>::value &&
                             Kokkos::is_view_mp_vector<SrcView>::value >::type >
  {
    typedef typename DstView::execution_space execution_space;
    typedef typename execution_space::size_type size_type;

    DstView dst;
    SrcView src;
    DstIdxView dst_idx;
    SrcIdxView src_idx;
    DstColView dst_col;
    SrcColView src_col;
    size_t numCols;
    Op op;

    PermuteArrayMultiColumnVariableStride(const DstView& dst_,
                                          const SrcView& src_,
                                          const DstIdxView& dst_idx_,
                                          const SrcIdxView& src_idx_,
                                          const DstColView& dst_col_,
                                          const SrcColView& src_col_,
                                          size_t numCols_,
                                          const Op& op_) :
      dst(dst_), src(src_), dst_idx(dst_idx_), src_idx(src_idx_),
      dst_col(dst_col_), src_col(src_col_),
      numCols(numCols_), op(op_) {}

    KOKKOS_INLINE_FUNCTION
    void operator()( const size_type k ) const {
      const typename DstIdxView::value_type toRow = dst_idx(k);
      const typename SrcIdxView::value_type fromRow = src_idx(k);
      nonatomic_tag tag;  // permute does not need atomics
      for (size_t j = 0; j < numCols; ++j)
        op(tag, dst(toRow, dst_col(j)),src(fromRow, src_col(j)));
    }

    KOKKOS_INLINE_FUNCTION
    void operator()( const size_type k, const size_type tidx ) const {
      const typename DstIdxView::value_type toRow = dst_idx(k);
      const typename SrcIdxView::value_type fromRow = src_idx(k);
      nonatomic_tag tag;  // permute does not need atomics
      for (size_t j = 0; j < numCols; ++j)
        op(tag, dst(toRow, dst_col(j)).fastAccessCoeff(tidx),
           src(fromRow, src_col(j)).fastAccessCoeff(tidx));
    }

    template <typename ExecutionSpace>
    static void permute(const ExecutionSpace &space,
                        const DstView& dst,
                        const SrcView& src,
                        const DstIdxView& dst_idx,
                        const SrcIdxView& src_idx,
                        const DstColView& dst_col,
                        const SrcColView& src_col,
                        size_t numCols,
                        const Op& op) {
      const size_type n = std::min( dst_idx.size(), src_idx.size() );
      Kokkos::parallel_for(
        Kokkos::MPVectorWorkConfig<execution_space>(
          space, n, Kokkos::dimension_scalar(dst) ),
        PermuteArrayMultiColumnVariableStride(
          dst,src,dst_idx,src_idx,dst_col,src_col,numCols,op) );
    }
  };

} // Details namespace
} // KokkosRefactor namespace
} // Tpetra namespace

#endif // TPETRA_KOKKOS_REFACTOR_DETAILS_MULTI_VECTOR_DIST_OBJECT_KERNELS_MP_VECTOR_HPP
