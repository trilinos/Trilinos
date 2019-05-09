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

#if defined( KOKKOS_ENABLE_CUDA )
  template< class D >
  struct device_is_cuda : public Kokkos::Impl::is_same<D,Kokkos::Cuda> {};
#else
  template< class D >
  struct device_is_cuda : public Kokkos::Impl::false_type {};
#endif

  // Functors for implementing packAndPrepare and unpackAndCombine
  // through parallel_for

  template <typename DS, typename ... DP,
            typename SS, typename ... SP,
            typename IdxView>
  struct PackArraySingleColumn<
    Kokkos::View<Sacado::MP::Vector<DS>*,DP...>,
    Kokkos::View<const Sacado::MP::Vector<SS>**,SP...>,
    IdxView >
  {
    typedef Kokkos::View<Sacado::MP::Vector<DS>*,DP...> DstView;
    typedef Kokkos::View<const Sacado::MP::Vector<SS>**,SP...> SrcView;
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
                     size_t col) {
      if ( Details::device_is_cuda<execution_space>::value )
        Kokkos::parallel_for(
          Kokkos::MPVectorWorkConfig<execution_space>( idx.size(), Kokkos::dimension_scalar(dst) ),
          PackArraySingleColumn(dst,src,idx,col) );
      else
         Kokkos::parallel_for( idx.size(),
                               PackArraySingleColumn(dst,src,idx,col) );
    }
  };

  template <typename DS, typename ... DP,
            typename SS, typename ... SP,
            typename IdxView>
  struct PackArrayMultiColumn<
    Kokkos::View<Sacado::MP::Vector<DS>*,DP...>,
    Kokkos::View<const Sacado::MP::Vector<SS>**,SP...>,
    IdxView >
  {
    typedef Kokkos::View<Sacado::MP::Vector<DS>*,DP...> DstView;
    typedef Kokkos::View<const Sacado::MP::Vector<SS>**,SP...> SrcView;
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
                     size_t numCols) {
      if ( Details::device_is_cuda<execution_space>::value )
        Kokkos::parallel_for(
          Kokkos::MPVectorWorkConfig<execution_space>( idx.size(), Kokkos::dimension_scalar(dst) ),
          PackArrayMultiColumn(dst,src,idx,numCols) );
      else
        Kokkos::parallel_for( idx.size(),
                              PackArrayMultiColumn(dst,src,idx,numCols) );
    }
  };

  template <typename DS, typename ... DP,
            typename SS, typename ... SP,
            typename IdxView,
            typename ColView>
  struct PackArrayMultiColumnVariableStride<
    Kokkos::View<Sacado::MP::Vector<DS>*,DP...>,
    Kokkos::View<const Sacado::MP::Vector<SS>**,SP...>,
    IdxView,
    ColView>
  {
    typedef Kokkos::View<Sacado::MP::Vector<DS>*,DP...> DstView;
    typedef Kokkos::View<const Sacado::MP::Vector<SS>**,SP...> SrcView;
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
                     size_t numCols) {
      if ( Details::device_is_cuda<execution_space>::value )
        Kokkos::parallel_for(
          Kokkos::MPVectorWorkConfig<execution_space>( idx.size(), Kokkos::dimension_scalar(dst) ),
          PackArrayMultiColumnVariableStride(dst,src,idx,col,numCols) );
      else
        Kokkos::parallel_for( idx.size(),
                              PackArrayMultiColumnVariableStride(
                                dst,src,idx,col,numCols) );
    }
  };

  template <typename ExecutionSpace,
            typename DS,
            typename ... DP,
            typename SS,
            typename ... SP,
            typename IdxView,
            typename Op>
  struct UnpackArrayMultiColumn<
    ExecutionSpace,
    Kokkos::View<Sacado::MP::Vector<DS>**,DP...>,
    Kokkos::View<const Sacado::MP::Vector<SS>*,SP...>,
    IdxView,
    Op >
  {
    typedef Kokkos::View<Sacado::MP::Vector<DS>**,DP...> DstView;
    typedef Kokkos::View<const Sacado::MP::Vector<SS>*,SP...> SrcView;
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

    KOKKOS_INLINE_FUNCTION void
    operator() (const size_type k) const
    {
      const typename IdxView::value_type localRow = idx(k);
      const size_t offset = k*numCols;
      for (size_t j = 0; j < numCols; ++j) {
        op (dst(localRow, j), src(offset+j));
      }
    }

    KOKKOS_INLINE_FUNCTION void
    operator() (const size_type k, const size_type tidx) const
    {
      const typename IdxView::value_type localRow = idx(k);
      const size_t offset = k*numCols;
      for (size_t j = 0; j < numCols; ++j) {
        op (dst(localRow, j).fastAccessCoeff(tidx),
            src(offset+j).fastAccessCoeff(tidx));
      }
    }

    static void
    unpack (const ExecutionSpace& execSpace,
            const DstView& dst,
            const SrcView& src,
            const IdxView& idx,
            const Op& op,
            const size_t numCols)
    {
      if (Details::device_is_cuda<execution_space>::value) {
        Kokkos::parallel_for
          ("Tpetra::MultiVector (Sacado::MP::Vector) unpack (constant stride)",
           Kokkos::MPVectorWorkConfig<execution_space> (idx.size (), Kokkos::dimension_scalar (dst)),
           UnpackArrayMultiColumn (execSpace, dst, src, idx, op, numCols));
      }
      else {
        Kokkos::parallel_for
          ("Tpetra::MultiVector (Sacado::MP::Vector) unpack (constant stride)",
           Kokkos::RangePolicy<execution_space, size_type> (0, idx.size ()),
           UnpackArrayMultiColumn (execSpace, dst, src, idx, op, numCols));
      }
    }
  };

  template <typename ExecutionSpace,
            typename DS,
            typename ... DP,
            typename SS,
            typename ... SP,
            typename IdxView,
            typename ColView,
            typename Op>
  struct UnpackArrayMultiColumnVariableStride<
    ExecutionSpace,
    Kokkos::View<Sacado::MP::Vector<DS>**,DP...>,
    Kokkos::View<const Sacado::MP::Vector<SS>*,SP...>,
    IdxView,
    ColView,
    Op>
  {
    typedef Kokkos::View<Sacado::MP::Vector<DS>**,DP...> DstView;
    typedef Kokkos::View<const Sacado::MP::Vector<SS>*,SP...> SrcView;
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

    KOKKOS_INLINE_FUNCTION void
    operator() (const size_type k) const
    {
      const typename IdxView::value_type localRow = idx(k);
      const size_t offset = k*numCols;
      for (size_t j = 0; j < numCols; ++j) {
        op (dst(localRow, col(j)), src(offset+j));
      }
    }

    KOKKOS_INLINE_FUNCTION void
    operator() (const size_type k, const size_type tidx) const
    {
      const typename IdxView::value_type localRow = idx(k);
      const size_t offset = k*numCols;
      for (size_t j = 0; j < numCols; ++j) {
        op (dst(localRow, col(j)).fastAccessCoeff(tidx),
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
            const size_t numCols)
    {
      if (Details::device_is_cuda<execution_space>::value) {
        Kokkos::parallel_for
          ("Tpetra::MultiVector unpack (Sacado::MP::Vector) (nonconstant stride)",
           Kokkos::MPVectorWorkConfig<execution_space> (idx.size (), Kokkos::dimension_scalar (dst)),
           UnpackArrayMultiColumnVariableStride (execSpace, dst, src, idx, col, op, numCols));
      }
      else {
        Kokkos::parallel_for
          ("Tpetra::MultiVector unpack (Sacado::MP::Vector) (nonconstant stride)",
           Kokkos::RangePolicy<execution_space, size_type> (0, idx.size ()),
           UnpackArrayMultiColumnVariableStride (execSpace, dst, src, idx, col, op, numCols));
      }
    }
  };

  template <typename DS, typename ... DP,
            typename SS, typename ... SP,
            typename DstIdxView, typename SrcIdxView>
  struct PermuteArrayMultiColumn<
    Kokkos::View<Sacado::MP::Vector<DS>**,DP...>,
    Kokkos::View<const Sacado::MP::Vector<SS>**,SP...>,
    DstIdxView,
    SrcIdxView>
  {
    typedef Kokkos::View<Sacado::MP::Vector<DS>**,DP...> DstView;
    typedef Kokkos::View<const Sacado::MP::Vector<SS>**,SP...> SrcView;
    typedef typename DstView::execution_space execution_space;
    typedef typename execution_space::size_type size_type;

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

    KOKKOS_INLINE_FUNCTION
    void operator()( const size_type k, const size_type tidx ) const {
      const typename DstIdxView::value_type toRow = dst_idx(k);
      const typename SrcIdxView::value_type fromRow = src_idx(k);
      for (size_t j = 0; j < numCols; ++j)
        dst(toRow, j).fastAccessCoeff(tidx) =
          src(fromRow, j).fastAccessCoeff(tidx);
    }

    static void permute(const DstView& dst,
                        const SrcView& src,
                        const DstIdxView& dst_idx,
                        const SrcIdxView& src_idx,
                        size_t numCols) {
      const size_type n = std::min( dst_idx.size(), src_idx.size() );
      if ( Details::device_is_cuda<execution_space>::value )
        Kokkos::parallel_for(
          Kokkos::MPVectorWorkConfig<execution_space>( n,  Kokkos::dimension_scalar(dst) ),
          PermuteArrayMultiColumn(dst,src,dst_idx,src_idx,numCols) );
      else
        Kokkos::parallel_for(
          n, PermuteArrayMultiColumn(dst,src,dst_idx,src_idx,numCols) );
    }
  };

  template <typename DS, typename ... DP,
            typename SS, typename ... SP,
            typename DstIdxView, typename SrcIdxView,
            typename DstColView, typename SrcColView>
  struct PermuteArrayMultiColumnVariableStride<
    Kokkos::View<Sacado::MP::Vector<DS>**,DP...>,
    Kokkos::View<const Sacado::MP::Vector<SS>**,SP...>,
    DstIdxView, SrcIdxView,
    DstColView, SrcColView >
  {
    typedef Kokkos::View<Sacado::MP::Vector<DS>**,DP...> DstView;
    typedef Kokkos::View<const Sacado::MP::Vector<SS>**,SP...> SrcView;
    typedef typename DstView::execution_space execution_space;
    typedef typename execution_space::size_type size_type;

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

    KOKKOS_INLINE_FUNCTION
    void operator()( const size_type k, const size_type tidx ) const {
      const typename DstIdxView::value_type toRow = dst_idx(k);
      const typename SrcIdxView::value_type fromRow = src_idx(k);
      for (size_t j = 0; j < numCols; ++j)
        dst(toRow, dst_col(j)).fastAccessCoeff(tidx) =
          src(fromRow, src_col(j)).fastAccessCoeff(tidx);
    }

    static void permute(const DstView& dst,
                        const SrcView& src,
                        const DstIdxView& dst_idx,
                        const SrcIdxView& src_idx,
                        const DstColView& dst_col,
                        const SrcColView& src_col,
                        size_t numCols) {
      const size_type n = std::min( dst_idx.size(), src_idx.size() );
      if ( Details::device_is_cuda<execution_space>::value )
        Kokkos::parallel_for(
          Kokkos::MPVectorWorkConfig<execution_space>( n, Kokkos::dimension_scalar(dst) ),
          PermuteArrayMultiColumnVariableStride(
            dst,src,dst_idx,src_idx,dst_col,src_col,numCols) );
      else
        Kokkos::parallel_for(
          n, PermuteArrayMultiColumnVariableStride(
            dst,src,dst_idx,src_idx,dst_col,src_col,numCols) );
    }
  };

} // Details namespace
} // KokkosRefactor namespace
} // Tpetra namespace

#endif // TPETRA_KOKKOS_REFACTOR_DETAILS_MULTI_VECTOR_DIST_OBJECT_KERNELS_MP_VECTOR_HPP
