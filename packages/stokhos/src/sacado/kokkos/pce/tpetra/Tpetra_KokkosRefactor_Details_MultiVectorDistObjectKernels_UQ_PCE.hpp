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

#ifndef TPETRA_KOKKOS_REFACTOR_DETAILS_MULTI_VECTOR_DIST_OBJECT_KERNELS_UQ_PCE_HPP
#define TPETRA_KOKKOS_REFACTOR_DETAILS_MULTI_VECTOR_DIST_OBJECT_KERNELS_UQ_PCE_HPP

//----------------------------------------------------------------------------
// Specializations of Tpetra::MultiVector pack/unpack kernels for UQ::PCE
//----------------------------------------------------------------------------

#include "Tpetra_KokkosRefactor_Details_MultiVectorDistObjectKernels.hpp"
#include "Stokhos_Sacado_Kokkos_UQ_PCE.hpp"
#include "Kokkos_Parallel_MP_Vector.hpp"

// For device_is_cuda<>
#include "Tpetra_KokkosRefactor_Details_MultiVectorDistObjectKernels_MP_Vector.hpp"

namespace Tpetra {
namespace KokkosRefactor {
namespace Details {

  // Functors for implementing packAndPrepare and unpackAndCombine
  // through parallel_for

  template <typename DT, typename DL, typename DD, typename DM,
            typename ST, typename SL, typename SD, typename SM,
            typename IdxView>
  struct PackArraySingleColumn<
    Kokkos::View<DT,DL,DD,DM,Kokkos::Impl::ViewPCEContiguous>,
    Kokkos::View<ST,SL,SD,SM,Kokkos::Impl::ViewPCEContiguous>,
    IdxView >
  {
    typedef Kokkos::Impl::ViewPCEContiguous Spec;
    typedef Kokkos::View<DT,DL,DD,DM,Spec> DstView;
    typedef Kokkos::View<ST,SL,SD,SM,Spec> SrcView;
    typedef typename DstView::execution_space execution_space;
    typedef typename execution_space::size_type size_type;

    static const unsigned BlockSize = 32;

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
      for (size_type i=tidx; i<dst.sacado_size(); i+=BlockSize)
        dst(k).fastAccessCoeff(i) = src(idx(k), col).fastAccessCoeff(i);
    }

    static void pack(const DstView& dst,
                     const SrcView& src,
                     const IdxView& idx,
                     size_t col) {
      if ( Details::device_is_cuda<execution_space>::value )
        Kokkos::parallel_for(
          Kokkos::MPVectorWorkConfig<execution_space>( idx.size(), BlockSize ),
          PackArraySingleColumn(dst,src,idx,col) );
      else
         Kokkos::parallel_for( idx.size(),
                               PackArraySingleColumn(dst,src,idx,col) );
    }
  };

  template <typename DT, typename DL, typename DD, typename DM,
            typename ST, typename SL, typename SD, typename SM,
            typename IdxView>
  struct PackArrayMultiColumn<
    Kokkos::View<DT,DL,DD,DM,Kokkos::Impl::ViewPCEContiguous>,
    Kokkos::View<ST,SL,SD,SM,Kokkos::Impl::ViewPCEContiguous>,
    IdxView >
  {
    typedef Kokkos::Impl::ViewPCEContiguous Spec;
    typedef Kokkos::View<DT,DL,DD,DM,Spec> DstView;
    typedef Kokkos::View<ST,SL,SD,SM,Spec> SrcView;
    typedef typename DstView::execution_space execution_space;
    typedef typename execution_space::size_type size_type;

    static const unsigned BlockSize = 32;

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
        for (size_type i=tidx; i<dst.sacado_size(); i+=BlockSize)
          dst(offset + j).fastAccessCoeff(i) =
            src(localRow, j).fastAccessCoeff(i);
    }

    static void pack(const DstView& dst,
                     const SrcView& src,
                     const IdxView& idx,
                     size_t numCols) {
      if ( Details::device_is_cuda<execution_space>::value )
        Kokkos::parallel_for(
          Kokkos::MPVectorWorkConfig<execution_space>( idx.size(), BlockSize ),
          PackArrayMultiColumn(dst,src,idx,numCols) );
      else
        Kokkos::parallel_for( idx.size(),
                              PackArrayMultiColumn(dst,src,idx,numCols) );
    }
  };

  template <typename DT, typename DL, typename DD, typename DM,
            typename ST, typename SL, typename SD, typename SM,
            typename IdxView,
            typename ColView>
  struct PackArrayMultiColumnVariableStride<
    Kokkos::View<DT,DL,DD,DM,Kokkos::Impl::ViewPCEContiguous>,
    Kokkos::View<ST,SL,SD,SM,Kokkos::Impl::ViewPCEContiguous>,
    IdxView,
    ColView>
  {
    typedef Kokkos::Impl::ViewPCEContiguous Spec;
    typedef Kokkos::View<DT,DL,DD,DM,Spec> DstView;
    typedef Kokkos::View<ST,SL,SD,SM,Spec> SrcView;
    typedef typename DstView::execution_space execution_space;
    typedef typename execution_space::size_type size_type;

    static const unsigned BlockSize = 32;

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
        for (size_type i=tidx; i<dst.sacado_size(); i+=BlockSize)
          dst(offset + j).fastAccessCoeff(i) =
            src(localRow, col(j)).fastAccessCoeff(i);
    }

    static void pack(const DstView& dst,
                     const SrcView& src,
                     const IdxView& idx,
                     const ColView& col,
                     size_t numCols) {
      if ( Details::device_is_cuda<execution_space>::value )
        Kokkos::parallel_for(
          Kokkos::MPVectorWorkConfig<execution_space>( idx.size(), BlockSize ),
          PackArrayMultiColumnVariableStride(dst,src,idx,col,numCols) );
      else
        Kokkos::parallel_for( idx.size(),
                              PackArrayMultiColumnVariableStride(
                                dst,src,idx,col,numCols) );
    }
  };

  template <typename DT, typename DL, typename DD, typename DM,
            typename ST, typename SL, typename SD, typename SM,
            typename IdxView, typename Op>
  struct UnpackArrayMultiColumn<
    Kokkos::View<DT,DL,DD,DM,Kokkos::Impl::ViewPCEContiguous>,
    Kokkos::View<ST,SL,SD,SM,Kokkos::Impl::ViewPCEContiguous>,
    IdxView,
    Op >
  {
    typedef Kokkos::Impl::ViewPCEContiguous Spec;
    typedef Kokkos::View<DT,DL,DD,DM,Spec> DstView;
    typedef Kokkos::View<ST,SL,SD,SM,Spec> SrcView;
    typedef typename DstView::execution_space execution_space;
    typedef typename execution_space::size_type size_type;

    static const unsigned BlockSize = 32;

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

    KOKKOS_INLINE_FUNCTION
    void operator()( const size_type k, const size_type tidx ) const {
      const typename IdxView::value_type localRow = idx(k);
      const size_t offset = k*numCols;
      for (size_t j = 0; j < numCols; ++j)
        for (size_type i=tidx; i<dst.sacado_size(); i+=BlockSize)
          op( dst(localRow,j).fastAccessCoeff(i),
              src(offset+j).fastAccessCoeff(i) );
    }

    static void unpack(const DstView& dst,
                       const SrcView& src,
                       const IdxView& idx,
                       const Op& op,
                       size_t numCols) {
      if ( Details::device_is_cuda<execution_space>::value )
        Kokkos::parallel_for(
          Kokkos::MPVectorWorkConfig<execution_space>( idx.size(), BlockSize ),
          UnpackArrayMultiColumn(dst,src,idx,op,numCols) );
      else
        Kokkos::parallel_for( idx.size(),
                              UnpackArrayMultiColumn(dst,src,idx,op,numCols) );
    }
  };

  template <typename DT, typename DL, typename DD, typename DM,
            typename ST, typename SL, typename SD, typename SM,
            typename IdxView, typename ColView, typename Op>
  struct UnpackArrayMultiColumnVariableStride<
    Kokkos::View<DT,DL,DD,DM,Kokkos::Impl::ViewPCEContiguous>,
    Kokkos::View<ST,SL,SD,SM,Kokkos::Impl::ViewPCEContiguous>,
    IdxView,
    ColView,
    Op>
  {
    typedef Kokkos::Impl::ViewPCEContiguous Spec;
    typedef Kokkos::View<DT,DL,DD,DM,Spec> DstView;
    typedef Kokkos::View<ST,SL,SD,SM,Spec> SrcView;
    typedef typename DstView::execution_space execution_space;
    typedef typename execution_space::size_type size_type;

    static const unsigned BlockSize = 32;

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

    KOKKOS_INLINE_FUNCTION
    void operator()( const size_type k, const size_type tidx ) const {
      const typename IdxView::value_type localRow = idx(k);
      const size_t offset = k*numCols;
      for (size_t j = 0; j < numCols; ++j)
        for (size_type i=tidx; i<dst.sacado_size(); i+=BlockSize)
          op( dst(localRow,col(j)).fastAccessCoeff(i),
              src(offset+j).fastAccessCoeff(i) );
    }

    static void unpack(const DstView& dst,
                       const SrcView& src,
                       const IdxView& idx,
                       const ColView& col,
                       const Op& op,
                       size_t numCols) {
      if ( Details::device_is_cuda<execution_space>::value )
        Kokkos::parallel_for(
          Kokkos::MPVectorWorkConfig<execution_space>( idx.size(), BlockSize ),
          UnpackArrayMultiColumnVariableStride(dst,src,idx,col,op,numCols) );
      else
        Kokkos::parallel_for( idx.size(),
                              UnpackArrayMultiColumnVariableStride(
                                dst,src,idx,col,op,numCols) );
    }
  };

  template <typename DT, typename DL, typename DD, typename DM,
            typename ST, typename SL, typename SD, typename SM,
            typename DstIdxView, typename SrcIdxView>
  struct PermuteArrayMultiColumn<
    Kokkos::View<DT,DL,DD,DM,Kokkos::Impl::ViewPCEContiguous>,
    Kokkos::View<ST,SL,SD,SM,Kokkos::Impl::ViewPCEContiguous>,
    DstIdxView,
    SrcIdxView>
  {
    typedef Kokkos::Impl::ViewPCEContiguous Spec;
    typedef Kokkos::View<DT,DL,DD,DM,Spec> DstView;
    typedef Kokkos::View<ST,SL,SD,SM,Spec> SrcView;
    typedef typename DstView::execution_space execution_space;
    typedef typename execution_space::size_type size_type;

    static const unsigned BlockSize = 32;

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
        for (size_type i=tidx; i<dst.sacado_size(); i+=BlockSize)
          dst(toRow, j).fastAccessCoeff(i) =
            src(fromRow, j).fastAccessCoeff(i);
    }

    static void permute(const DstView& dst,
                        const SrcView& src,
                        const DstIdxView& dst_idx,
                        const SrcIdxView& src_idx,
                        size_t numCols) {
      const size_type n = std::min( dst_idx.size(), src_idx.size() );
      if ( Details::device_is_cuda<execution_space>::value )
        Kokkos::parallel_for(
          Kokkos::MPVectorWorkConfig<execution_space>( n, BlockSize ),
          PermuteArrayMultiColumn(dst,src,dst_idx,src_idx,numCols) );
      else
        Kokkos::parallel_for(
          n, PermuteArrayMultiColumn(dst,src,dst_idx,src_idx,numCols) );
    }
  };

  template <typename DT, typename DL, typename DD, typename DM,
            typename ST, typename SL, typename SD, typename SM,
            typename DstIdxView, typename SrcIdxView,
            typename DstColView, typename SrcColView>
  struct PermuteArrayMultiColumnVariableStride<
    Kokkos::View<DT,DL,DD,DM,Kokkos::Impl::ViewPCEContiguous>,
    Kokkos::View<ST,SL,SD,SM,Kokkos::Impl::ViewPCEContiguous>,
    DstIdxView, SrcIdxView,
    DstColView, SrcColView >
  {
    typedef Kokkos::Impl::ViewPCEContiguous Spec;
    typedef Kokkos::View<DT,DL,DD,DM,Spec> DstView;
    typedef Kokkos::View<ST,SL,SD,SM,Spec> SrcView;
    typedef typename DstView::execution_space execution_space;
    typedef typename execution_space::size_type size_type;

    static const unsigned BlockSize = 32;

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
        for (size_type i=tidx; i<dst.sacado_size(); i+=BlockSize)
          dst(toRow, dst_col(j)).fastAccessCoeff(i) =
            src(fromRow, src_col(j)).fastAccessCoeff(i);
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
          Kokkos::MPVectorWorkConfig<execution_space>( n, BlockSize ),
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

#endif // TPETRA_KOKKOS_REFACTOR_DETAILS_MULTI_VECTOR_DIST_OBJECT_KERNELS_UQ_PCE_HPP
