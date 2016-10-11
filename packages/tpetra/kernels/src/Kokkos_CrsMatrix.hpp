/*
//@HEADER
// ************************************************************************
//
//          Kokkos: Node API and Parallel Node Kernels
//              Copyright (2008) Sandia Corporation
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
//@HEADER
*/
#ifndef KOKKOS_CRSMATRIX_H_
#define KOKKOS_CRSMATRIX_H_

/// \file Kokkos_CrsMatrix.hpp
/// \brief Local sparse matrix interface
/// \warning Do NOT include this file.  This file is DEPRECATED!!!
///   Include Kokkos_Sparse_CrsMatrix.hpp instead.

#include <algorithm>

#include <Kokkos_Core.hpp>
#include <Kokkos_ArithTraits.hpp>
#include <Kokkos_StaticCrsGraph.hpp>

#ifdef KOKKOS_USE_CUSPARSE
#  include <cusparse.h>
#  include <Kokkos_CrsMatrix_CuSparse.hpp>
#endif // KOKKOS_USE_CUSPARSE

#ifdef KOKKOS_USE_MKL
#  include <mkl.h>
#  include <mkl_spblas.h>
#  include <Kokkos_CrsMatrix_MKL.hpp>
#endif // KOKKOS_USE_MKL

//#include <Kokkos_Vectorization.hpp>
#include <impl/Kokkos_Error.hpp>
#include <Kokkos_Sparse_CrsMatrix.hpp>
namespace Kokkos {
#if true
  using KokkosSparse::CrsMatrix;
  using KokkosSparse::RowsPerThread;
  using KokkosSparse::SparseRowView;
  using KokkosSparse::SparseRowViewConst;
  using KokkosSparse::DeviceConfig;
#endif // true

template<class DeviceType, typename ScalarType, int NNZPerRow=27>
struct MV_MultiplyShflThreadsPerRow {
private:
  typedef typename Kokkos::Impl::remove_const< ScalarType >::type value_type;

  // The shuffle operation only works with CUDA, and only works for
  // certain ScalarType types.
#ifdef KOKKOS_HAVE_CUDA
  enum { shfl_possible =
    Kokkos::Impl::is_same< DeviceType , Kokkos::Cuda >::value &&
    (
      Kokkos::Impl::is_same< value_type , unsigned int >::value ||
      Kokkos::Impl::is_same< value_type , int >::value ||
      Kokkos::Impl::is_same< value_type , float >::value ||
      Kokkos::Impl::is_same< value_type , double >::value
    )};
#else // NOT KOKKOS_HAVE_CUDA
  enum { shfl_possible = 0 };
#endif // KOKKOS_HAVE_CUDA

public:

#if defined( __CUDA_ARCH__ )
  enum { device_value = shfl_possible && ( 300 <= __CUDA_ARCH__ ) ?
         (NNZPerRow<8?2:
         (NNZPerRow<16?4:
         (NNZPerRow<32?8:
         (NNZPerRow<64?16:
         32))))
         :1 };
#else
  enum { device_value = 1 };
#endif

#ifdef KOKKOS_HAVE_CUDA
  inline static int host_value()
    { return shfl_possible && ( 300 <= Kokkos::Cuda::device_arch() ) ?
         (NNZPerRow<8?2:
         (NNZPerRow<16?4:
         (NNZPerRow<32?8:
         (NNZPerRow<64?16:
         32))))
         :1; }
#else // NOT KOKKOS_HAVE_CUDA
  inline static int host_value() { return 1; }
#endif // KOKKOS_HAVE_CUDA
};

//----------------------------------------------------------------------------

template<class RangeVector,
         class CrsMatrix,
         class DomainVector,
         class CoeffVector1,
         class CoeffVector2,
         int doalpha,
         int dobeta>
struct MV_MultiplyFunctor {
  typedef typename CrsMatrix::execution_space execution_space;
  typedef typename CrsMatrix::size_type       size_type;
  typedef typename CrsMatrix::ordinal_type    ordinal_type;
  typedef typename CrsMatrix::non_const_value_type value_type;
  typedef typename Kokkos::View<value_type*, execution_space> range_values;
  typedef typename Kokkos::TeamPolicy<execution_space>        team_policy;
  typedef typename team_policy::member_type                   team_member;

  CoeffVector1 beta;
  CoeffVector2 alpha;
  CrsMatrix  m_A;
  DomainVector  m_x;
  RangeVector  m_y;
  /// \brief The number of columns in the input and output MultiVectors.
  ///
  /// Its approxpriate type is therefore size_type, but we don't
  /// expect the input and output MultiVectors to have more columns
  /// than the sparse matrix has rows or columns.  Thus, we prefer the
  /// (likely both smaller and signed, vs. the larger and likely
  /// unsigned) size_type.
  ordinal_type n;
  int rows_per_thread;

  MV_MultiplyFunctor (const CoeffVector1 beta_,
                      const CoeffVector2 alpha_,
                      const CrsMatrix m_A_,
                      const DomainVector m_x_,
                      const RangeVector m_y_,
                      const ordinal_type n_,
                      const int rows_per_thread_) :
      beta (beta_), alpha (alpha_),
      m_A (m_A_), m_x (m_x_), m_y (m_y_), n (n_),
      rows_per_thread (rows_per_thread_)
  {}

  template<int UNROLL>
  KOKKOS_INLINE_FUNCTION void
  strip_mine (const team_member& dev, const size_type& iRow, const size_type& kk) const
  {
    value_type sum[UNROLL];

#ifdef KOKKOS_HAVE_PRAGMA_IVDEP
#pragma ivdep
#endif
#ifdef KOKKOS_HAVE_PRAGMA_UNROLL
#pragma unroll
#endif
    for (int k = 0; k < UNROLL; ++k) {
      // NOTE (mfh 09 Aug 2013) This requires that assignment from int
      // (in this case, 0) to value_type be defined.  It's not for
      // types like arprec and dd_real.
      //
      // mfh 29 Sep 2013: On the other hand, arprec and dd_real won't
      // work on CUDA devices anyway, since their methods aren't
      // device functions.  arprec has other issues (e.g., dynamic
      // memory allocation, and the array-of-structs memory layout
      // which is unfavorable to GPUs), but could be dealt with in the
      // same way as Sacado's AD types.
      sum[k] = 0;
    }

    const SparseRowViewConst<CrsMatrix> row = m_A.rowConst (iRow);

    // NOTE (mfh 20 Mar 2015) Unfortunately, Kokkos::Vectorization
    // lacks a typedef for determining the type of the return value of
    // begin().  I know that it returns int now, but this may change
    // at some point.
    //
    // The correct type of iEntry is ordinal_type.  This is because we
    // assume that rows have no duplicate entries.  As a result, a row
    // cannot have more entries than the number of columns in the
    // matrix.

#ifdef KOKKOS_HAVE_PRAGMA_IVDEP
#pragma ivdep
#endif
#ifdef KOKKOS_HAVE_PRAGMA_UNROLL
#pragma unroll
#endif
#ifdef KOKKOS_HAVE_PRAGMA_LOOPCOUNT
#pragma loop count (15)
#endif
#ifdef __CUDA_ARCH__
        for (ordinal_type iEntry = static_cast<ordinal_type> (threadIdx.x);
             iEntry < static_cast<ordinal_type> (row.length);
             iEntry += static_cast<ordinal_type> (blockDim.x)) {
#else
        for (ordinal_type iEntry = 0;
             iEntry < static_cast<ordinal_type> (row.length);
             iEntry ++) {
#endif
      const value_type val = row.value(iEntry);
      const ordinal_type ind = row.colidx(iEntry);

#ifdef KOKKOS_HAVE_PRAGMA_UNROLL
#pragma unroll
#endif
      for (int k = 0; k < UNROLL; ++k) {
        sum[k] +=  val * m_x(ind, kk + k);
      }
    }

    if (doalpha == -1) {
      for (int ii=0; ii < UNROLL; ++ii) {
        value_type sumt=sum[ii];
#if defined(__CUDA_ARCH__) && defined(KOKKOS_HAVE_CUDA)
        if (blockDim.x > 1)
          sumt += shfl_down(sumt, 1,blockDim.x);
        if (blockDim.x > 2)
          sumt += shfl_down(sumt, 2,blockDim.x);
        if (blockDim.x > 4)
          sumt += shfl_down(sumt, 4,blockDim.x);
        if (blockDim.x > 8)
          sumt += shfl_down(sumt, 8,blockDim.x);
        if (blockDim.x > 16)
          sumt += shfl_down(sumt, 16,blockDim.x);
#endif // defined(__CUDA_ARCH__) && defined(KOKKOS_HAVE_CUDA)
        sum[ii] = - sumt;
      }
    }
    else {
      for (int ii=0; ii < UNROLL; ++ii) {
        value_type sumt = sum[ii];
#if defined(__CUDA_ARCH__) && defined(KOKKOS_HAVE_CUDA)
        if (blockDim.x > 1)
          sumt += shfl_down(sumt, 1,blockDim.x);
        if (blockDim.x > 2)
          sumt += shfl_down(sumt, 2,blockDim.x);
        if (blockDim.x > 4)
          sumt += shfl_down(sumt, 4,blockDim.x);
        if (blockDim.x > 8)
          sumt += shfl_down(sumt, 8,blockDim.x);
        if (blockDim.x > 16)
          sumt += shfl_down(sumt, 16,blockDim.x);
#endif // defined(__CUDA_ARCH__) && defined(KOKKOS_HAVE_CUDA)
        sum[ii] = sumt;
      }
    }

#if defined(__CUDA_ARCH__) && defined(KOKKOS_HAVE_CUDA)
    if (threadIdx.x==0) {
#else
    if (true) {
#endif // defined(__CUDA_ARCH__) && defined(KOKKOS_HAVE_CUDA)
      if (doalpha * doalpha != 1) {
#ifdef KOKKOS_HAVE_PRAGMA_IVDEP
#pragma ivdep
#endif
#ifdef KOKKOS_HAVE_PRAGMA_UNROLL
#pragma unroll
#endif
        for (int k = 0; k < UNROLL; ++k) {
          sum[k] *= alpha(kk + k);
        }
      }

      if (dobeta == 0) {
#ifdef KOKKOS_HAVE_PRAGMA_IVDEP
#pragma ivdep
#endif
#ifdef KOKKOS_HAVE_PRAGMA_UNROLL
#pragma unroll
#endif
        for (int k = 0; k < UNROLL; ++k) {
          m_y(iRow, kk + k) = sum[k];
        }
      } else if (dobeta == 1) {
#ifdef KOKKOS_HAVE_PRAGMA_IVDEP
#pragma ivdep
#endif
#ifdef KOKKOS_HAVE_PRAGMA_UNROLL
#pragma unroll
#endif
        for (int k = 0; k < UNROLL; ++k) {
          m_y(iRow, kk + k) += sum[k];
        }
      } else if (dobeta == -1) {
#ifdef KOKKOS_HAVE_PRAGMA_IVDEP
#pragma ivdep
#endif
#ifdef KOKKOS_HAVE_PRAGMA_UNROLL
#pragma unroll
#endif
        for (int k = 0; k < UNROLL; ++k) {
          m_y(iRow, kk + k) = -m_y(iRow, kk + k) +  sum[k];
        }
      } else {
#ifdef KOKKOS_HAVE_PRAGMA_IVDEP
#pragma ivdep
#endif
#ifdef KOKKOS_HAVE_PRAGMA_UNROLL
#pragma unroll
#endif
        for (int k = 0; k < UNROLL; ++k) {
          m_y(iRow, kk + k) = beta(kk + k) * m_y(iRow, kk + k) + sum[k] ;
        }
      }
    }
  }

  KOKKOS_INLINE_FUNCTION void
  strip_mine_1 (const team_member& dev, const size_type& iRow) const
  {
    value_type sum = 0;

    const SparseRowViewConst<CrsMatrix> row = m_A.rowConst (iRow);

    // NOTE (mfh 20 Mar 2015) Unfortunately, Kokkos::Vectorization
    // lacks a typedef for determining the type of the return value of
    // begin().  I know that it returns int now, but this may change
    // at some point.
    //
    // The correct type of iEntry is ordinal_type.  This is because we
    // assume that rows have no duplicate entries.  As a result, a row
    // cannot have more entries than the number of columns in the
    // matrix.

#ifdef KOKKOS_HAVE_PRAGMA_IVDEP
#pragma ivdep
#endif
#ifdef KOKKOS_HAVE_PRAGMA_UNROLL
#pragma unroll
#endif
#ifdef KOKKOS_HAVE_PRAGMA_LOOPCOUNT
#pragma loop count (15)
#endif
#ifdef __CUDA_ARCH__
    for (ordinal_type iEntry = static_cast<ordinal_type> (threadIdx.x);
       iEntry < static_cast<ordinal_type> (row.length);
       iEntry += static_cast<ordinal_type> (blockDim.x)) {
#else
    for (ordinal_type iEntry = 0;
       iEntry < static_cast<ordinal_type> (row.length);
       iEntry ++) {
#endif
      sum += row.value(iEntry) * m_x(row.colidx(iEntry),0);
    }
#if defined(__CUDA_ARCH__) && defined(KOKKOS_HAVE_CUDA)
    if (blockDim.x > 1)
      sum += shfl_down(sum, 1,blockDim.x);
    if (blockDim.x > 2)
      sum += shfl_down(sum, 2,blockDim.x);
    if (blockDim.x > 4)
      sum += shfl_down(sum, 4,blockDim.x);
    if (blockDim.x > 8)
      sum += shfl_down(sum, 8,blockDim.x);
    if (blockDim.x > 16)
      sum += shfl_down(sum, 16,blockDim.x);
#endif // defined(__CUDA_ARCH__) && defined(KOKKOS_HAVE_CUDA)

#ifdef __CUDA_ARCH__
    if (threadIdx.x==0) {
#else
    if (true) {
#endif
      if (doalpha == -1) {
        sum *= value_type(-1);
      } else if (doalpha * doalpha != 1) {
        sum *= alpha(0);
      }

      if (dobeta == 0) {
        m_y(iRow, 0) = sum ;
      } else if (dobeta == 1) {
        m_y(iRow, 0) += sum ;
      } else if (dobeta == -1) {
        m_y(iRow, 0) = -m_y(iRow, 0) +  sum;
      } else {
        m_y(iRow, 0) = beta(0) * m_y(iRow, 0) + sum;
      }
    }
  }


  KOKKOS_INLINE_FUNCTION void
  operator() (const team_member& dev) const
  {
    for (int loop = 0; loop < rows_per_thread; ++loop) {

      // NOTE (mfh 20 Mar 2015) Unfortunately, Kokkos::Vectorization
      // lacks a typedef for determining the type of the return value
      // of global_thread_rank().  I know that it returns int now, but
      // this may change at some point.
      //
      // iRow represents a row of the matrix, so its correct type is
      // ordinal_type.

      const ordinal_type iRow = (dev.league_rank() * dev.team_size() + dev.team_rank())
                                * rows_per_thread + loop;
      if (iRow >= m_A.numRows ()) {
        return;
      }

      // mfh 20 Mar 2015: This relates to n, so its correct type is
      // ordinal_type.  Once we can use C++11 without protection, the
      // right thing to do would be to use decltype to pick up n's
      // type here, rather than assuming that it's ordinal_type.
      ordinal_type kk = 0;

#ifdef KOKKOS_FAST_COMPILE
      for (; kk + 4 <= n; kk += 4) {
        strip_mine<4>(dev, iRow, kk);
      }
      for( ; kk < n; ++kk) {
        strip_mine<1>(dev, iRow, kk);
      }
#else
#  ifdef __CUDA_ARCH__
      if ((n > 8) && (n % 8 == 1)) {
        strip_mine<9>(dev, iRow, kk);
        kk += 9;
      }
      for(; kk + 8 <= n; kk += 8)
        strip_mine<8>(dev, iRow, kk);
      if(kk < n)
        switch(n - kk) {
#  else // NOT a CUDA device
          if ((n > 16) && (n % 16 == 1)) {
            strip_mine<17>(dev, iRow, kk);
            kk += 17;
          }

          for (; kk + 16 <= n; kk += 16) {
            strip_mine<16>(dev, iRow, kk);
          }

          if(kk < n)
            switch(n - kk) {
            case 15:
              strip_mine<15>(dev, iRow, kk);
              break;

            case 14:
              strip_mine<14>(dev, iRow, kk);
              break;

            case 13:
              strip_mine<13>(dev, iRow, kk);
              break;

            case 12:
              strip_mine<12>(dev, iRow, kk);
              break;

            case 11:
              strip_mine<11>(dev, iRow, kk);
              break;

            case 10:
              strip_mine<10>(dev, iRow, kk);
              break;

            case 9:
              strip_mine<9>(dev, iRow, kk);
              break;

            case 8:
              strip_mine<8>(dev, iRow, kk);
              break;
#  endif // __CUDA_ARCH__
            case 7:
              strip_mine<7>(dev, iRow, kk);
              break;

            case 6:
              strip_mine<6>(dev, iRow, kk);
              break;

            case 5:
              strip_mine<5>(dev, iRow, kk);
              break;

            case 4:
              strip_mine<4>(dev, iRow, kk);
              break;

            case 3:
              strip_mine<3>(dev, iRow, kk);
              break;

            case 2:
              strip_mine<2>(dev, iRow, kk);
              break;

            case 1:
              strip_mine_1(dev, iRow);
              break;
            }
#endif // KOKKOS_FAST_COMPILE
        }
    }
  };


  template<class RangeVector,
           class CrsMatrix,
           class DomainVector,
           class CoeffVector1,
           class CoeffVector2,
           int doalpha,
           int dobeta>
  struct MV_MultiplySingleFunctor {
    typedef typename CrsMatrix::execution_space      execution_space;
    typedef typename CrsMatrix::ordinal_type         ordinal_type;
    typedef typename CrsMatrix::non_const_value_type value_type;
    typedef typename Kokkos::View<value_type*, typename CrsMatrix::execution_space> range_values;
    typedef typename Kokkos::TeamPolicy<execution_space> team_policy;
    typedef typename team_policy::member_type            team_member;

    CoeffVector1 beta;
    CoeffVector2 alpha;
    CrsMatrix  m_A;
    DomainVector m_x;
    RangeVector m_y;

    ordinal_type rows_per_thread;

    MV_MultiplySingleFunctor (const CoeffVector1 beta_,
                              const CoeffVector2 alpha_,
                              const CrsMatrix m_A_,
                              const DomainVector m_x_,
                              const RangeVector m_y_,
                              const int rows_per_thread_) :
      beta (beta_), alpha (alpha_),
      m_A (m_A_), m_x (m_x_), m_y (m_y_),
      rows_per_thread (rows_per_thread_)
    {}

    KOKKOS_INLINE_FUNCTION void
    operator() (const team_member& dev) const
    {
      // This should be a thread loop as soon as we can use C++11.
      //
      // FIXME (mfh 20 Mar 2015, 11 Apr 2015) The correct type of
      // 'loop' should be ordinal_type, not int.  Ditto for
      // rows_per_thread.  The cast avoids a build warning.
      for (int loop = 0; loop < static_cast<int> (rows_per_thread); ++loop) {
        // NOTE (mfh 20 Mar 2015) Unfortunately, Kokkos::Vectorization
        // lacks a typedef for determining the type of the return
        // value of global_thread_rank().  I know that it returns int
        // now, but this may change at some point.
        //
        // iRow represents a row of the matrix, so its correct type is
        // ordinal_type.
        const ordinal_type iRow = (dev.league_rank() * dev.team_size() + dev.team_rank())
                                  * rows_per_thread + loop;
        if (iRow >= m_A.numRows ()) {
          return;
        }
        const SparseRowViewConst<CrsMatrix> row = m_A.rowConst (iRow);
        const ordinal_type row_length = static_cast<ordinal_type> (row.length);
        value_type sum = 0;

        // Use explicit Cuda below to avoid C++11 for now. This should be a vector reduce loop !
        #ifdef KOKKOS_HAVE_PRAGMA_IVDEP
        #pragma ivdep
        #endif
        #ifdef KOKKOS_HAVE_PRAGMA_UNROLL
        #pragma unroll
        #endif
        #ifdef KOKKOS_HAVE_PRAGMA_LOOPCOUNT
        #pragma loop count (15)
        #endif
#ifdef __CUDA_ARCH__
        for (ordinal_type iEntry = static_cast<ordinal_type> (threadIdx.x);
             iEntry < static_cast<ordinal_type> (row_length);
             iEntry += static_cast<ordinal_type> (blockDim.x)) {
#else
        for (ordinal_type iEntry = 0;
             iEntry < static_cast<ordinal_type> (row_length);
             iEntry ++) {
#endif
          sum += row.value(iEntry) * m_x(row.colidx(iEntry));
        }

#if defined(__CUDA_ARCH__) && defined(KOKKOS_HAVE_CUDA)
        if (blockDim.x > 1)
          sum += shfl_down(sum, 1,blockDim.x);
        if (blockDim.x > 2)
          sum += shfl_down(sum, 2,blockDim.x);
        if (blockDim.x > 4)
          sum += shfl_down(sum, 4,blockDim.x);
        if (blockDim.x > 8)
          sum += shfl_down(sum, 8,blockDim.x);
        if (blockDim.x > 16)
          sum += shfl_down(sum, 16,blockDim.x);

        if (threadIdx.x==0) {
#else
        if (true) {
#endif // defined(__CUDA_ARCH__) && defined(KOKKOS_HAVE_CUDA)
          if (doalpha == -1) {
            sum *= value_type(-1);
          } else if (doalpha * doalpha != 1) {
            sum *= alpha(0);
          }

          if (dobeta == 0) {
            m_y(iRow) = sum ;
          } else if (dobeta == 1) {
            m_y(iRow) += sum ;
          } else if (dobeta == -1) {
            m_y(iRow) = -m_y(iRow) +  sum;
          } else {
            m_y(iRow) = beta(0) * m_y(iRow) + sum;
          }
        }
      }
    }
  };

  namespace Impl {

    template <class RangeVector,
              class CrsMatrix,
              class DomainVector,
              class CoeffVector1,
              class CoeffVector2>
    void
    MV_Multiply_Check_Compatibility (const CoeffVector1 &betav,
                                     const RangeVector &y,
                                     const CoeffVector2 &alphav,
                                     const CrsMatrix &A,
                                     const DomainVector &x,
                                     const int& doalpha,
                                     const int& dobeta)
    {
      typename DomainVector::size_type numVecs = x.dimension_1();
      typename DomainVector::size_type numRows = A.numRows();
      typename DomainVector::size_type numCols = A.numCols();

      if (y.dimension_1() != numVecs) {
        std::ostringstream msg;
        msg << "Error in CRSMatrix - Vector Multiply (y = by + aAx): 2nd dimensions of y and x do not match\n";
        msg << "\t Labels are: y(" << y.tracker().label() << ") b("
            << betav.tracker().label() << ") a("
            << alphav.tracker().label() << ") x("
            << A.values.tracker().label() << ") x("
            << x.tracker().label() << ")\n";
        msg << "\t Dimensions are: y(" << y.dimension_0() << "," << y.dimension_1() << ") x(" << x.dimension_0() << "," << x.dimension_1() << ")\n";
        Impl::throw_runtime_exception( msg.str() );
      }
      if (numRows > y.dimension_0()) {
        std::ostringstream msg;
        msg << "Error in CRSMatrix - Vector Multiply (y = by + aAx): dimensions of y and A do not match\n";
        msg << "\t Labels are: y(" << y.tracker().label() << ") b("
            << betav.tracker().label() << ") a("
            << alphav.tracker().label() << ") x("
            << A.values.tracker().label() << ") x("
            << x.tracker().label() << ")\n";
        msg << "\t Dimensions are: y(" << y.dimension_0() << "," << y.dimension_1() << ") A(" << A.numRows() << "," << A.numCols() << ")\n";
        Impl::throw_runtime_exception( msg.str() );
      }
      if (numCols > x.dimension_0()) {
        std::ostringstream msg;
        msg << "Error in CRSMatrix - Vector Multiply (y = by + aAx): dimensions of x and A do not match\n";
        msg << "\t Labels are: y(" << y.tracker().label() << ") b("
            << betav.tracker().label() << ") a("
            << alphav.tracker().label() << ") x("
            << A.values.tracker().label() << ") x("
            << x.tracker().label() << ")\n";
        msg << "\t Dimensions are: x(" << x.dimension_0() << "," << x.dimension_1() << ") A(" << A.numRows() << "," << A.numCols() << ")\n";
        Impl::throw_runtime_exception( msg.str() );
      }
      if (dobeta==2) {
        if (betav.dimension_0()!=numVecs) {
          std::ostringstream msg;
          msg << "Error in CRSMatrix - Vector Multiply (y = by + aAx): 2nd dimensions of y and b do not match\n";
          msg << "\t Labels are: y(" << y.tracker().label() << ") b("
              << betav.tracker().label() << ") a("
              << alphav.tracker().label() << ") x("
              << A.values.tracker().label() << ") x("
              << x.tracker().label() << ")\n";
          msg << "\t Dimensions are: y(" << y.dimension_0() << "," << y.dimension_1() << ") b(" << betav.dimension_0() << ")\n";
          Impl::throw_runtime_exception( msg.str() );
        }
      }
      if(doalpha==2) {
        if(alphav.dimension_0()!=numVecs) {
          std::ostringstream msg;
          msg << "Error in CRSMatrix - Vector Multiply (y = by + aAx): 2nd dimensions of x and b do not match\n";
          msg << "\t Labels are: y(" << y.tracker().label() << ") b("
              << betav.tracker().label() << ") a("
              << alphav.tracker().label() << ") x("
              << A.values.tracker().label() << ") x("
              << x.tracker().label() << ")\n";
          msg << "\t Dimensions are: x(" << x.dimension_0() << "," << x.dimension_1() << ") b(" << betav.dimension_0() << ")\n";
          Impl::throw_runtime_exception( msg.str() );
        }
      }
    }
  } // namespace Impl

  // This TransposeFunctor is functional, but not necessarily performant.
  template<class RangeVector,
           class CrsMatrix,
           class DomainVector,
           class CoeffVector1,
           class CoeffVector2,
           int doalpha,
           int dobeta,
           bool conjugate = false,
           int NNZPerRow = 27>
  struct MV_MultiplyTransposeFunctor {
    typedef typename CrsMatrix::execution_space                 execution_space;
    typedef typename CrsMatrix::ordinal_type                    ordinal_type;
    typedef typename CrsMatrix::size_type                       size_type;
    typedef typename CrsMatrix::non_const_value_type            value_type;
    typedef typename Kokkos::View<value_type*, execution_space> range_values;

    typedef MV_MultiplyShflThreadsPerRow<execution_space, value_type, NNZPerRow> ShflThreadsPerRow;

    CoeffVector1 beta;
    CoeffVector2 alpha;
    CrsMatrix  m_A ;
    DomainVector  m_x ;
    RangeVector  m_y ;
    ordinal_type n;

    // This is an iteration over rows of the matrix (modulo the
    // shuffle width), so the correct type of i is ordinal_type.
    KOKKOS_INLINE_FUNCTION
    void operator() (const ordinal_type i) const {
      typedef Kokkos::Details::ArithTraits<value_type> ATV;

      const ordinal_type iRow = i / ShflThreadsPerRow::device_value;
      const int lane = static_cast<int> (i) % ShflThreadsPerRow::device_value;
      const SparseRowViewConst<CrsMatrix> row = m_A.rowConst (iRow);

      for (ordinal_type iEntry = static_cast<ordinal_type> (lane);
           iEntry < static_cast<ordinal_type> (row.length);
           iEntry += static_cast<ordinal_type> (ShflThreadsPerRow::device_value)) {
        const value_type val = conjugate ?
          ATV::conj (row.value(iEntry)) :
          row.value(iEntry);
        const ordinal_type ind = row.colidx(iEntry);

        if (doalpha != 1) {
#ifdef KOKKOS_HAVE_PRAGMA_UNROLL
#pragma unroll
#endif
          for (ordinal_type k = 0; k < n; ++k) {
            atomic_add (&m_y(ind,k), value_type(alpha(k) * val * m_x(iRow, k)));
          }
        } else {
#ifdef KOKKOS_HAVE_PRAGMA_UNROLL
#pragma unroll
#endif
          for (ordinal_type k = 0; k < n; ++k) {
            atomic_add (&m_y(ind,k), value_type(val * m_x(iRow, k)));
          }
        }
      }
    }
  };

  // This TansposeFunctor is functional, but not necessarily performant.
  template<class RangeVector,
           class CrsMatrix,
           class DomainVector,
           class CoeffVector1,
           class CoeffVector2,
           int doalpha,
           int dobeta,
           bool conjugate = false,
           int NNZPerRow = 27 >
  struct MV_MultiplyTransposeSingleFunctor {
    typedef typename CrsMatrix::execution_space      execution_space;
    typedef typename CrsMatrix::ordinal_type         ordinal_type;
    typedef typename CrsMatrix::non_const_value_type value_type;
    typedef typename Kokkos::View<value_type*, execution_space> range_values;

    typedef MV_MultiplyShflThreadsPerRow< execution_space , value_type , NNZPerRow > ShflThreadsPerRow ;

    CoeffVector1 beta;
    CoeffVector2 alpha;
    CrsMatrix  m_A ;
    DomainVector  m_x ;
    RangeVector  m_y ;
    ordinal_type n;

    KOKKOS_INLINE_FUNCTION
    void operator() (const ordinal_type i) const {
      typedef Kokkos::Details::ArithTraits<value_type> ATV;

      const ordinal_type iRow = i / ShflThreadsPerRow::device_value;
      const int lane = static_cast<int> (i) % ShflThreadsPerRow::device_value;
      const SparseRowViewConst<CrsMatrix> row = m_A.rowConst (iRow);

      for (ordinal_type iEntry = static_cast<ordinal_type> (lane);
           iEntry < static_cast<ordinal_type> (row.length);
           iEntry += static_cast<ordinal_type> (ShflThreadsPerRow::device_value)) {
        const value_type val = conjugate ?
          ATV::conj (row.value(iEntry)) :
          row.value(iEntry);
        const ordinal_type ind = row.colidx(iEntry);

        if (doalpha != 1) {
          atomic_add (&m_y(ind), value_type(alpha(0) * val * m_x(iRow)));
        } else {
          atomic_add (&m_y(ind), value_type(val * m_x(iRow)));
        }
      }
    }
  };

template <class RangeVector,
          class TCrsMatrix,
          class DomainVector,
          class CoeffVector1,
          class CoeffVector2,
          int doalpha,
          int dobeta>
void
MV_MultiplyTranspose (typename Kokkos::Impl::enable_if<DomainVector::Rank == 2, const CoeffVector1>::type& betav,
                      const RangeVector &y,
                      const CoeffVector2 &alphav,
                      const TCrsMatrix &A,
                      const DomainVector &x,
                      const bool conjugate = false)
{
  typedef typename TCrsMatrix::ordinal_type ordinal_type;

  // FIXME (mfh 02 Jan 2015) Is numRows() always signed?  More
  // importantly, if the calling process owns zero rows in the row
  // Map, numRows() should return 0, not -1.
  //
  //Special case for zero Rows RowMap
  if (A.numRows () == static_cast<ordinal_type> (-1)) {
    return;
  }

  if (doalpha == 0) {
    if (dobeta == 2) {
      MV_MulScalar (y, betav, y);
    } else {
      MV_MulScalar (y, static_cast<typename RangeVector::const_value_type> (dobeta), y);
    }
    return;
  } else {
    typedef View< typename RangeVector::non_const_data_type ,
                  typename RangeVector::array_layout ,
                  typename RangeVector::execution_space ,
                  typename RangeVector::memory_traits >
    RangeVectorType;

    typedef View< typename DomainVector::const_data_type ,
                  typename DomainVector::array_layout ,
                  typename DomainVector::execution_space ,
                  Kokkos::MemoryRandomAccess >
    DomainVectorType;

    typedef View< typename CoeffVector1::const_data_type ,
                  typename CoeffVector1::array_layout ,
                  typename CoeffVector1::execution_space ,
                  Kokkos::MemoryRandomAccess >
    CoeffVector1Type;

    typedef View< typename CoeffVector2::const_data_type ,
                  typename CoeffVector2::array_layout ,
                  typename CoeffVector2::execution_space ,
                  Kokkos::MemoryRandomAccess >
    CoeffVector2Type;

    typedef CrsMatrix<typename TCrsMatrix::const_value_type,
                      typename TCrsMatrix::ordinal_type,
                      typename TCrsMatrix::execution_space,
                      typename TCrsMatrix::memory_traits,
                      typename TCrsMatrix::size_type> CrsMatrixType;

    // FIXME (mfh 20 Mar 2015) The dimension doesn't have type int.
    // On the other hand, this is the number of columns in the input
    // (and output) MultiVectors, so it's not likely to be large (the
    // typical case is < 100).
    int numVecs = x.dimension_1();
    CoeffVector1 beta = betav;
    CoeffVector2 alpha = alphav;

    if (doalpha != 2) {
      alpha = CoeffVector2("CrsMatrix::auto_a", numVecs);
      typename CoeffVector2::HostMirror h_a = Kokkos::create_mirror_view(alpha);
      typename CoeffVector2::value_type s_a = (typename CoeffVector2::value_type) doalpha;

      for (int i = 0; i < numVecs; ++i) {
        h_a(i) = s_a;
      }

      Kokkos::deep_copy (alpha, h_a);
    }

    if (dobeta != 2) {
      beta = CoeffVector1("CrsMatrix::auto_b", numVecs);
      typename CoeffVector1::HostMirror h_b = Kokkos::create_mirror_view(beta);
      typename CoeffVector1::value_type s_b = (typename CoeffVector1::value_type) dobeta;

      for (int i = 0; i < numVecs; ++i) {
        h_b(i) = s_b;
      }

      Kokkos::deep_copy (beta, h_b);
    }

    if (dobeta == 2) {
      MV_MulScalar (y, betav, y);
    } else {
      if (dobeta != 1) {
        MV_MulScalar (y, static_cast<typename RangeVector::const_value_type> (dobeta), y);
      }
    }

    const ordinal_type nrow = A.numRows ();

    if (conjugate) {
      typedef MV_MultiplyTransposeFunctor<RangeVectorType, CrsMatrixType,
                                          DomainVectorType, CoeffVector1Type,
                                          CoeffVector2Type, 2, 2, true> OpType;
      OpType op ;
      op.m_A = A;
      op.m_x = x;
      op.m_y = y;
      op.beta = beta;
      op.alpha = alpha;
      op.n = x.dimension_1();
      Kokkos::parallel_for (nrow * OpType::ShflThreadsPerRow::host_value (), op);
    }
    else {
      typedef MV_MultiplyTransposeFunctor<RangeVectorType, CrsMatrixType,
                                          DomainVectorType, CoeffVector1Type,
                                          CoeffVector2Type, 2, 2, false> OpType;
      OpType op ;
      op.m_A = A;
      op.m_x = x;
      op.m_y = y;
      op.beta = beta;
      op.alpha = alpha;
      op.n = x.dimension_1();
      Kokkos::parallel_for (nrow * OpType::ShflThreadsPerRow::host_value (), op);
    }

//#endif // KOKKOS_FAST_COMPILE
  }
}

template<class RangeVector,
         class CrsMatrix,
         class DomainVector,
         class CoeffVector1,
         class CoeffVector2>
void
MV_MultiplyTranspose (const CoeffVector1& betav,
                      const RangeVector& y,
                      const CoeffVector2& alphav,
                      const CrsMatrix& A,
                      const DomainVector& x,
                      int beta,
                      int alpha,
                      const bool conjugate = false)
{
  if (beta == 0) {
    if (alpha == 0) {
      MV_MultiplyTranspose<RangeVector, CrsMatrix, DomainVector, CoeffVector1,
                           CoeffVector2, 0, 0 > (betav, y, alphav, A, x, conjugate);
    }
    else if (alpha == 1) {
      MV_MultiplyTranspose<RangeVector, CrsMatrix, DomainVector, CoeffVector1,
                           CoeffVector2, 1, 0 > (betav, y, alphav, A, x, conjugate);
    }
    else if (alpha == -1) {
      MV_MultiplyTranspose<RangeVector, CrsMatrix, DomainVector, CoeffVector1,
                           CoeffVector2, -1, 0 > (betav, y, alphav, A, x, conjugate);
    }
    else {
      MV_MultiplyTranspose<RangeVector, CrsMatrix, DomainVector, CoeffVector1,
                           CoeffVector2, 2, 0 > (betav, y, alphav, A, x, conjugate);
    }
  } else if (beta == 1) {
    if (alpha == 0) {
      return;
    }
    else if (alpha == 1) {
      MV_MultiplyTranspose<RangeVector, CrsMatrix, DomainVector, CoeffVector1,
                           CoeffVector2, 1, 1 > (betav, y, alphav, A, x, conjugate);
    }
    else if (alpha == -1) {
      MV_MultiplyTranspose<RangeVector, CrsMatrix, DomainVector, CoeffVector1,
                           CoeffVector2, -1, 1 > (betav, y, alphav, A, x, conjugate);
    }
    else {
      MV_MultiplyTranspose<RangeVector, CrsMatrix, DomainVector, CoeffVector1,
                           CoeffVector2, 2, 1 > (betav, y, alphav, A, x, conjugate);
    }
  } else if (beta == -1) {
    if (alpha == 0) {
      MV_MultiplyTranspose<RangeVector, CrsMatrix, DomainVector, CoeffVector1,
                           CoeffVector2, 0, -1 > (betav, y, alphav, A, x, conjugate);
    }
    else if (alpha == 1) {
      MV_MultiplyTranspose<RangeVector, CrsMatrix, DomainVector, CoeffVector1,
                           CoeffVector2, 1, -1 > (betav, y, alphav, A, x, conjugate);
    }
    else if (alpha == -1) {
      MV_MultiplyTranspose<RangeVector, CrsMatrix, DomainVector, CoeffVector1,
                           CoeffVector2, -1, -1 > (betav, y, alphav, A, x, conjugate);
    }
    else {
      MV_MultiplyTranspose<RangeVector, CrsMatrix, DomainVector, CoeffVector1,
                           CoeffVector2, 2, -1 > (betav, y, alphav, A, x, conjugate);
    }
  } else {
    if (alpha == 0) {
      MV_MultiplyTranspose<RangeVector, CrsMatrix, DomainVector, CoeffVector1,
                           CoeffVector2, 0, 2 > (betav, y, alphav, A, x, conjugate);
    }
    else if (alpha == 1) {
      MV_MultiplyTranspose<RangeVector, CrsMatrix, DomainVector, CoeffVector1,
                           CoeffVector2, 1, 2 > (betav, y, alphav, A, x, conjugate);
    }
    else if (alpha == -1) {
      MV_MultiplyTranspose<RangeVector, CrsMatrix, DomainVector, CoeffVector1,
                           CoeffVector2, -1, 2 > (betav, y, alphav, A, x, conjugate);
    }
    else {
      MV_MultiplyTranspose<RangeVector, CrsMatrix, DomainVector, CoeffVector1,
                           CoeffVector2, 2, 2> (betav, y, alphav, A, x, conjugate);
    }
  }
}

template<class RangeVector, class CrsMatrix, class DomainVector>
void
MV_MultiplyTranspose (typename RangeVector::const_value_type s_b,
                      const RangeVector& y,
                      typename DomainVector::const_value_type s_a,
                      const CrsMatrix& A,
                      const DomainVector& x,
                      const bool conjugate = false)
{
/*#ifdef KOKKOS_USE_CUSPARSE
  if (MV_Multiply_Try_CuSparse (s_b, y, s_a, A, x, conjugate)) {
    return;
  }
#endif // KOKKOSE_USE_CUSPARSE
#ifdef KOKKOS_USE_MKL
  if (MV_Multiply_Try_MKL (s_b, y, s_a, A, x, conjugate)) {
    return;
  }
#endif // KOKKOS_USE_MKL*/
  typedef Kokkos::View<typename RangeVector::value_type*,
                       typename RangeVector::execution_space> aVector;
  aVector a;
  aVector b;
  int numVecs = x.dimension_1();

  if (s_b == 0) {
    if (s_a == 0)
      return MV_MultiplyTranspose (a, y, a, A, x, 0, 0, conjugate);
    else if (s_a == 1)
      return MV_MultiplyTranspose (a, y, a, A, x, 0, 1, conjugate);
    else if (s_a == static_cast<typename DomainVector::const_value_type> (-1))
      return MV_MultiplyTranspose (a, y, a, A, x, 0, -1, conjugate);
    else {
      a = aVector("a", numVecs);
      typename aVector::HostMirror h_a = Kokkos::create_mirror_view (a);
      for (int i = 0; i < numVecs; ++i) {
        h_a(i) = s_a;
      }
      Kokkos::deep_copy (a, h_a);
      return MV_MultiplyTranspose (a, y, a, A, x, 0, 2, conjugate);
    }
  } else if (s_b == 1) {
    if (s_a == 0)
      return MV_MultiplyTranspose (a, y, a, A, x, 1, 0, conjugate);
    else if (s_a == 1)
      return MV_MultiplyTranspose (a, y, a, A, x, 1, 1, conjugate);
    else if (s_a == static_cast<typename DomainVector::const_value_type> (-1))
      return MV_MultiplyTranspose (a, y, a, A, x, 1, -1, conjugate);
    else {
      a = aVector("a", numVecs);
      typename aVector::HostMirror h_a = Kokkos::create_mirror_view (a);
      for (int i = 0; i < numVecs; ++i) {
        h_a(i) = s_a;
      }
      Kokkos::deep_copy (a, h_a);
      return MV_MultiplyTranspose (a, y, a, A, x, 1, 2, conjugate);
    }
  } else if (s_b == static_cast<typename RangeVector::const_value_type> (-1)) {
    if (s_a == 0)
      return MV_MultiplyTranspose (a, y, a, A, x, -1, 0, conjugate);
    else if (s_a == 1)
      return MV_MultiplyTranspose (a, y, a, A, x, -1, 1, conjugate);
    else if (s_a == static_cast<typename DomainVector::const_value_type> (-1))
      return MV_MultiplyTranspose (a, y, a, A, x, -1, -1, conjugate);
    else {
      a = aVector("a", numVecs);
      typename aVector::HostMirror h_a = Kokkos::create_mirror_view (a);
      for (int i = 0; i < numVecs; ++i) {
        h_a(i) = s_a;
      }
      Kokkos::deep_copy (a, h_a);
      return MV_MultiplyTranspose (a, y, a, A, x, -1, 2, conjugate);
    }
  } else {
    b = aVector("b", numVecs);
    typename aVector::HostMirror h_b = Kokkos::create_mirror_view (b);
    for (int i = 0; i < numVecs; ++i) {
      h_b(i) = s_b;
    }
    Kokkos::deep_copy (b, h_b);

    if (s_a == 0)
      return MV_MultiplyTranspose (b, y, a, A, x, 2, 0, conjugate);
    else if (s_a == 1)
      return MV_MultiplyTranspose (b, y, a, A, x, 2, 1, conjugate);
    else if (s_a == static_cast<typename DomainVector::const_value_type> (-1))
      return MV_MultiplyTranspose (b, y, a, A, x, 2, -1, conjugate);
    else {
      a = aVector("a", numVecs);
      typename aVector::HostMirror h_a = Kokkos::create_mirror_view (a);
      for (int i = 0; i < numVecs; ++i) {
        h_a(i) = s_a;
      }
      Kokkos::deep_copy (a, h_a);
      return MV_MultiplyTranspose (b, y, a, A, x, 2, 2, conjugate);
    }
  }
}

  template<class RangeVector,
           class TCrsMatrix,
           class DomainVector,
           class CoeffVector1,
           class CoeffVector2,
           int doalpha,
           int dobeta>
  void
  MV_MultiplySingle (typename Kokkos::Impl::enable_if<DomainVector::Rank == 1, const CoeffVector1>::type& betav,
               const RangeVector &y,
               const CoeffVector2 &alphav,
               const TCrsMatrix& A,
               const DomainVector& x)
  {
    typedef typename TCrsMatrix::ordinal_type ordinal_type;

    if (A.numRows () <= static_cast<ordinal_type> (0)) {
      return;
    }
    if (doalpha == 0) {
      if (dobeta==2) {
        V_MulScalar (y, betav, y);
      }
      else {
        V_MulScalar (y, typename RangeVector::value_type (dobeta), y);
      }
      return;
    } else {
      typedef View< typename RangeVector::non_const_data_type ,
                    typename RangeVector::array_layout ,
                    typename RangeVector::execution_space ,
                    typename RangeVector::memory_traits >
      RangeVectorType;

      typedef View< typename DomainVector::const_data_type ,
                    typename DomainVector::array_layout ,
                    typename DomainVector::execution_space ,
                    //typename DomainVector::memory_traits >
                    Kokkos::MemoryRandomAccess >
      DomainVectorType;

      typedef View< typename CoeffVector1::const_data_type ,
                    typename CoeffVector1::array_layout ,
                    typename CoeffVector1::execution_space ,
                    Kokkos::MemoryRandomAccess >
      CoeffVector1Type;

      typedef View< typename CoeffVector2::const_data_type ,
                    typename CoeffVector2::array_layout ,
                    typename CoeffVector2::execution_space ,
                    Kokkos::MemoryRandomAccess >
      CoeffVector2Type;

      typedef CrsMatrix<typename TCrsMatrix::const_value_type,
                        typename TCrsMatrix::ordinal_type,
                        typename TCrsMatrix::execution_space,
                        typename TCrsMatrix::memory_traits,
                        typename TCrsMatrix::size_type>
      CrsMatrixType;
      typedef typename CrsMatrixType::size_type size_type;

      Impl::MV_Multiply_Check_Compatibility(betav,y,alphav,A,x,doalpha,dobeta);

      // NNZPerRow could be anywhere from 0, to A.numRows()*A.numCols().
      // Thus, the appropriate type is size_type.
      const size_type NNZPerRow = A.nnz () / A.numRows ();

      int vector_length = 1;
      while( (static_cast<size_type> (vector_length*2*3) <= NNZPerRow) && (vector_length<32) ) vector_length*=2;

#ifndef KOKKOS_FAST_COMPILE // This uses templated fucntions on doalpha and dobeta and will produce 16 kernels

      typedef MV_MultiplySingleFunctor<RangeVectorType, CrsMatrixType, DomainVectorType,
                                       CoeffVector1Type, CoeffVector2Type, doalpha, dobeta > OpType ;

      const typename CrsMatrixType::ordinal_type nrow = A.numRows();

      OpType op(betav,alphav,A,x,y,RowsPerThread<typename RangeVector::execution_space >(NNZPerRow)) ;

      const int rows_per_thread = RowsPerThread<typename RangeVector::execution_space >(NNZPerRow);
      const int team_size = Kokkos::TeamPolicy< typename RangeVector::execution_space >::team_size_recommended(op,vector_length);
      const int rows_per_team = rows_per_thread * team_size;
      const int nteams = (nrow+rows_per_team-1)/rows_per_team;
      Kokkos::parallel_for( Kokkos::TeamPolicy< typename RangeVector::execution_space >
         ( nteams , team_size , vector_length ) , op );
#else // NOT KOKKOS_FAST_COMPILE

      typedef MV_MultiplySingleFunctor<RangeVectorType, CrsMatrixType, DomainVectorType,
                               CoeffVector1Type, CoeffVector2Type, 2, 2> OpType ;

      int numVecs = x.dimension_1(); // == 1
      CoeffVector1 beta = betav;
      CoeffVector2 alpha = alphav;

      if(doalpha!=2) {
              alpha = CoeffVector2("CrsMatrix::auto_a", numVecs);
              typename CoeffVector2::HostMirror h_a = Kokkos::create_mirror_view(alpha);
              typename CoeffVector2::value_type s_a = (typename CoeffVector2::value_type) doalpha;

              for(int i = 0; i < numVecs; i++)
                h_a(i) = s_a;

              Kokkos::deep_copy(alpha, h_a);
      }
      if(dobeta!=2) {
              beta = CoeffVector1("CrsMatrix::auto_b", numVecs);
              typename CoeffVector1::HostMirror h_b = Kokkos::create_mirror_view(beta);
              typename CoeffVector1::value_type s_b = (typename CoeffVector1::value_type) dobeta;

              for(int i = 0; i < numVecs; i++)
                h_b(i) = s_b;

              Kokkos::deep_copy(beta, h_b);
      }

      const typename CrsMatrixType::ordinal_type nrow = A.numRows();

      OpType op(beta,alpha,A,x,y,RowsPerThread<typename RangeVector::execution_space >(NNZPerRow)) ;
      const int team_size = Kokkos::TeamPolicy< typename RangeVector::execution_space >::team_size_recommended(op,vector_length);
      const int rows_per_team = rows_per_thread * team_size;
      const int nteams = (nrow+rows_per_team-1)/rows_per_team;
      Kokkos::parallel_for( Kokkos::TeamPolicy< typename RangeVector::execution_space >
         ( nteams , team_size , vector_length ) , op );

#endif // KOKKOS_FAST_COMPILE
    }
  }
} // namespace Kokkos

#endif /* KOKKOS_CRSMATRIX_H_ */
