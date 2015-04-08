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
#ifndef KOKKOS_BLAS1_MV_IMPL_HPP_
#define KOKKOS_BLAS1_MV_IMPL_HPP_

#include <TpetraKernels_config.h>
#include <Kokkos_Core.hpp>
#include <Kokkos_InnerProductSpaceTraits.hpp>

namespace KokkosBlas {
namespace Impl {

/// \brief Dot product functor for single vectors.
///
/// \tparam RV 0-D output View
/// \tparam XV 1-D input View
/// \tparam YV 1-D input View
/// \tparam SizeType Index type.  Use int (32 bits) if possible.
template<class RV, class XV, class YV, class SizeType = typename XV::size_type>
struct V_Dot_Functor
{
  typedef typename XV::execution_space              execution_space;
  typedef SizeType                                        size_type;
  typedef typename XV::non_const_value_type             xvalue_type;
  typedef Kokkos::Details::InnerProductSpaceTraits<xvalue_type> IPT;
  typedef Kokkos::Details::ArithTraits<typename IPT::dot_type>   AT;
  typedef typename IPT::dot_type                         value_type;

  RV m_r;
  XV m_x;
  YV m_y;

  V_Dot_Functor (const RV& r, const XV& x, const YV& y) :
    m_r (r), m_x (x), m_y (y)
  {}

  KOKKOS_INLINE_FUNCTION
  void operator() (const size_type& i, value_type& sum) const
  {
    sum += IPT::dot (m_x(i), m_y(i)); // m_x(i) * m_y(i)
  }

  KOKKOS_INLINE_FUNCTION void init (value_type& update) const
  {
    update = AT::zero ();
  }

  KOKKOS_INLINE_FUNCTION void
  join (volatile value_type& update,
        const volatile value_type& source) const
  {
    update += source;
  }

  // On device, write the reduction result to the output View.
  KOKKOS_INLINE_FUNCTION void
  final (const value_type& dst) const
  {
    m_r() = dst;
  }
};

/// \brief Column-wise dot product functor for multivectors; works for
///   any layout, but best performance with LayoutRight.
///
/// \tparam RV 1-D output View
/// \tparam XMV 2-D input View
/// \tparam YMV 2-D input View
/// \tparam SizeType Index type.  Use int (32 bits) if possible.
template<class RV, class XMV, class YMV, class SizeType = typename XMV::size_type>
struct MV_Dot_Right_FunctorVector
{
  typedef typename XMV::execution_space             execution_space;
  typedef SizeType                                        size_type;
  typedef typename XMV::non_const_value_type            xvalue_type;
  typedef Kokkos::Details::InnerProductSpaceTraits<xvalue_type> IPT;
  typedef Kokkos::Details::ArithTraits<typename IPT::dot_type>   AT;
  typedef typename IPT::dot_type                       value_type[];

  size_type value_count;
  RV m_r;
  typename XMV::const_type m_x;
  typename YMV::const_type m_y;

  MV_Dot_Right_FunctorVector (const RV& r, const XMV& x, const YMV& y) :
    value_count (x.dimension_1 ()), m_r (r), m_x (x), m_y (y)
  {}

  KOKKOS_INLINE_FUNCTION
  void operator() (const size_type& i, value_type sum) const
  {
    const size_type numVecs = value_count;
#ifdef KOKKOS_HAVE_PRAGMA_IVDEP
#pragma ivdep
#endif
#ifdef KOKKOS_HAVE_PRAGMA_VECTOR
#pragma vector always
#endif
    for (size_type k = 0; k < numVecs; ++k) {
      sum[k] += IPT::dot (m_x(i,k), m_y(i,k)); // m_x(i,k) * m_y(i,k)
    }
  }

  KOKKOS_INLINE_FUNCTION void init (value_type update) const
  {
    const size_type numVecs = value_count;
#ifdef KOKKOS_HAVE_PRAGMA_IVDEP
#pragma ivdep
#endif
#ifdef KOKKOS_HAVE_PRAGMA_VECTOR
#pragma vector always
#endif
    for (size_type k = 0; k < numVecs; ++k) {
      update[k] = AT::zero ();
    }
  }

  KOKKOS_INLINE_FUNCTION void
  join (volatile value_type update,
        const volatile value_type source) const
  {
    const size_type numVecs = value_count;
#ifdef KOKKOS_HAVE_PRAGMA_IVDEP
#pragma ivdep
#endif
#ifdef KOKKOS_HAVE_PRAGMA_VECTOR
#pragma vector always
#endif
    for (size_type k = 0; k < numVecs; ++k) {
      update[k] += source[k];
    }
  }

  // On device, write the reduction result to the output View.
  KOKKOS_INLINE_FUNCTION void
  final (const value_type dst) const
  {
    const size_type numVecs = value_count;
#ifdef KOKKOS_HAVE_PRAGMA_IVDEP
#pragma ivdep
#endif
#ifdef KOKKOS_HAVE_PRAGMA_VECTOR
#pragma vector always
#endif
    for (size_type k = 0; k < numVecs; ++k) {
      m_r(k) = dst[k];
    }
  }
};

/// \brief Column-wise dot product functor for multivectors with
///   number of columns known at compile time; works for any layout,
///   but best performance with LayoutRight.
///
/// \tparam RV 1-D output View
/// \tparam XMV 2-D input View
/// \tparam YMV 2-D input View
/// \tparam UNROLL Number of columns (vectors)
/// \tparam SizeType Index type.  Use int (32 bits) if possible.
template<class RV, class XMV, class YMV, int UNROLL, class SizeType = typename XMV::size_type>
struct MV_Dot_Right_FunctorUnroll
{
  typedef typename XMV::execution_space             execution_space;
  typedef SizeType                                        size_type;
  typedef typename XMV::non_const_value_type            xvalue_type;
  typedef Kokkos::Details::InnerProductSpaceTraits<xvalue_type> IPT;
  typedef Kokkos::Details::ArithTraits<typename IPT::dot_type>   AT;
  typedef typename IPT::dot_type                       value_type[];

  size_type value_count;
  RV m_r;
  typename XMV::const_type m_x;
  typename YMV::const_type m_y;

  MV_Dot_Right_FunctorUnroll (const RV& r, const XMV& x, const YMV& y) :
    value_count (x.dimension_1 ()), m_r (r), m_x (x), m_y (y)
  {}

  KOKKOS_INLINE_FUNCTION
  void operator() (const size_type& i, value_type sum) const
  {
#ifdef KOKKOS_HAVE_PRAGMA_UNROLL
#pragma unroll
#endif
    for (int k = 0; k < UNROLL; ++k) {
      sum[k] += IPT::dot (m_x(i,k), m_y(i,k)); // m_x(i,k) * m_y(i,k)
    }
  }

  KOKKOS_INLINE_FUNCTION void init (volatile value_type update) const
  {
#ifdef KOKKOS_HAVE_PRAGMA_UNROLL
#pragma unroll
#endif
    for (int k = 0; k < UNROLL; ++k) {
      update[k] = AT::zero ();
    }
  }

  KOKKOS_INLINE_FUNCTION void
  join (volatile value_type update,
        const volatile value_type source) const
  {
#ifdef KOKKOS_HAVE_PRAGMA_UNROLL
#pragma unroll
#endif
    for (int k = 0; k < UNROLL; ++k) {
      update[k] += source[k];
    }
  }

  // On device, write the reduction result to the output View.
  KOKKOS_INLINE_FUNCTION void
  final (const value_type dst) const
  {
#ifdef KOKKOS_HAVE_PRAGMA_UNROLL
#pragma unroll
#endif
    for (int k = 0; k < UNROLL; ++k) {
      m_r(k) = dst[k];
    }
  }
};


//! Implementation of KokkosBlas::dot for multivectors.
template<class RT, class RL, class RD, class RM, class RS,
         class XT, class XL, class XD, class XM, class XS,
         class YT, class YL, class YD, class YM, class YS>
struct Dot_MV {
  typedef Kokkos::View<RT,RL,RD,RM,RS> RV;
  typedef Kokkos::View<XT,XL,XD,XM,XS> XMV;
  typedef Kokkos::View<YT,YL,YD,YM,YS> YMV;
  typedef typename XMV::size_type size_type;

  /// \brief Compute the dot product(s) of the column(s) of the
  ///   multivectors (2-D views) x and y, and store result(s) in r.
  static void dot (const RV& r, const XMV& X, const YMV& Y)
  {
    const size_type numRows = X.dimension_0 ();
    const size_type numVecs = X.dimension_1 ();

    if (numVecs > 16) {
      MV_Dot_Right_FunctorVector<RV, XMV, YMV> op (r, X, Y);
      Kokkos::parallel_reduce (numRows, op);
    }
    else {
      switch (numVecs) {
      case 16: {
        MV_Dot_Right_FunctorUnroll<RV, XMV, YMV, 16> op (r, X, Y);
        Kokkos::parallel_reduce (numRows, op);
        break;
      }
      case 15: {
        MV_Dot_Right_FunctorUnroll<RV, XMV, YMV, 15> op (r, X, Y);
        Kokkos::parallel_reduce (numRows, op);
        break;
      }
      case 14: {
        MV_Dot_Right_FunctorUnroll<RV, XMV, YMV, 14> op (r, X, Y);
        Kokkos::parallel_reduce (numRows, op);
        break;
      }
      case 13: {
        MV_Dot_Right_FunctorUnroll<RV, XMV, YMV, 13> op (r, X, Y);
        Kokkos::parallel_reduce (numRows, op);
        break;
      }
      case 12: {
        MV_Dot_Right_FunctorUnroll<RV, XMV, YMV, 12> op (r, X, Y);
        Kokkos::parallel_reduce (numRows, op);
        break;
      }
      case 11: {
        MV_Dot_Right_FunctorUnroll<RV, XMV, YMV, 11> op (r, X, Y);
        Kokkos::parallel_reduce (numRows, op);
        break;
      }
      case 10: {
        MV_Dot_Right_FunctorUnroll<RV, XMV, YMV, 10> op (r, X, Y);
        Kokkos::parallel_reduce (numRows, op);
        break;
      }
      case 9: {
        MV_Dot_Right_FunctorUnroll<RV, XMV, YMV, 9> op (r, X, Y);
        Kokkos::parallel_reduce (numRows, op);
        break;
      }
      case 8: {
        MV_Dot_Right_FunctorUnroll<RV, XMV, YMV, 8> op (r, X, Y);
        Kokkos::parallel_reduce (numRows, op);
        break;
      }
      case 7: {
        MV_Dot_Right_FunctorUnroll<RV, XMV, YMV, 7> op (r, X, Y);
        Kokkos::parallel_reduce (numRows, op);
        break;
      }
      case 6: {
        MV_Dot_Right_FunctorUnroll<RV, XMV, YMV, 6> op (r, X, Y);
        Kokkos::parallel_reduce (numRows, op);
        break;
      }
      case 5: {
        MV_Dot_Right_FunctorUnroll<RV, XMV, YMV, 5> op (r, X, Y);
        Kokkos::parallel_reduce (numRows, op);
        break;
      }
      case 4: {
        MV_Dot_Right_FunctorUnroll<RV, XMV, YMV, 4> op (r, X, Y);
        Kokkos::parallel_reduce (numRows, op);
        break;
      }
      case 3: {
        MV_Dot_Right_FunctorUnroll<RV, XMV, YMV, 3> op (r, X, Y);
        Kokkos::parallel_reduce (numRows, op);
        break;
      }
      case 2: {
        MV_Dot_Right_FunctorUnroll<RV, XMV, YMV, 2> op (r, X, Y);
        Kokkos::parallel_reduce (numRows, op);
        break;
      }
      case 1: {
        // RV needs to turn 0-D, and XMV and YMV need to turn 1-D.
        using Kokkos::ALL;
        using Kokkos::subview;
        typedef Kokkos::View<typename RV::value_type, typename RV::array_layout,
            typename RV::device_type, typename RV::memory_traits,
            typename RV::specialize> RV1D;
        typedef Kokkos::View<typename XMV::const_value_type*,
            XL,XD,XM,XS> XMV1D;
        typedef Kokkos::View<typename YMV::const_value_type*,
            YL,YD,YM,YS> YMV1D;
        typedef V_Dot_Functor<RV1D, XMV1D, YMV1D> op_type;
        op_type op (subview (r, 0), subview (X, ALL (), 0), subview (Y, ALL (), 0));
        Kokkos::parallel_reduce (numRows, op);
        break;
      }
      } // switch
    } // if-else
  }

  /// \brief Compute the dot product of X(:,X_col) and Y(:,Y_col), and
  ///   store result in r(r_col).
  static void
  dot (const RV& r, const size_t r_col,
       const XMV& X, const size_t X_col,
       const YMV& Y, const size_t Y_col)
  {
    // RV needs to turn 0-D, and XMV and YMV need to turn 1-D.
    using Kokkos::ALL;
    using Kokkos::subview;
    typedef Kokkos::View<typename RV::value_type, typename RV::array_layout,
      typename RV::device_type, typename RV::memory_traits,
      typename RV::specialize> RV0D;
    typedef Kokkos::View<typename XMV::const_value_type*, XL, XD, XM, XS> XMV1D;
    typedef Kokkos::View<typename YMV::const_value_type*, YL, YD, YM, YS> YMV1D;

    const size_type numRows = X.dimension_0 ();
    const size_type numCols = X.dimension_1 ();
    if (numRows < static_cast<size_type> (INT_MAX) &&
        numRows * numCols < static_cast<size_type> (INT_MAX)) {
      typedef V_Dot_Functor<RV0D, XMV1D, YMV1D, int> op_type;
      op_type op (subview (r, r_col), subview (X, ALL (), X_col),
                  subview (Y, ALL (), Y_col));
      Kokkos::parallel_reduce (numRows, op);
    }
    else {
      typedef V_Dot_Functor<RV0D, XMV1D, YMV1D, size_type> op_type;
      op_type op (subview (r, r_col), subview (X, ALL (), X_col),
                  subview (Y, ALL (), Y_col));
      Kokkos::parallel_reduce (numRows, op);
    }
  }
};

// Full specializations for cases of interest for Tpetra::MultiVector.
//
// Currently, we include specializations for Scalar = double,
// LayoutLeft (which is what Tpetra::MultiVector uses at the moment),
// and all execution spaces.  This may change in the future.

#ifdef KOKKOS_HAVE_SERIAL
template<>
struct Dot_MV<double*,
              Kokkos::Serial::array_layout,
              Kokkos::Device<Kokkos::Serial, Kokkos::HostSpace>,
              Kokkos::MemoryTraits<Kokkos::Unmanaged>,
              Kokkos::Impl::ViewDefault,
              const double**,
              Kokkos::LayoutLeft,
              Kokkos::Device<Kokkos::Serial, Kokkos::HostSpace>,
              Kokkos::MemoryTraits<Kokkos::Unmanaged>,
              Kokkos::Impl::ViewDefault,
              const double**,
              Kokkos::LayoutLeft,
              Kokkos::Device<Kokkos::Serial, Kokkos::HostSpace>,
              Kokkos::MemoryTraits<Kokkos::Unmanaged>,
              Kokkos::Impl::ViewDefault> {
  typedef double* RT;
  typedef Kokkos::Serial::array_layout RL;
  typedef Kokkos::Device<Kokkos::Serial, Kokkos::HostSpace> RD;
  typedef Kokkos::MemoryTraits<Kokkos::Unmanaged> RM;
  typedef Kokkos::Impl::ViewDefault RS;
  typedef Kokkos::View<RT,RL,RD,RM,RS> RV;

  typedef const double** XT;
  typedef Kokkos::LayoutLeft XL;
  typedef Kokkos::Device<Kokkos::Serial, Kokkos::HostSpace> XD;
  typedef Kokkos::MemoryTraits<Kokkos::Unmanaged> XM;
  typedef Kokkos::Impl::ViewDefault XS;
  typedef Kokkos::View<XT,XL,XD,XM,XS> XMV;

  typedef Kokkos::View<XT,XL,XD,XM,XS> YMV;
  typedef XMV::size_type size_type;

  static void dot (const RV& r, const XMV& X, const YMV& Y);

  static void
  dot (const RV& r, const size_t r_col,
       const XMV& X, const size_t X_col,
       const YMV& Y, const size_t Y_col);
};
#endif // KOKKOS_HAVE_SERIAL

#ifdef KOKKOS_HAVE_OPENMP
template<>
struct Dot_MV<double*,
              Kokkos::OpenMP::array_layout,
              Kokkos::Device<Kokkos::OpenMP, Kokkos::HostSpace>,
              Kokkos::MemoryTraits<Kokkos::Unmanaged>,
              Kokkos::Impl::ViewDefault,
              const double**,
              Kokkos::LayoutLeft,
              Kokkos::Device<Kokkos::OpenMP, Kokkos::HostSpace>,
              Kokkos::MemoryTraits<Kokkos::Unmanaged>,
              Kokkos::Impl::ViewDefault,
              const double**,
              Kokkos::LayoutLeft,
              Kokkos::Device<Kokkos::OpenMP, Kokkos::HostSpace>,
              Kokkos::MemoryTraits<Kokkos::Unmanaged>,
              Kokkos::Impl::ViewDefault> {
  typedef double* RT;
  typedef Kokkos::OpenMP::array_layout RL;
  typedef Kokkos::Device<Kokkos::OpenMP, Kokkos::HostSpace> RD;
  typedef Kokkos::MemoryTraits<Kokkos::Unmanaged> RM;
  typedef Kokkos::Impl::ViewDefault RS;
  typedef Kokkos::View<RT,RL,RD,RM,RS> RV;

  typedef const double** XT;
  typedef Kokkos::LayoutLeft XL;
  typedef Kokkos::Device<Kokkos::OpenMP, Kokkos::HostSpace> XD;
  typedef Kokkos::MemoryTraits<Kokkos::Unmanaged> XM;
  typedef Kokkos::Impl::ViewDefault XS;
  typedef Kokkos::View<XT,XL,XD,XM,XS> XMV;

  typedef Kokkos::View<XT,XL,XD,XM,XS> YMV;
  typedef XMV::size_type size_type;

  static void dot (const RV& r, const XMV& X, const YMV& Y);

  static void
  dot (const RV& r, const size_t r_col,
       const XMV& X, const size_t X_col,
       const YMV& Y, const size_t Y_col);
};
#endif // KOKKOS_HAVE_OPENMP

#ifdef KOKKOS_HAVE_PTHREAD
template<>
struct Dot_MV<double*,
              Kokkos::Threads::array_layout,
              Kokkos::Device<Kokkos::Threads, Kokkos::HostSpace>,
              Kokkos::MemoryTraits<Kokkos::Unmanaged>,
              Kokkos::Impl::ViewDefault,
              const double**,
              Kokkos::LayoutLeft,
              Kokkos::Device<Kokkos::Threads, Kokkos::HostSpace>,
              Kokkos::MemoryTraits<Kokkos::Unmanaged>,
              Kokkos::Impl::ViewDefault,
              const double**,
              Kokkos::LayoutLeft,
              Kokkos::Device<Kokkos::Threads, Kokkos::HostSpace>,
              Kokkos::MemoryTraits<Kokkos::Unmanaged>,
              Kokkos::Impl::ViewDefault> {
  typedef double* RT;
  typedef Kokkos::Threads::array_layout RL;
  typedef Kokkos::Device<Kokkos::Threads, Kokkos::HostSpace> RD;
  typedef Kokkos::MemoryTraits<Kokkos::Unmanaged> RM;
  typedef Kokkos::Impl::ViewDefault RS;
  typedef Kokkos::View<RT,RL,RD,RM,RS> RV;

  typedef const double** XT;
  typedef Kokkos::LayoutLeft XL;
  typedef Kokkos::Device<Kokkos::Threads, Kokkos::HostSpace> XD;
  typedef Kokkos::MemoryTraits<Kokkos::Unmanaged> XM;
  typedef Kokkos::Impl::ViewDefault XS;
  typedef Kokkos::View<XT,XL,XD,XM,XS> XMV;

  typedef Kokkos::View<XT,XL,XD,XM,XS> YMV;
  typedef XMV::size_type size_type;

  static void dot (const RV& r, const XMV& X, const YMV& Y);

  static void
  dot (const RV& r, const size_t r_col,
       const XMV& X, const size_t X_col,
       const YMV& Y, const size_t Y_col);
};
#endif // KOKKOS_HAVE_PTHREAD

#ifdef KOKKOS_HAVE_CUDA
template<>
struct Dot_MV<double*,
              Kokkos::Cuda::array_layout,
              Kokkos::Device<Kokkos::Cuda, Kokkos::CudaSpace>,
              Kokkos::MemoryTraits<Kokkos::Unmanaged>,
              Kokkos::Impl::ViewDefault,
              const double**,
              Kokkos::LayoutLeft,
              Kokkos::Device<Kokkos::Cuda, Kokkos::CudaSpace>,
              Kokkos::MemoryTraits<Kokkos::Unmanaged>,
              Kokkos::Impl::ViewDefault,
              const double**,
              Kokkos::LayoutLeft,
              Kokkos::Device<Kokkos::Cuda, Kokkos::CudaSpace>,
              Kokkos::MemoryTraits<Kokkos::Unmanaged>,
              Kokkos::Impl::ViewDefault> {
  typedef double* RT;
  typedef Kokkos::Cuda::array_layout RL;
  typedef Kokkos::Device<Kokkos::Cuda, Kokkos::CudaSpace> RD;
  typedef Kokkos::MemoryTraits<Kokkos::Unmanaged> RM;
  typedef Kokkos::Impl::ViewDefault RS;
  typedef Kokkos::View<RT,RL,RD,RM,RS> RV;

  typedef const double** XT;
  typedef Kokkos::LayoutLeft XL;
  typedef Kokkos::Device<Kokkos::Cuda, Kokkos::CudaSpace> XD;
  typedef Kokkos::MemoryTraits<Kokkos::Unmanaged> XM;
  typedef Kokkos::Impl::ViewDefault XS;
  typedef Kokkos::View<XT,XL,XD,XM,XS> XMV;

  typedef Kokkos::View<XT,XL,XD,XM,XS> YMV;
  typedef XMV::size_type size_type;

  static void dot (const RV& r, const XMV& X, const YMV& Y);

  static void
  dot (const RV& r, const size_t r_col,
       const XMV& X, const size_t X_col,
       const YMV& Y, const size_t Y_col);
};
#endif // KOKKOS_HAVE_CUDA

#ifdef KOKKOS_HAVE_CUDA
template<>
struct Dot_MV<double*,
              Kokkos::Cuda::array_layout,
              Kokkos::Device<Kokkos::Cuda, Kokkos::CudaUVMSpace>,
              Kokkos::MemoryTraits<Kokkos::Unmanaged>,
              Kokkos::Impl::ViewDefault,
              const double**,
              Kokkos::LayoutLeft,
              Kokkos::Device<Kokkos::Cuda, Kokkos::CudaUVMSpace>,
              Kokkos::MemoryTraits<Kokkos::Unmanaged>,
              Kokkos::Impl::ViewDefault,
              const double**,
              Kokkos::LayoutLeft,
              Kokkos::Device<Kokkos::Cuda, Kokkos::CudaUVMSpace>,
              Kokkos::MemoryTraits<Kokkos::Unmanaged>,
              Kokkos::Impl::ViewDefault> {
  typedef double* RT;
  typedef Kokkos::Cuda::array_layout RL;
  typedef Kokkos::Device<Kokkos::Cuda, Kokkos::CudaUVMSpace> RD;
  typedef Kokkos::MemoryTraits<Kokkos::Unmanaged> RM;
  typedef Kokkos::Impl::ViewDefault RS;
  typedef Kokkos::View<RT,RL,RD,RM,RS> RV;

  typedef const double** XT;
  typedef Kokkos::LayoutLeft XL;
  typedef Kokkos::Device<Kokkos::Cuda, Kokkos::CudaUVMSpace> XD;
  typedef Kokkos::MemoryTraits<Kokkos::Unmanaged> XM;
  typedef Kokkos::Impl::ViewDefault XS;
  typedef Kokkos::View<XT,XL,XD,XM,XS> XMV;

  typedef Kokkos::View<XT,XL,XD,XM,XS> YMV;
  typedef XMV::size_type size_type;

  static void dot (const RV& r, const XMV& X, const YMV& Y);

  static void
  dot (const RV& r, const size_t r_col,
       const XMV& X, const size_t X_col,
       const YMV& Y, const size_t Y_col);
};
#endif // KOKKOS_HAVE_CUDA

//
// fill
//

template<class XMV, class SizeType = typename XMV::size_type>
struct MV_FillFunctor {
  typedef typename XMV::execution_space             execution_space;
  typedef SizeType                                        size_type;
  typedef typename XMV::non_const_value_type            xvalue_type;

  const size_type numCols_;
  const xvalue_type val_;
  XMV X_;

  MV_FillFunctor (const XMV& X, const xvalue_type& val) :
    numCols_ (X.dimension_1 ()), val_ (val), X_ (X)
  {}

  KOKKOS_INLINE_FUNCTION
  void operator() (const size_type& i) const
  {
    for (size_type j = 0; j < numCols_; ++j) {
      X_(i,j) = val_;
    }
  }
};

template<class XV, class SizeType = typename XV::size_type>
struct V_FillFunctor {
  typedef typename XV::execution_space            execution_space;
  typedef SizeType                                      size_type;
  typedef typename XV::non_const_value_type           xvalue_type;

  const xvalue_type val_;
  XV x_;

  V_FillFunctor (const XV& x, const xvalue_type& val) :
    val_ (val), x_ (x)
  {}

  KOKKOS_INLINE_FUNCTION
  void operator() (const size_type& i) const
  {
    x_(i) = val_;
  }
};

//! Implementation of KokkosBlas::fill for multivectors.
template<class XT, class XL, class XD, class XM, class XS>
struct Fill_MV {
  typedef Kokkos::View<XT,XL,XD,XM,XS> XMV;
  typedef typename XMV::execution_space execution_space;
  typedef typename XMV::size_type size_type;

  static void fill (const XMV& X, const typename XMV::non_const_value_type& val)
  {
    const size_type numRows = X.dimension_0 ();
    const size_type numCols = X.dimension_1 ();

    // The first condition helps avoid overflow with the
    // multiplication in the second condition.
    if (numRows < static_cast<size_type> (INT_MAX) &&
        numRows * numCols < static_cast<size_type> (INT_MAX)) {
      Kokkos::RangePolicy<execution_space, int> policy (0, numRows);
      MV_FillFunctor<XMV, int> op (X, val);
      Kokkos::parallel_for (policy, op);
    }
    else {
      Kokkos::RangePolicy<execution_space, size_type> policy (0, numRows);
      MV_FillFunctor<XMV, size_type> op (X, val);
      Kokkos::parallel_for (policy, op);
    }
  }
};

#ifdef KOKKOS_HAVE_SERIAL
template<>
struct Fill_MV<double**, Kokkos::LayoutLeft, Kokkos::Device<Kokkos::Serial, Kokkos::HostSpace>,
               Kokkos::MemoryTraits<Kokkos::Unmanaged>, Kokkos::Impl::ViewDefault> {
  typedef Kokkos::View<double**, Kokkos::LayoutLeft, Kokkos::Device<Kokkos::Serial, Kokkos::HostSpace>,
                       Kokkos::MemoryTraits<Kokkos::Unmanaged>, Kokkos::Impl::ViewDefault> XMV;
  typedef Kokkos::Serial execution_space;
  typedef XMV::size_type size_type;

  static void fill (const XMV& X, const double& val);
};
#endif // KOKKOS_HAVE_SERIAL

#ifdef KOKKOS_HAVE_OPENMP
template<>
struct Fill_MV<double**, Kokkos::LayoutLeft, Kokkos::Device<Kokkos::OpenMP, Kokkos::HostSpace>,
               Kokkos::MemoryTraits<Kokkos::Unmanaged>, Kokkos::Impl::ViewDefault> {
  typedef Kokkos::View<double**, Kokkos::LayoutLeft, Kokkos::Device<Kokkos::OpenMP, Kokkos::HostSpace>,
                       Kokkos::MemoryTraits<Kokkos::Unmanaged>, Kokkos::Impl::ViewDefault> XMV;
  typedef Kokkos::OpenMP execution_space;
  typedef XMV::size_type size_type;

  static void fill (const XMV& X, const double& val);
};
#endif // KOKKOS_HAVE_OPENMP

#ifdef KOKKOS_HAVE_PTHREAD
template<>
struct Fill_MV<double**, Kokkos::LayoutLeft, Kokkos::Device<Kokkos::Threads, Kokkos::HostSpace>,
               Kokkos::MemoryTraits<Kokkos::Unmanaged>, Kokkos::Impl::ViewDefault> {
  typedef Kokkos::View<double**, Kokkos::LayoutLeft, Kokkos::Device<Kokkos::Threads, Kokkos::HostSpace>,
                       Kokkos::MemoryTraits<Kokkos::Unmanaged>, Kokkos::Impl::ViewDefault> XMV;
  typedef Kokkos::Threads execution_space;
  typedef XMV::size_type size_type;

  static void fill (const XMV& X, const double& val);
};
#endif // KOKKOS_HAVE_PTHREAD

#ifdef KOKKOS_HAVE_CUDA
template<>
struct Fill_MV<double**, Kokkos::LayoutLeft, Kokkos::Device<Kokkos::Cuda, Kokkos::CudaSpace>,
               Kokkos::MemoryTraits<Kokkos::Unmanaged>, Kokkos::Impl::ViewDefault> {
  typedef Kokkos::View<double**, Kokkos::LayoutLeft, Kokkos::Device<Kokkos::Cuda, Kokkos::CudaSpace>,
                       Kokkos::MemoryTraits<Kokkos::Unmanaged>, Kokkos::Impl::ViewDefault> XMV;
  typedef Kokkos::Cuda execution_space;
  typedef XMV::size_type size_type;

  static void fill (const XMV& X, const double& val);
};
#endif // KOKKOS_HAVE_CUDA

#ifdef KOKKOS_HAVE_CUDA
template<>
struct Fill_MV<double**, Kokkos::LayoutLeft, Kokkos::Device<Kokkos::Cuda, Kokkos::CudaUVMSpace>,
               Kokkos::MemoryTraits<Kokkos::Unmanaged>, Kokkos::Impl::ViewDefault> {
  typedef Kokkos::View<double**, Kokkos::LayoutLeft, Kokkos::Device<Kokkos::Cuda, Kokkos::CudaUVMSpace>,
                       Kokkos::MemoryTraits<Kokkos::Unmanaged>, Kokkos::Impl::ViewDefault> XMV;
  typedef Kokkos::Cuda execution_space;
  typedef XMV::size_type size_type;

  static void fill (const XMV& X, const double& val);
};
#endif // KOKKOS_HAVE_CUDA

//
// nrm2_squared
//

/// \brief 2-norm (squared) functor for single vectors.
///
/// \tparam RV 0-D output View
/// \tparam XV 1-D input View
/// \tparam SizeType Index type.  Use int (32 bits) if possible.
template<class RV, class XV, class SizeType = typename XV::size_type>
struct V_Nrm2Squared_Functor
{
  typedef typename XV::execution_space              execution_space;
  typedef SizeType                                        size_type;
  typedef typename XV::non_const_value_type             xvalue_type;
  typedef Kokkos::Details::InnerProductSpaceTraits<xvalue_type> IPT;
  typedef Kokkos::Details::ArithTraits<typename IPT::mag_type>   AT;
  typedef typename IPT::mag_type                         value_type;

  RV m_r;
  typename XV::const_type m_x;

  V_Nrm2Squared_Functor (const RV& r, const XV& x) :
    m_r (r), m_x (x)
  {}

  KOKKOS_INLINE_FUNCTION
  void operator() (const size_type& i, value_type& sum) const
  {
    const typename IPT::mag_type tmp = IPT::norm (m_x(i));
    sum += tmp * tmp;
  }

  KOKKOS_INLINE_FUNCTION void init (value_type& update) const
  {
    update = AT::zero ();
  }

  KOKKOS_INLINE_FUNCTION void
  join (volatile value_type& update,
        const volatile value_type& source) const
  {
    update += source;
  }

  // On device, write the reduction result to the output View.
  KOKKOS_INLINE_FUNCTION void
  final (const value_type& dst) const
  {
    m_r() = dst;
  }
};

/// \brief Column-wise 2-norm functor for multivectors; works for
///   any layout, but best performance with LayoutRight.
///
/// \tparam RV 1-D output View
/// \tparam XMV 2-D input View
/// \tparam SizeType Index type.  Use int (32 bits) if possible.
template<class RV, class XMV, class SizeType = typename XMV::size_type>
struct MV_Nrm2Squared_Right_FunctorVector
{
  typedef typename XMV::execution_space             execution_space;
  typedef SizeType                                        size_type;
  typedef typename XMV::non_const_value_type            xvalue_type;
  typedef Kokkos::Details::InnerProductSpaceTraits<xvalue_type> IPT;
  typedef Kokkos::Details::ArithTraits<typename IPT::mag_type>   AT;
  typedef typename IPT::mag_type                       value_type[];

  size_type value_count;
  RV m_r;
  typename XMV::const_type m_x;

  MV_Nrm2Squared_Right_FunctorVector (const RV& r, const XMV& x) :
    value_count (x.dimension_1 ()), m_r (r), m_x (x)
  {}

  KOKKOS_INLINE_FUNCTION void
  operator() (const size_type i, value_type sum) const
  {
    const size_type numVecs = value_count;
#ifdef KOKKOS_HAVE_PRAGMA_IVDEP
#pragma ivdep
#endif
#ifdef KOKKOS_HAVE_PRAGMA_VECTOR
#pragma vector always
#endif
    for (size_type j = 0; j < numVecs; ++j) {
      const typename IPT::mag_type tmp = IPT::norm (m_x(i,j));
      sum[j] += tmp * tmp;
    }
  }

  KOKKOS_INLINE_FUNCTION void
  init (value_type update) const
  {
    const size_type numVecs = value_count;
#ifdef KOKKOS_HAVE_PRAGMA_IVDEP
#pragma ivdep
#endif
#ifdef KOKKOS_HAVE_PRAGMA_VECTOR
#pragma vector always
#endif
    for (size_type j = 0; j < numVecs; ++j) {
      update[j] = AT::zero ();
    }
  }

  KOKKOS_INLINE_FUNCTION void
  join (volatile value_type update,
        const volatile value_type source) const
  {
    const size_type numVecs = value_count;
#ifdef KOKKOS_HAVE_PRAGMA_IVDEP
#pragma ivdep
#endif
#ifdef KOKKOS_HAVE_PRAGMA_VECTOR
#pragma vector always
#endif
    for (size_type j = 0; j < numVecs; ++j) {
      update[j] += source[j];
    }
  }

  // On device, write the reduction result to the output View.
  KOKKOS_INLINE_FUNCTION void
  final (const value_type dst) const
  {
    const size_type numVecs = value_count;
#ifdef KOKKOS_HAVE_PRAGMA_IVDEP
#pragma ivdep
#endif
#ifdef KOKKOS_HAVE_PRAGMA_VECTOR
#pragma vector always
#endif
    for (size_type j = 0; j < numVecs; ++j) {
      m_r(j) = dst[j];
    }
  }
};

//! Implementation of KokkosBlas::nrm2_squared for multivectors.
template<class RT, class RL, class RD, class RM, class RS,
         class XT, class XL, class XD, class XM, class XS>
struct Nrm2_MV {
  typedef Kokkos::View<RT,RL,RD,RM,RS> RV;
  typedef Kokkos::View<XT,XL,XD,XM,XS> XMV;
  typedef typename XMV::execution_space execution_space;
  typedef typename XMV::size_type size_type;

  /// \brief Compute the square of the 2-norm(s) of the column(s) of
  ///   the multivector (2-D View) X, and store result(s) in r.
  static void nrm2_squared (const RV& r, const XMV& X)
  {
    const size_type numRows = X.dimension_0 ();
    const size_type numCols = X.dimension_1 ();

    // int is generally faster than size_t, but check for overflow first.
    if (numRows < static_cast<size_type> (INT_MAX) &&
        numRows * numCols < static_cast<size_type> (INT_MAX)) {
      typedef MV_Nrm2Squared_Right_FunctorVector<RV, XMV, int> functor_type;
      Kokkos::RangePolicy<execution_space, int> policy (0, numRows);
      functor_type op (r, X);
      Kokkos::parallel_reduce (policy, op);
    }
    else {
      typedef MV_Nrm2Squared_Right_FunctorVector<RV, XMV, size_type> functor_type;
      Kokkos::RangePolicy<execution_space, size_type> policy (0, numRows);
      functor_type op (r, X);
      Kokkos::parallel_reduce (policy, op);
    }
  }

  /// \brief Compute the square of the 2-norm of X(:,X_col), and store
  ///   result in r(r_col).
  static void nrm2_squared (const RV& r, const size_t r_col, const XMV& X, const size_t X_col)
  {
    using Kokkos::ALL;
    using Kokkos::subview;
    typedef Kokkos::View<typename RV::value_type,
      typename RV::array_layout,
      typename RV::device_type,
      typename RV::memory_traits,
      typename RV::specialize> RV0D;
    typedef Kokkos::View<typename XMV::const_value_type*,
      typename XMV::array_layout,
      typename XMV::device_type,
      typename XMV::memory_traits,
      typename XMV::specialize> XMV1D;

    const size_type numRows = X.dimension_0 ();
    const size_type numCols = X.dimension_1 ();

    // int is generally faster than size_t, but check for overflow first.
    if (numRows < static_cast<size_type> (INT_MAX) &&
        numRows * numCols < static_cast<size_type> (INT_MAX)) {
      typedef V_Nrm2Squared_Functor<RV0D, XMV1D, int> functor_type;
      Kokkos::RangePolicy<execution_space, int> policy (0, numRows);
      functor_type op (subview (r, r_col), subview (X, ALL (), X_col));
      Kokkos::parallel_reduce (policy, op);
    }
    else {
      typedef V_Nrm2Squared_Functor<RV0D, XMV1D, size_type> functor_type;
      Kokkos::RangePolicy<execution_space, size_type> policy (0, numRows);
      functor_type op (subview (r, r_col), subview (X, ALL (), X_col));
      Kokkos::parallel_reduce (policy, op);
    }
  }
};

// Full specializations for cases of interest for Tpetra::MultiVector.
//
// Currently, we include specializations for Scalar = double,
// LayoutLeft (which is what Tpetra::MultiVector uses at the moment),
// and all execution spaces.  This may change in the future.  The
// output 1-D View _always_ uses the execution space's default array
// layout, which is what Tpetra::MultiVector wants for the output
// argument of norm2().

#ifdef KOKKOS_HAVE_SERIAL
template<>
struct Nrm2_MV<double*,
               Kokkos::Serial::array_layout,
               Kokkos::Device<Kokkos::Serial, Kokkos::HostSpace>,
               Kokkos::MemoryTraits<Kokkos::Unmanaged>,
               Kokkos::Impl::ViewDefault,
               const double**,
               Kokkos::LayoutLeft,
               Kokkos::Device<Kokkos::Serial, Kokkos::HostSpace>,
               Kokkos::MemoryTraits<Kokkos::Unmanaged>,
               Kokkos::Impl::ViewDefault> {
  typedef double* RT;
  typedef Kokkos::Serial::array_layout RL;
  typedef Kokkos::Device<Kokkos::Serial, Kokkos::HostSpace> RD;
  typedef Kokkos::MemoryTraits<Kokkos::Unmanaged> RM;
  typedef Kokkos::Impl::ViewDefault RS;
  typedef Kokkos::View<RT,RL,RD,RM,RS> RV;

  typedef const double** XT;
  typedef Kokkos::LayoutLeft XL;
  typedef Kokkos::Device<Kokkos::Serial, Kokkos::HostSpace> XD;
  typedef Kokkos::MemoryTraits<Kokkos::Unmanaged> XM;
  typedef Kokkos::Impl::ViewDefault XS;
  typedef Kokkos::View<XT,XL,XD,XM,XS> XMV;

  typedef XMV::execution_space execution_space;
  typedef XMV::size_type size_type;

  static void nrm2_squared (const RV& r, const XMV& X);
  static void nrm2_squared (const RV& r, const size_t r_col, const XMV& X, const size_t X_col);
};
#endif // KOKKOS_HAVE_SERIAL

#ifdef KOKKOS_HAVE_OPENMP
template<>
struct Nrm2_MV<double*,
               Kokkos::OpenMP::array_layout,
               Kokkos::Device<Kokkos::OpenMP, Kokkos::HostSpace>,
               Kokkos::MemoryTraits<Kokkos::Unmanaged>,
               Kokkos::Impl::ViewDefault,
               const double**,
               Kokkos::LayoutLeft,
               Kokkos::Device<Kokkos::OpenMP, Kokkos::HostSpace>,
               Kokkos::MemoryTraits<Kokkos::Unmanaged>,
               Kokkos::Impl::ViewDefault> {
  typedef double* RT;
  typedef Kokkos::OpenMP::array_layout RL;
  typedef Kokkos::Device<Kokkos::OpenMP, Kokkos::HostSpace> RD;
  typedef Kokkos::MemoryTraits<Kokkos::Unmanaged> RM;
  typedef Kokkos::Impl::ViewDefault RS;
  typedef Kokkos::View<RT,RL,RD,RM,RS> RV;

  typedef const double** XT;
  typedef Kokkos::LayoutLeft XL;
  typedef Kokkos::Device<Kokkos::OpenMP, Kokkos::HostSpace> XD;
  typedef Kokkos::MemoryTraits<Kokkos::Unmanaged> XM;
  typedef Kokkos::Impl::ViewDefault XS;
  typedef Kokkos::View<XT,XL,XD,XM,XS> XMV;

  typedef XMV::execution_space execution_space;
  typedef XMV::size_type size_type;

  static void nrm2_squared (const RV& r, const XMV& X);
  static void nrm2_squared (const RV& r, const size_t r_col, const XMV& X, const size_t X_col);
};
#endif // KOKKOS_HAVE_OPENMP

#ifdef KOKKOS_HAVE_PTHREAD
template<>
struct Nrm2_MV<double*,
               Kokkos::Threads::array_layout,
               Kokkos::Device<Kokkos::Threads, Kokkos::HostSpace>,
               Kokkos::MemoryTraits<Kokkos::Unmanaged>,
               Kokkos::Impl::ViewDefault,
               const double**,
               Kokkos::LayoutLeft,
               Kokkos::Device<Kokkos::Threads, Kokkos::HostSpace>,
               Kokkos::MemoryTraits<Kokkos::Unmanaged>,
               Kokkos::Impl::ViewDefault> {
  typedef double* RT;
  typedef Kokkos::Threads::array_layout RL;
  typedef Kokkos::Device<Kokkos::Threads, Kokkos::HostSpace> RD;
  typedef Kokkos::MemoryTraits<Kokkos::Unmanaged> RM;
  typedef Kokkos::Impl::ViewDefault RS;
  typedef Kokkos::View<RT,RL,RD,RM,RS> RV;

  typedef const double** XT;
  typedef Kokkos::LayoutLeft XL;
  typedef Kokkos::Device<Kokkos::Threads, Kokkos::HostSpace> XD;
  typedef Kokkos::MemoryTraits<Kokkos::Unmanaged> XM;
  typedef Kokkos::Impl::ViewDefault XS;
  typedef Kokkos::View<XT,XL,XD,XM,XS> XMV;

  typedef XMV::execution_space execution_space;
  typedef XMV::size_type size_type;

  static void nrm2_squared (const RV& r, const XMV& X);
  static void nrm2_squared (const RV& r, const size_t r_col, const XMV& X, const size_t X_col);
};
#endif // KOKKOS_HAVE_PTHREAD

#ifdef KOKKOS_HAVE_CUDA
template<>
struct Nrm2_MV<double*,
               Kokkos::Cuda::array_layout,
               Kokkos::Device<Kokkos::Cuda, Kokkos::CudaSpace>,
               Kokkos::MemoryTraits<Kokkos::Unmanaged>,
               Kokkos::Impl::ViewDefault,
               const double**,
               Kokkos::LayoutLeft,
               Kokkos::Device<Kokkos::Cuda, Kokkos::CudaSpace>,
               Kokkos::MemoryTraits<Kokkos::Unmanaged>,
               Kokkos::Impl::ViewDefault> {
  typedef double* RT;
  typedef Kokkos::Cuda::array_layout RL;
  typedef Kokkos::Device<Kokkos::Cuda, Kokkos::CudaSpace> RD;
  typedef Kokkos::MemoryTraits<Kokkos::Unmanaged> RM;
  typedef Kokkos::Impl::ViewDefault RS;
  typedef Kokkos::View<RT,RL,RD,RM,RS> RV;

  typedef const double** XT;
  typedef Kokkos::LayoutLeft XL;
  typedef Kokkos::Device<Kokkos::Cuda, Kokkos::CudaSpace> XD;
  typedef Kokkos::MemoryTraits<Kokkos::Unmanaged> XM;
  typedef Kokkos::Impl::ViewDefault XS;
  typedef Kokkos::View<XT,XL,XD,XM,XS> XMV;

  typedef XMV::execution_space execution_space;
  typedef XMV::size_type size_type;

  static void nrm2_squared (const RV& r, const XMV& X);
  static void nrm2_squared (const RV& r, const size_t r_col, const XMV& X, const size_t X_col);
};
#endif // KOKKOS_HAVE_CUDA

#ifdef KOKKOS_HAVE_CUDA
template<>
struct Nrm2_MV<double*,
               Kokkos::Cuda::array_layout,
               Kokkos::Device<Kokkos::Cuda, Kokkos::CudaUVMSpace>,
               Kokkos::MemoryTraits<Kokkos::Unmanaged>,
               Kokkos::Impl::ViewDefault,
               const double**,
               Kokkos::LayoutLeft,
               Kokkos::Device<Kokkos::Cuda, Kokkos::CudaUVMSpace>,
               Kokkos::MemoryTraits<Kokkos::Unmanaged>,
               Kokkos::Impl::ViewDefault> {
  typedef double* RT;
  typedef Kokkos::Cuda::array_layout RL;
  typedef Kokkos::Device<Kokkos::Cuda, Kokkos::CudaUVMSpace> RD;
  typedef Kokkos::MemoryTraits<Kokkos::Unmanaged> RM;
  typedef Kokkos::Impl::ViewDefault RS;
  typedef Kokkos::View<RT,RL,RD,RM,RS> RV;

  typedef const double** XT;
  typedef Kokkos::LayoutLeft XL;
  typedef Kokkos::Device<Kokkos::Cuda, Kokkos::CudaUVMSpace> XD;
  typedef Kokkos::MemoryTraits<Kokkos::Unmanaged> XM;
  typedef Kokkos::Impl::ViewDefault XS;
  typedef Kokkos::View<XT,XL,XD,XM,XS> XMV;

  typedef XMV::execution_space execution_space;
  typedef XMV::size_type size_type;

  static void nrm2_squared (const RV& r, const XMV& X);
  static void nrm2_squared (const RV& r, const size_t r_col, const XMV& X, const size_t X_col);
};
#endif // KOKKOS_HAVE_CUDA

//
// nrm1
//

/// \brief 1-norm functor for single vectors.
///
/// \tparam RV 0-D output View
/// \tparam XV 1-D input View
/// \tparam SizeType Index type.  Use int (32 bits) if possible.
template<class RV, class XV, class SizeType = typename XV::size_type>
struct V_Nrm1_Functor
{
  typedef typename XV::execution_space              execution_space;
  typedef SizeType                                        size_type;
  typedef typename XV::non_const_value_type             xvalue_type;
  typedef Kokkos::Details::InnerProductSpaceTraits<xvalue_type> IPT;
  typedef Kokkos::Details::ArithTraits<typename IPT::mag_type>   AT;
  typedef typename IPT::mag_type                         value_type;

  RV m_r;
  XV m_x;

  V_Nrm1_Functor (const RV& r, const XV& x) :
    m_r (r), m_x (x)
  {}

  KOKKOS_INLINE_FUNCTION
  void operator() (const size_type& i, value_type& sum) const
  {
    sum += IPT::norm (m_x(i)); // absolute value
  }

  KOKKOS_INLINE_FUNCTION void init (value_type& update) const
  {
    update = AT::zero ();
  }

  KOKKOS_INLINE_FUNCTION void
  join (volatile value_type& update,
        const volatile value_type& source) const
  {
    update += source;
  }

  // On device, write the reduction result to the output View.
  KOKKOS_INLINE_FUNCTION void final (const value_type& dst) const
  {
    m_r() = dst;
  }
};

/// \brief 1-norm functor for multivectors.
///
/// \tparam RV 1-D output View
/// \tparam XMV 2-D input View
/// \tparam SizeType Index type.  Use int (32 bits) if possible.
template<class RV, class XMV, class SizeType = typename XMV::size_type>
struct MV_Nrm1_Functor {
  typedef typename XMV::execution_space             execution_space;
  typedef SizeType                                        size_type;
  typedef typename XMV::non_const_value_type            xvalue_type;
  typedef Kokkos::Details::InnerProductSpaceTraits<xvalue_type> IPT;
  typedef Kokkos::Details::ArithTraits<typename IPT::mag_type>   AT;
  typedef typename IPT::mag_type                       value_type[];

  const size_type value_count;
  RV norms_;
  XMV X_;

  MV_Nrm1_Functor (const RV& norms, const XMV& X) :
    value_count (X.dimension_1 ()), norms_ (norms), X_ (X)
  {}

  KOKKOS_INLINE_FUNCTION void
  operator() (const size_type& i, value_type sum) const
  {
#ifdef KOKKOS_HAVE_PRAGMA_IVDEP
#pragma ivdep
#endif
#ifdef KOKKOS_HAVE_PRAGMA_VECTOR
#pragma vector always
#endif
    for (size_type k = 0; k < value_count; ++k) {
      sum[k] += IPT::norm (X_(i,k)); // absolute value
    }
  }

  KOKKOS_INLINE_FUNCTION void
  init (value_type update) const
  {
#ifdef KOKKOS_HAVE_PRAGMA_IVDEP
#pragma ivdep
#endif
#ifdef KOKKOS_HAVE_PRAGMA_VECTOR
#pragma vector always
#endif
    for (size_type k = 0; k < value_count; ++k) {
      update[k] = AT::zero ();
    }
  }

  KOKKOS_INLINE_FUNCTION void
  join (volatile value_type update,
        const volatile value_type source) const
  {
#ifdef KOKKOS_HAVE_PRAGMA_IVDEP
#pragma ivdep
#endif
#ifdef KOKKOS_HAVE_PRAGMA_VECTOR
#pragma vector always
#endif
    for (size_type k = 0; k < value_count; ++k) {
      update[k] += source[k];
    }
  }

  // On device, write the reduction result to the output View.
  KOKKOS_INLINE_FUNCTION void
  final (const value_type dst) const
  {
#ifdef KOKKOS_HAVE_PRAGMA_IVDEP
#pragma ivdep
#endif
#ifdef KOKKOS_HAVE_PRAGMA_VECTOR
#pragma vector always
#endif
    for (size_type k = 0; k < value_count; ++k) {
      norms_(k) = dst[k];
    }
  }
};

//! Implementation of KokkosBlas::nrm1 for multivectors.
template<class RT, class RL, class RD, class RM, class RS,
         class XT, class XL, class XD, class XM, class XS>
struct Nrm1_MV {
  typedef Kokkos::View<RT,RL,RD,RM,RS> RV;
  typedef Kokkos::View<XT,XL,XD,XM,XS> XMV;
  typedef typename XMV::execution_space execution_space;
  typedef typename XMV::size_type size_type;

  /// \brief Compute the 1-norm(s) of the column(s) of the multivector
  ///   (2-D View) X, and store result(s) in r.
  static void nrm1 (const RV& r, const XMV& X)
  {
    const size_type numRows = X.dimension_0 ();
    const size_type numCols = X.dimension_1 ();

    // int is generally faster than size_t, but check for overflow first.
    if (numRows < static_cast<size_type> (INT_MAX) &&
        numRows * numCols < static_cast<size_type> (INT_MAX)) {
      typedef MV_Nrm1_Functor<RV, XMV, int> functor_type;
      Kokkos::RangePolicy<execution_space, int> policy (0, numRows);
      functor_type op (r, X);
      Kokkos::parallel_reduce (policy, op);
    }
    else {
      typedef MV_Nrm1_Functor<RV, XMV, size_type> functor_type;
      Kokkos::RangePolicy<execution_space, size_type> policy (0, numRows);
      functor_type op (r, X);
      Kokkos::parallel_reduce (policy, op);
    }
  }

  //! Compute the 1-norm of X(:,X_col), and store result in r(r_col).
  static void nrm1 (const RV& r, const size_t r_col, const XMV& X, const size_t X_col)
  {
    using Kokkos::ALL;
    using Kokkos::subview;
    typedef Kokkos::View<typename RV::value_type,
      typename RV::array_layout,
      typename RV::device_type,
      typename RV::memory_traits,
      typename RV::specialize> RV0D;
    typedef Kokkos::View<typename XMV::const_value_type*,
      typename XMV::array_layout,
      typename XMV::device_type,
      typename XMV::memory_traits,
      typename XMV::specialize> XMV1D;

    const size_type numRows = X.dimension_0 ();
    const size_type numCols = X.dimension_1 ();

    // int is generally faster than size_t, but check for overflow first.
    if (numRows < static_cast<size_type> (INT_MAX) &&
        numRows * numCols < static_cast<size_type> (INT_MAX)) {
      typedef V_Nrm1_Functor<RV0D, XMV1D, int> functor_type;
      Kokkos::RangePolicy<execution_space, int> policy (0, numRows);
      functor_type op (subview (r, r_col), subview (X, ALL (), X_col));
      Kokkos::parallel_reduce (policy, op);
    }
    else {
      typedef V_Nrm1_Functor<RV0D, XMV1D, size_type> functor_type;
      Kokkos::RangePolicy<execution_space, size_type> policy (0, numRows);
      functor_type op (subview (r, r_col), subview (X, ALL (), X_col));
      Kokkos::parallel_reduce (policy, op);
    }
  }
};

// Full specializations for cases of interest for Tpetra::MultiVector.
//
// Currently, we include specializations for Scalar = double,
// LayoutLeft (which is what Tpetra::MultiVector uses at the moment),
// and all execution spaces.  This may change in the future.  The
// output 1-D View _always_ uses the execution space's default array
// layout, which is what Tpetra::MultiVector wants for the output
// argument of norm1().

#ifdef KOKKOS_HAVE_SERIAL
template<>
struct Nrm1_MV<double*,
               Kokkos::Serial::array_layout,
               Kokkos::Device<Kokkos::Serial, Kokkos::HostSpace>,
               Kokkos::MemoryTraits<Kokkos::Unmanaged>,
               Kokkos::Impl::ViewDefault,
               const double**,
               Kokkos::LayoutLeft,
               Kokkos::Device<Kokkos::Serial, Kokkos::HostSpace>,
               Kokkos::MemoryTraits<Kokkos::Unmanaged>,
               Kokkos::Impl::ViewDefault> {
  typedef double* RT;
  typedef Kokkos::Serial::array_layout RL;
  typedef Kokkos::Device<Kokkos::Serial, Kokkos::HostSpace> RD;
  typedef Kokkos::MemoryTraits<Kokkos::Unmanaged> RM;
  typedef Kokkos::Impl::ViewDefault RS;
  typedef Kokkos::View<RT,RL,RD,RM,RS> RV;

  typedef const double** XT;
  typedef Kokkos::LayoutLeft XL;
  typedef Kokkos::Device<Kokkos::Serial, Kokkos::HostSpace> XD;
  typedef Kokkos::MemoryTraits<Kokkos::Unmanaged> XM;
  typedef Kokkos::Impl::ViewDefault XS;
  typedef Kokkos::View<XT,XL,XD,XM,XS> XMV;

  typedef XMV::execution_space execution_space;
  typedef XMV::size_type size_type;

  static void nrm1 (const RV& r, const XMV& X);
  static void nrm1 (const RV& r, const size_t r_col, const XMV& X, const size_t X_col);
};
#endif // KOKKOS_HAVE_SERIAL

#ifdef KOKKOS_HAVE_OPENMP
template<>
struct Nrm1_MV<double*,
               Kokkos::OpenMP::array_layout,
               Kokkos::Device<Kokkos::OpenMP, Kokkos::HostSpace>,
               Kokkos::MemoryTraits<Kokkos::Unmanaged>,
               Kokkos::Impl::ViewDefault,
               const double**,
               Kokkos::LayoutLeft,
               Kokkos::Device<Kokkos::OpenMP, Kokkos::HostSpace>,
               Kokkos::MemoryTraits<Kokkos::Unmanaged>,
               Kokkos::Impl::ViewDefault> {
  typedef double* RT;
  typedef Kokkos::OpenMP::array_layout RL;
  typedef Kokkos::Device<Kokkos::OpenMP, Kokkos::HostSpace> RD;
  typedef Kokkos::MemoryTraits<Kokkos::Unmanaged> RM;
  typedef Kokkos::Impl::ViewDefault RS;
  typedef Kokkos::View<RT,RL,RD,RM,RS> RV;

  typedef const double** XT;
  typedef Kokkos::LayoutLeft XL;
  typedef Kokkos::Device<Kokkos::OpenMP, Kokkos::HostSpace> XD;
  typedef Kokkos::MemoryTraits<Kokkos::Unmanaged> XM;
  typedef Kokkos::Impl::ViewDefault XS;
  typedef Kokkos::View<XT,XL,XD,XM,XS> XMV;

  typedef XMV::execution_space execution_space;
  typedef XMV::size_type size_type;

  static void nrm1 (const RV& r, const XMV& X);
  static void nrm1 (const RV& r, const size_t r_col, const XMV& X, const size_t X_col);
};
#endif // KOKKOS_HAVE_OPENMP

#ifdef KOKKOS_HAVE_PTHREAD
template<>
struct Nrm1_MV<double*,
               Kokkos::Threads::array_layout,
               Kokkos::Device<Kokkos::Threads, Kokkos::HostSpace>,
               Kokkos::MemoryTraits<Kokkos::Unmanaged>,
               Kokkos::Impl::ViewDefault,
               const double**,
               Kokkos::LayoutLeft,
               Kokkos::Device<Kokkos::Threads, Kokkos::HostSpace>,
               Kokkos::MemoryTraits<Kokkos::Unmanaged>,
               Kokkos::Impl::ViewDefault> {
  typedef double* RT;
  typedef Kokkos::Threads::array_layout RL;
  typedef Kokkos::Device<Kokkos::Threads, Kokkos::HostSpace> RD;
  typedef Kokkos::MemoryTraits<Kokkos::Unmanaged> RM;
  typedef Kokkos::Impl::ViewDefault RS;
  typedef Kokkos::View<RT,RL,RD,RM,RS> RV;

  typedef const double** XT;
  typedef Kokkos::LayoutLeft XL;
  typedef Kokkos::Device<Kokkos::Threads, Kokkos::HostSpace> XD;
  typedef Kokkos::MemoryTraits<Kokkos::Unmanaged> XM;
  typedef Kokkos::Impl::ViewDefault XS;
  typedef Kokkos::View<XT,XL,XD,XM,XS> XMV;

  typedef XMV::execution_space execution_space;
  typedef XMV::size_type size_type;

  static void nrm1 (const RV& r, const XMV& X);
  static void nrm1 (const RV& r, const size_t r_col, const XMV& X, const size_t X_col);
};
#endif // KOKKOS_HAVE_PTHREAD

#ifdef KOKKOS_HAVE_CUDA
template<>
struct Nrm1_MV<double*,
               Kokkos::Cuda::array_layout,
               Kokkos::Device<Kokkos::Cuda, Kokkos::CudaSpace>,
               Kokkos::MemoryTraits<Kokkos::Unmanaged>,
               Kokkos::Impl::ViewDefault,
               const double**,
               Kokkos::LayoutLeft,
               Kokkos::Device<Kokkos::Cuda, Kokkos::CudaSpace>,
               Kokkos::MemoryTraits<Kokkos::Unmanaged>,
               Kokkos::Impl::ViewDefault> {
  typedef double* RT;
  typedef Kokkos::Cuda::array_layout RL;
  typedef Kokkos::Device<Kokkos::Cuda, Kokkos::CudaSpace> RD;
  typedef Kokkos::MemoryTraits<Kokkos::Unmanaged> RM;
  typedef Kokkos::Impl::ViewDefault RS;
  typedef Kokkos::View<RT,RL,RD,RM,RS> RV;

  typedef const double** XT;
  typedef Kokkos::LayoutLeft XL;
  typedef Kokkos::Device<Kokkos::Cuda, Kokkos::CudaSpace> XD;
  typedef Kokkos::MemoryTraits<Kokkos::Unmanaged> XM;
  typedef Kokkos::Impl::ViewDefault XS;
  typedef Kokkos::View<XT,XL,XD,XM,XS> XMV;

  typedef XMV::execution_space execution_space;
  typedef XMV::size_type size_type;

  static void nrm1 (const RV& r, const XMV& X);
  static void nrm1 (const RV& r, const size_t r_col, const XMV& X, const size_t X_col);
};
#endif // KOKKOS_HAVE_CUDA

#ifdef KOKKOS_HAVE_CUDA
template<>
struct Nrm1_MV<double*,
               Kokkos::Cuda::array_layout,
               Kokkos::Device<Kokkos::Cuda, Kokkos::CudaUVMSpace>,
               Kokkos::MemoryTraits<Kokkos::Unmanaged>,
               Kokkos::Impl::ViewDefault,
               const double**,
               Kokkos::LayoutLeft,
               Kokkos::Device<Kokkos::Cuda, Kokkos::CudaUVMSpace>,
               Kokkos::MemoryTraits<Kokkos::Unmanaged>,
               Kokkos::Impl::ViewDefault> {
  typedef double* RT;
  typedef Kokkos::Cuda::array_layout RL;
  typedef Kokkos::Device<Kokkos::Cuda, Kokkos::CudaUVMSpace> RD;
  typedef Kokkos::MemoryTraits<Kokkos::Unmanaged> RM;
  typedef Kokkos::Impl::ViewDefault RS;
  typedef Kokkos::View<RT,RL,RD,RM,RS> RV;

  typedef const double** XT;
  typedef Kokkos::LayoutLeft XL;
  typedef Kokkos::Device<Kokkos::Cuda, Kokkos::CudaUVMSpace> XD;
  typedef Kokkos::MemoryTraits<Kokkos::Unmanaged> XM;
  typedef Kokkos::Impl::ViewDefault XS;
  typedef Kokkos::View<XT,XL,XD,XM,XS> XMV;

  typedef XMV::execution_space execution_space;
  typedef XMV::size_type size_type;

  static void nrm1 (const RV& r, const XMV& X);
  static void nrm1 (const RV& r, const size_t r_col, const XMV& X, const size_t X_col);
};
#endif // KOKKOS_HAVE_CUDA

//
// nrmInf
//

/// \brief Inf-norm functor for single vectors.
///
/// \tparam RV 0-D output View
/// \tparam XV 1-D input View
/// \tparam SizeType Index type.  Use int (32 bits) if possible.
template<class RV, class XV, class SizeType = typename XV::size_type>
struct V_NrmInf_Functor
{
  typedef typename XV::execution_space              execution_space;
  typedef SizeType                                        size_type;
  typedef typename XV::non_const_value_type             xvalue_type;
  typedef Kokkos::Details::InnerProductSpaceTraits<xvalue_type> IPT;
  typedef Kokkos::Details::ArithTraits<typename IPT::mag_type>   AT;
  typedef typename IPT::mag_type                         value_type;

  RV m_r;
  XV m_x;

  V_NrmInf_Functor (const RV& r, const XV& x) :
    m_r (r), m_x (x)
  {}

  KOKKOS_INLINE_FUNCTION
  void operator() (const size_type& i, value_type& update) const
  {
    const typename IPT::mag_type tmp = IPT::norm (m_x(i));
    if (update < tmp) {
      update = tmp;
    }
  }

  KOKKOS_INLINE_FUNCTION void init (value_type& update) const
  {
    update = AT::zero ();
  }

  KOKKOS_INLINE_FUNCTION void
  join (volatile value_type& update,
        const volatile value_type& source) const
  {
    if (update < source) {
      update = source;
    }
  }

  // On device, write the reduction result to the output View.
  KOKKOS_INLINE_FUNCTION void final (const value_type& dst) const
  {
    m_r() = dst;
  }
};

/// \brief Inf-norm functor for multivectors.
///
/// \tparam RV 1-D output View
/// \tparam XMV 2-D input View
/// \tparam SizeType Index type.  Use int (32 bits) if possible.
template<class RV, class XMV, class SizeType = typename XMV::size_type>
struct MV_NrmInf_Functor {
  typedef typename XMV::execution_space             execution_space;
  typedef SizeType                                        size_type;
  typedef typename XMV::non_const_value_type            xvalue_type;
  typedef Kokkos::Details::InnerProductSpaceTraits<xvalue_type> IPT;
  typedef Kokkos::Details::ArithTraits<typename IPT::mag_type>   AT;
  typedef typename IPT::mag_type                       value_type[];

  const size_type value_count;
  RV norms_;
  XMV X_;

  MV_NrmInf_Functor (const RV& norms, const XMV& X) :
    value_count (X.dimension_1 ()), norms_ (norms), X_ (X)
  {}

  KOKKOS_INLINE_FUNCTION void
  operator() (const size_type& i, value_type update) const
  {
#ifdef KOKKOS_HAVE_PRAGMA_IVDEP
#pragma ivdep
#endif
#ifdef KOKKOS_HAVE_PRAGMA_VECTOR
#pragma vector always
#endif
    for (size_type j = 0; j < value_count; ++j) {
      const typename IPT::mag_type tmp = IPT::norm (X_(i,j));
      if (update[j] < tmp) {
        update[j] = tmp;
      }
    }
  }

  KOKKOS_INLINE_FUNCTION void
  init (value_type update) const
  {
#ifdef KOKKOS_HAVE_PRAGMA_IVDEP
#pragma ivdep
#endif
#ifdef KOKKOS_HAVE_PRAGMA_VECTOR
#pragma vector always
#endif
    for (size_type k = 0; k < value_count; ++k) {
      update[k] = AT::zero ();
    }
  }

  KOKKOS_INLINE_FUNCTION void
  join (volatile value_type update,
        const volatile value_type source) const
  {
#ifdef KOKKOS_HAVE_PRAGMA_IVDEP
#pragma ivdep
#endif
#ifdef KOKKOS_HAVE_PRAGMA_VECTOR
#pragma vector always
#endif
    for (size_type k = 0; k < value_count; ++k) {
      if (update[k] < source[k]) {
        update[k] = source[k];
      }
    }
  }

  // On device, write the reduction result to the output View.
  KOKKOS_INLINE_FUNCTION void
  final (const value_type dst) const
  {
#ifdef KOKKOS_HAVE_PRAGMA_IVDEP
#pragma ivdep
#endif
#ifdef KOKKOS_HAVE_PRAGMA_VECTOR
#pragma vector always
#endif
    for (size_type k = 0; k < value_count; ++k) {
      norms_(k) = dst[k];
    }
  }
};

//! Implementation of KokkosBlas::nrmInf for multivectors.
template<class RT, class RL, class RD, class RM, class RS,
         class XT, class XL, class XD, class XM, class XS>
struct NrmInf_MV {
  typedef Kokkos::View<RT,RL,RD,RM,RS> RV;
  typedef Kokkos::View<XT,XL,XD,XM,XS> XMV;
  typedef typename XMV::execution_space execution_space;
  typedef typename XMV::size_type size_type;

  /// \brief Compute the inf-norm(s) of the column(s) of the
  ///   multivector (2-D View) X, and store result(s) in r.
  static void nrmInf (const RV& r, const XMV& X)
  {
    const size_type numRows = X.dimension_0 ();
    const size_type numCols = X.dimension_1 ();

    // int is generally faster than size_t, but check for overflow first.
    if (numRows < static_cast<size_type> (INT_MAX) &&
        numRows * numCols < static_cast<size_type> (INT_MAX)) {
      typedef MV_NrmInf_Functor<RV, XMV, int> functor_type;
      Kokkos::RangePolicy<execution_space, int> policy (0, numRows);
      functor_type op (r, X);
      Kokkos::parallel_reduce (policy, op);
    }
    else {
      typedef MV_NrmInf_Functor<RV, XMV, size_type> functor_type;
      Kokkos::RangePolicy<execution_space, size_type> policy (0, numRows);
      functor_type op (r, X);
      Kokkos::parallel_reduce (policy, op);
    }
  }

  //! Compute the inf-norm of X(:,X_col), and store result in r(r_col).
  static void nrmInf (const RV& r, const size_t r_col, const XMV& X, const size_t X_col)
  {
    using Kokkos::ALL;
    using Kokkos::subview;
    typedef Kokkos::View<typename RV::value_type,
      typename RV::array_layout,
      typename RV::device_type,
      typename RV::memory_traits,
      typename RV::specialize> RV0D;
    typedef Kokkos::View<typename XMV::const_value_type*,
      typename XMV::array_layout,
      typename XMV::device_type,
      typename XMV::memory_traits,
      typename XMV::specialize> XMV1D;

    const size_type numRows = X.dimension_0 ();
    const size_type numCols = X.dimension_1 ();

    // int is generally faster than size_t, but check for overflow first.
    if (numRows < static_cast<size_type> (INT_MAX) &&
        numRows * numCols < static_cast<size_type> (INT_MAX)) {
      typedef V_NrmInf_Functor<RV0D, XMV1D, int> functor_type;
      Kokkos::RangePolicy<execution_space, int> policy (0, numRows);
      functor_type op (subview (r, r_col), subview (X, ALL (), X_col));
      Kokkos::parallel_reduce (policy, op);
    }
    else {
      typedef V_NrmInf_Functor<RV0D, XMV1D, size_type> functor_type;
      Kokkos::RangePolicy<execution_space, size_type> policy (0, numRows);
      functor_type op (subview (r, r_col), subview (X, ALL (), X_col));
      Kokkos::parallel_reduce (policy, op);
    }
  }
};

// Full specializations for cases of interest for Tpetra::MultiVector.
//
// Currently, we include specializations for Scalar = double,
// LayoutLeft (which is what Tpetra::MultiVector uses at the moment),
// and all execution spaces.  This may change in the future.  The
// output 1-D View _always_ uses the execution space's default array
// layout, which is what Tpetra::MultiVector wants for the output
// argument of normInf().

#ifdef KOKKOS_HAVE_SERIAL
template<>
struct NrmInf_MV<double*,
                 Kokkos::Serial::array_layout,
                 Kokkos::Device<Kokkos::Serial, Kokkos::HostSpace>,
                 Kokkos::MemoryTraits<Kokkos::Unmanaged>,
                 Kokkos::Impl::ViewDefault,
                 const double**,
                 Kokkos::LayoutLeft,
                 Kokkos::Device<Kokkos::Serial, Kokkos::HostSpace>,
                 Kokkos::MemoryTraits<Kokkos::Unmanaged>,
                 Kokkos::Impl::ViewDefault> {
  typedef double* RT;
  typedef Kokkos::Serial::array_layout RL;
  typedef Kokkos::Device<Kokkos::Serial, Kokkos::HostSpace> RD;
  typedef Kokkos::MemoryTraits<Kokkos::Unmanaged> RM;
  typedef Kokkos::Impl::ViewDefault RS;
  typedef Kokkos::View<RT,RL,RD,RM,RS> RV;

  typedef const double** XT;
  typedef Kokkos::LayoutLeft XL;
  typedef Kokkos::Device<Kokkos::Serial, Kokkos::HostSpace> XD;
  typedef Kokkos::MemoryTraits<Kokkos::Unmanaged> XM;
  typedef Kokkos::Impl::ViewDefault XS;
  typedef Kokkos::View<XT,XL,XD,XM,XS> XMV;

  typedef XMV::execution_space execution_space;
  typedef XMV::size_type size_type;

  static void nrmInf (const RV& r, const XMV& X);
  static void nrmInf (const RV& r, const size_t r_col, const XMV& X, const size_t X_col);
};
#endif // KOKKOS_HAVE_SERIAL

#ifdef KOKKOS_HAVE_OPENMP
template<>
struct NrmInf_MV<double*,
                 Kokkos::OpenMP::array_layout,
                 Kokkos::Device<Kokkos::OpenMP, Kokkos::HostSpace>,
                 Kokkos::MemoryTraits<Kokkos::Unmanaged>,
                 Kokkos::Impl::ViewDefault,
                 const double**,
                 Kokkos::LayoutLeft,
                 Kokkos::Device<Kokkos::OpenMP, Kokkos::HostSpace>,
                 Kokkos::MemoryTraits<Kokkos::Unmanaged>,
                 Kokkos::Impl::ViewDefault> {
  typedef double* RT;
  typedef Kokkos::OpenMP::array_layout RL;
  typedef Kokkos::Device<Kokkos::OpenMP, Kokkos::HostSpace> RD;
  typedef Kokkos::MemoryTraits<Kokkos::Unmanaged> RM;
  typedef Kokkos::Impl::ViewDefault RS;
  typedef Kokkos::View<RT,RL,RD,RM,RS> RV;

  typedef const double** XT;
  typedef Kokkos::LayoutLeft XL;
  typedef Kokkos::Device<Kokkos::OpenMP, Kokkos::HostSpace> XD;
  typedef Kokkos::MemoryTraits<Kokkos::Unmanaged> XM;
  typedef Kokkos::Impl::ViewDefault XS;
  typedef Kokkos::View<XT,XL,XD,XM,XS> XMV;

  typedef XMV::execution_space execution_space;
  typedef XMV::size_type size_type;

  static void nrmInf (const RV& r, const XMV& X);
  static void nrmInf (const RV& r, const size_t r_col, const XMV& X, const size_t X_col);
};
#endif // KOKKOS_HAVE_OPENMP

#ifdef KOKKOS_HAVE_PTHREAD
template<>
struct NrmInf_MV<double*,
                 Kokkos::Threads::array_layout,
                 Kokkos::Device<Kokkos::Threads, Kokkos::HostSpace>,
                 Kokkos::MemoryTraits<Kokkos::Unmanaged>,
                 Kokkos::Impl::ViewDefault,
                 const double**,
                 Kokkos::LayoutLeft,
                 Kokkos::Device<Kokkos::Threads, Kokkos::HostSpace>,
                 Kokkos::MemoryTraits<Kokkos::Unmanaged>,
                 Kokkos::Impl::ViewDefault> {
  typedef double* RT;
  typedef Kokkos::Threads::array_layout RL;
  typedef Kokkos::Device<Kokkos::Threads, Kokkos::HostSpace> RD;
  typedef Kokkos::MemoryTraits<Kokkos::Unmanaged> RM;
  typedef Kokkos::Impl::ViewDefault RS;
  typedef Kokkos::View<RT,RL,RD,RM,RS> RV;

  typedef const double** XT;
  typedef Kokkos::LayoutLeft XL;
  typedef Kokkos::Device<Kokkos::Threads, Kokkos::HostSpace> XD;
  typedef Kokkos::MemoryTraits<Kokkos::Unmanaged> XM;
  typedef Kokkos::Impl::ViewDefault XS;
  typedef Kokkos::View<XT,XL,XD,XM,XS> XMV;

  typedef XMV::execution_space execution_space;
  typedef XMV::size_type size_type;

  static void nrmInf (const RV& r, const XMV& X);
  static void nrmInf (const RV& r, const size_t r_col, const XMV& X, const size_t X_col);
};
#endif // KOKKOS_HAVE_PTHREAD

#ifdef KOKKOS_HAVE_CUDA
template<>
struct NrmInf_MV<double*,
                 Kokkos::Cuda::array_layout,
                 Kokkos::Device<Kokkos::Cuda, Kokkos::CudaSpace>,
                 Kokkos::MemoryTraits<Kokkos::Unmanaged>,
                 Kokkos::Impl::ViewDefault,
                 const double**,
                 Kokkos::LayoutLeft,
                 Kokkos::Device<Kokkos::Cuda, Kokkos::CudaSpace>,
                 Kokkos::MemoryTraits<Kokkos::Unmanaged>,
                 Kokkos::Impl::ViewDefault> {
  typedef double* RT;
  typedef Kokkos::Cuda::array_layout RL;
  typedef Kokkos::Device<Kokkos::Cuda, Kokkos::CudaSpace> RD;
  typedef Kokkos::MemoryTraits<Kokkos::Unmanaged> RM;
  typedef Kokkos::Impl::ViewDefault RS;
  typedef Kokkos::View<RT,RL,RD,RM,RS> RV;

  typedef const double** XT;
  typedef Kokkos::LayoutLeft XL;
  typedef Kokkos::Device<Kokkos::Cuda, Kokkos::CudaSpace> XD;
  typedef Kokkos::MemoryTraits<Kokkos::Unmanaged> XM;
  typedef Kokkos::Impl::ViewDefault XS;
  typedef Kokkos::View<XT,XL,XD,XM,XS> XMV;

  typedef XMV::execution_space execution_space;
  typedef XMV::size_type size_type;

  static void nrmInf (const RV& r, const XMV& X);
  static void nrmInf (const RV& r, const size_t r_col, const XMV& X, const size_t X_col);
};
#endif // KOKKOS_HAVE_CUDA

#ifdef KOKKOS_HAVE_CUDA
template<>
struct NrmInf_MV<double*,
                 Kokkos::Cuda::array_layout,
                 Kokkos::Device<Kokkos::Cuda, Kokkos::CudaUVMSpace>,
                 Kokkos::MemoryTraits<Kokkos::Unmanaged>,
                 Kokkos::Impl::ViewDefault,
                 const double**,
                 Kokkos::LayoutLeft,
                 Kokkos::Device<Kokkos::Cuda, Kokkos::CudaUVMSpace>,
                 Kokkos::MemoryTraits<Kokkos::Unmanaged>,
                 Kokkos::Impl::ViewDefault> {
  typedef double* RT;
  typedef Kokkos::Cuda::array_layout RL;
  typedef Kokkos::Device<Kokkos::Cuda, Kokkos::CudaUVMSpace> RD;
  typedef Kokkos::MemoryTraits<Kokkos::Unmanaged> RM;
  typedef Kokkos::Impl::ViewDefault RS;
  typedef Kokkos::View<RT,RL,RD,RM,RS> RV;

  typedef const double** XT;
  typedef Kokkos::LayoutLeft XL;
  typedef Kokkos::Device<Kokkos::Cuda, Kokkos::CudaUVMSpace> XD;
  typedef Kokkos::MemoryTraits<Kokkos::Unmanaged> XM;
  typedef Kokkos::Impl::ViewDefault XS;
  typedef Kokkos::View<XT,XL,XD,XM,XS> XMV;

  typedef XMV::execution_space execution_space;
  typedef XMV::size_type size_type;

  static void nrmInf (const RV& r, const XMV& X);
  static void nrmInf (const RV& r, const size_t r_col, const XMV& X, const size_t X_col);
};
#endif // KOKKOS_HAVE_CUDA

// Functor for multivectors X and Y and 1-D views a and b, that
// computes any of the following:
//
// 1. R(i,j) = alpha*X(i,j) + beta*Y(i,j) for alpha,beta in -1,0,1
// 2. R(i,j) = a(j)*X(i,j) + beta*Y(i,j) for beta in -1,0,1
// 3. R(i,j) = alpha*X(i,j) + beta*Y(i,j) for alpha in -1,0,1
// 4. R(i,j) = a(j)*X(i,j) + b(j)*Y(i,j)
//
// The template parameters scalar_x and scalar_y correspond to alpha
// resp. beta in the operation y = alpha*x + beta*y.  The values -1,
// 0, and -1 correspond to literal values of those coefficients.  The
// value 2 tells the functor to use the corresponding vector of
// coefficients.  Any literal coefficient of zero has BLAS semantics
// of ignoring the corresponding (multi)vector entry.  This does not
// apply to coefficients in the a and b vectors, if they are used.
template<class RV, class aVector, class XMV, class bVector, class YMV,
         int scalar_x, int scalar_y, class SizeType = typename RV::size_type>
struct MV_Axpby_Functor
{
  typedef typename RV::execution_space execution_space;
  typedef SizeType                           size_type;
  typedef Kokkos::Details::ArithTraits<typename RV::non_const_value_type> ATS;

  const size_type numCols;
  RV m_r;
  XMV m_x;
  YMV m_y;
  aVector m_a;
  bVector m_b;

  MV_Axpby_Functor (const RV& r, const XMV& x, const YMV& y,
                    const aVector& a, const bVector& b) :
    numCols (x.dimension_1 ()), m_r (r), m_x (x), m_y (y), m_a (a), m_b (b)
  {}

  KOKKOS_INLINE_FUNCTION
  void operator() (const size_type& i) const
  {
    // scalar_x and scalar_y are compile-time constants (since they
    // are template parameters), so the compiler should evaluate these
    // branches at compile time.
    if (scalar_x == 0 && scalar_y == 0) {
#ifdef KOKKOS_HAVE_PRAGMA_IVDEP
#pragma ivdep
#endif
#ifdef KOKKOS_HAVE_PRAGMA_VECTOR
#pragma vector always
#endif
      for (size_type k = 0; k < numCols; ++k) {
        m_r(i,k) = ATS::zero ();
      }
    }
    if (scalar_x == 0 && scalar_y == -1) {
#ifdef KOKKOS_HAVE_PRAGMA_IVDEP
#pragma ivdep
#endif
#ifdef KOKKOS_HAVE_PRAGMA_VECTOR
#pragma vector always
#endif
      for (size_type k = 0; k < numCols; ++k) {
        m_r(i,k) = -m_y(i,k);
      }
    }
    if (scalar_x == 0 && scalar_y == 1) {
#ifdef KOKKOS_HAVE_PRAGMA_IVDEP
#pragma ivdep
#endif
#ifdef KOKKOS_HAVE_PRAGMA_VECTOR
#pragma vector always
#endif
      for (size_type k = 0; k < numCols; ++k) {
        m_r(i,k) = m_y(i,k);
      }
    }
    if (scalar_x == 0 && scalar_y == 2) {
#ifdef KOKKOS_HAVE_PRAGMA_IVDEP
#pragma ivdep
#endif
#ifdef KOKKOS_HAVE_PRAGMA_VECTOR
#pragma vector always
#endif
      for (size_type k = 0; k < numCols; ++k) {
        m_r(i,k) = m_b(k)*m_y(i,k);
      }
    }
    if (scalar_x == -1 && scalar_y == 0) {
#ifdef KOKKOS_HAVE_PRAGMA_IVDEP
#pragma ivdep
#endif
#ifdef KOKKOS_HAVE_PRAGMA_VECTOR
#pragma vector always
#endif
      for (size_type k = 0; k < numCols; ++k) {
        m_r(i,k) = -m_x(i,k);
      }
    }
    if (scalar_x == -1 && scalar_y == -1) {
#ifdef KOKKOS_HAVE_PRAGMA_IVDEP
#pragma ivdep
#endif
#ifdef KOKKOS_HAVE_PRAGMA_VECTOR
#pragma vector always
#endif
      for (size_type k = 0; k < numCols; ++k) {
        m_r(i,k) = -m_x(i,k) - m_y(i,k);
      }
    }
    if (scalar_x == -1 && scalar_y == 1) {
#ifdef KOKKOS_HAVE_PRAGMA_IVDEP
#pragma ivdep
#endif
#ifdef KOKKOS_HAVE_PRAGMA_VECTOR
#pragma vector always
#endif
      for (size_type k = 0; k < numCols; ++k) {
        m_r(i,k) = -m_x(i,k) + m_y(i,k);
      }
    }
    if (scalar_x == -1 && scalar_y == 2) {
#ifdef KOKKOS_HAVE_PRAGMA_IVDEP
#pragma ivdep
#endif
#ifdef KOKKOS_HAVE_PRAGMA_VECTOR
#pragma vector always
#endif
      for (size_type k = 0; k < numCols; ++k) {
        m_r(i,k) = -m_x(i,k) + m_b(k)*m_y(i,k);
      }
    }
    if (scalar_x == 1 && scalar_y == 0) {
#ifdef KOKKOS_HAVE_PRAGMA_IVDEP
#pragma ivdep
#endif
#ifdef KOKKOS_HAVE_PRAGMA_VECTOR
#pragma vector always
#endif
      for (size_type k = 0; k < numCols; ++k) {
        m_r(i,k) = m_x(i,k);
      }
    }
    if (scalar_x == 1 && scalar_y == -1) {
#ifdef KOKKOS_HAVE_PRAGMA_IVDEP
#pragma ivdep
#endif
#ifdef KOKKOS_HAVE_PRAGMA_VECTOR
#pragma vector always
#endif
      for (size_type k = 0; k < numCols; ++k) {
        m_r(i,k) = m_x(i,k) - m_y(i,k);
      }
    }
    if (scalar_x == 1 && scalar_y == 1) {
#ifdef KOKKOS_HAVE_PRAGMA_IVDEP
#pragma ivdep
#endif
#ifdef KOKKOS_HAVE_PRAGMA_VECTOR
#pragma vector always
#endif
      for (size_type k = 0; k < numCols; ++k) {
        m_r(i,k) = m_x(i,k) + m_y(i,k);
      }
    }
    if (scalar_x == 1 && scalar_y == 2) {
#ifdef KOKKOS_HAVE_PRAGMA_IVDEP
#pragma ivdep
#endif
#ifdef KOKKOS_HAVE_PRAGMA_VECTOR
#pragma vector always
#endif
      for (size_type k = 0; k < numCols; ++k) {
        m_r(i,k) = m_x(i,k) + m_b(k)*m_y(i,k);
      }
    }
    if (scalar_x == 2 && scalar_y == 0) {
#ifdef KOKKOS_HAVE_PRAGMA_IVDEP
#pragma ivdep
#endif
#ifdef KOKKOS_HAVE_PRAGMA_VECTOR
#pragma vector always
#endif
      for (size_type k = 0; k < numCols; ++k) {
        m_r(i,k) = m_a(k)*m_x(i,k);
      }
    }
    if (scalar_x == 2 && scalar_y == -1) {
#ifdef KOKKOS_HAVE_PRAGMA_IVDEP
#pragma ivdep
#endif
#ifdef KOKKOS_HAVE_PRAGMA_VECTOR
#pragma vector always
#endif
      for (size_type k = 0; k < numCols; ++k) {
        m_r(i,k) = m_a(k)*m_x(i,k) - m_y(i,k);
      }
    }
    if (scalar_x == 2 && scalar_y == 1) {
#ifdef KOKKOS_HAVE_PRAGMA_IVDEP
#pragma ivdep
#endif
#ifdef KOKKOS_HAVE_PRAGMA_VECTOR
#pragma vector always
#endif
      for (size_type k = 0; k < numCols; ++k) {
        m_r(i,k) = m_a(k)*m_x(i,k) + m_y(i,k);
      }
    }
    if (scalar_x == 2 && scalar_y == 2) {
#ifdef KOKKOS_HAVE_PRAGMA_IVDEP
#pragma ivdep
#endif
#ifdef KOKKOS_HAVE_PRAGMA_VECTOR
#pragma vector always
#endif
      for (size_type k = 0; k < numCols; ++k) {
        m_r(i,k) = m_a(k)*m_x(i,k) + m_b(k)*m_y(i,k);
      }
    }
  }
};

// Variant of MV_Axpby_Functor, where a and b are scalars.
// This functor computes any of the following:
//
// 1. R(i,j) = alpha*X(i,j) + beta*Y(i,j) for alpha,beta in -1,0,1
// 2. R(i,j) = a*X(i,j) + beta*Y(i,j) for beta in -1,0,1
// 3. R(i,j) = alpha*X(i,j) + beta*Y(i,j) for alpha in -1,0,1
// 4. R(i,j) = a*X(i,j) + b*Y(i,j)
//
// The template parameters scalar_x and scalar_y correspond to alpha
// resp. beta in the operation y = alpha*x + beta*y.  The values -1,
// 0, and -1 correspond to literal values of those coefficients.  The
// value 2 tells the functor to use the corresponding vector of
// coefficients.  Any literal coefficient of zero has BLAS semantics
// of ignoring the corresponding (multi)vector entry.  This does not
// apply to coefficients in the a and b vectors, if they are used.
//
// This version works by partial specialization on aVector and bVector.
// In this partial specialization, both aVector and bVector are scalars.
template<class RV, class XMV, class YMV, int scalar_x, int scalar_y,
         class SizeType>
struct MV_Axpby_Functor<RV, typename XMV::non_const_value_type, XMV,
                        typename YMV::non_const_value_type, YMV,
                        scalar_x, scalar_y, SizeType>
{
  typedef typename RV::execution_space execution_space;
  typedef SizeType                           size_type;
  typedef Kokkos::Details::ArithTraits<typename RV::non_const_value_type> ATS;

  const size_type numCols;
  RV m_r;
  XMV m_x;
  YMV m_y;
  const typename XMV::non_const_value_type m_a;
  const typename YMV::non_const_value_type m_b;

  MV_Axpby_Functor (const RV& r, const XMV& x, const YMV& y,
                  const typename XMV::non_const_value_type& a,
                  const typename YMV::non_const_value_type& b) :
    numCols (x.dimension_1 ()), m_r (r), m_x (x), m_y (y), m_a (a), m_b (b)
  {}

  KOKKOS_INLINE_FUNCTION
  void operator() (const size_type& i) const
  {
    // scalar_x and scalar_y are compile-time constants (since they
    // are template parameters), so the compiler should evaluate these
    // branches at compile time.
    if (scalar_x == 0 && scalar_y == 0) {
#ifdef KOKKOS_HAVE_PRAGMA_IVDEP
#pragma ivdep
#endif
#ifdef KOKKOS_HAVE_PRAGMA_VECTOR
#pragma vector always
#endif
      for (size_type k = 0; k < numCols; ++k) {
        m_r(i,k) = ATS::zero ();
      }
    }
    if (scalar_x == 0 && scalar_y == -1) {
#ifdef KOKKOS_HAVE_PRAGMA_IVDEP
#pragma ivdep
#endif
#ifdef KOKKOS_HAVE_PRAGMA_VECTOR
#pragma vector always
#endif
      for (size_type k = 0; k < numCols; ++k) {
        m_r(i,k) = -m_y(i,k);
      }
    }
    if (scalar_x == 0 && scalar_y == 1) {
#ifdef KOKKOS_HAVE_PRAGMA_IVDEP
#pragma ivdep
#endif
#ifdef KOKKOS_HAVE_PRAGMA_VECTOR
#pragma vector always
#endif
      for (size_type k = 0; k < numCols; ++k) {
        m_r(i,k) = m_y(i,k);
      }
    }
    if (scalar_x == 0 && scalar_y == 2) {
#ifdef KOKKOS_HAVE_PRAGMA_IVDEP
#pragma ivdep
#endif
#ifdef KOKKOS_HAVE_PRAGMA_VECTOR
#pragma vector always
#endif
      for (size_type k = 0; k < numCols; ++k) {
        m_r(i,k) = m_b*m_y(i,k);
      }
    }
    if (scalar_x == -1 && scalar_y == 0) {
#ifdef KOKKOS_HAVE_PRAGMA_IVDEP
#pragma ivdep
#endif
#ifdef KOKKOS_HAVE_PRAGMA_VECTOR
#pragma vector always
#endif
      for (size_type k = 0; k < numCols; ++k) {
        m_r(i,k) = -m_x(i,k);
      }
    }
    if (scalar_x == -1 && scalar_y == -1) {
#ifdef KOKKOS_HAVE_PRAGMA_IVDEP
#pragma ivdep
#endif
#ifdef KOKKOS_HAVE_PRAGMA_VECTOR
#pragma vector always
#endif
      for (size_type k = 0; k < numCols; ++k) {
        m_r(i,k) = -m_x(i,k) - m_y(i,k);
      }
    }
    if (scalar_x == -1 && scalar_y == 1) {
#ifdef KOKKOS_HAVE_PRAGMA_IVDEP
#pragma ivdep
#endif
#ifdef KOKKOS_HAVE_PRAGMA_VECTOR
#pragma vector always
#endif
      for (size_type k = 0; k < numCols; ++k) {
        m_r(i,k) = -m_x(i,k) + m_y(i,k);
      }
    }
    if (scalar_x == -1 && scalar_y == 2) {
#ifdef KOKKOS_HAVE_PRAGMA_IVDEP
#pragma ivdep
#endif
#ifdef KOKKOS_HAVE_PRAGMA_VECTOR
#pragma vector always
#endif
      for (size_type k = 0; k < numCols; ++k) {
        m_r(i,k) = -m_x(i,k) + m_b*m_y(i,k);
      }
    }
    if (scalar_x == 1 && scalar_y == 0) {
#ifdef KOKKOS_HAVE_PRAGMA_IVDEP
#pragma ivdep
#endif
#ifdef KOKKOS_HAVE_PRAGMA_VECTOR
#pragma vector always
#endif
      for (size_type k = 0; k < numCols; ++k) {
        m_r(i,k) = m_x(i,k);
      }
    }
    if (scalar_x == 1 && scalar_y == -1) {
#ifdef KOKKOS_HAVE_PRAGMA_IVDEP
#pragma ivdep
#endif
#ifdef KOKKOS_HAVE_PRAGMA_VECTOR
#pragma vector always
#endif
      for (size_type k = 0; k < numCols; ++k) {
        m_r(i,k) = m_x(i,k) - m_y(i,k);
      }
    }
    if (scalar_x == 1 && scalar_y == 1) {
#ifdef KOKKOS_HAVE_PRAGMA_IVDEP
#pragma ivdep
#endif
#ifdef KOKKOS_HAVE_PRAGMA_VECTOR
#pragma vector always
#endif
      for (size_type k = 0; k < numCols; ++k) {
        m_r(i,k) = m_x(i,k) + m_y(i,k);
      }
    }
    if (scalar_x == 1 && scalar_y == 2) {
#ifdef KOKKOS_HAVE_PRAGMA_IVDEP
#pragma ivdep
#endif
#ifdef KOKKOS_HAVE_PRAGMA_VECTOR
#pragma vector always
#endif
      for (size_type k = 0; k < numCols; ++k) {
        m_r(i,k) = m_x(i,k) + m_b*m_y(i,k);
      }
    }
    if (scalar_x == 2 && scalar_y == 0) {
#ifdef KOKKOS_HAVE_PRAGMA_IVDEP
#pragma ivdep
#endif
#ifdef KOKKOS_HAVE_PRAGMA_VECTOR
#pragma vector always
#endif
      for (size_type k = 0; k < numCols; ++k) {
        m_r(i,k) = m_a*m_x(i,k);
      }
    }
    if (scalar_x == 2 && scalar_y == -1) {
#ifdef KOKKOS_HAVE_PRAGMA_IVDEP
#pragma ivdep
#endif
#ifdef KOKKOS_HAVE_PRAGMA_VECTOR
#pragma vector always
#endif
      for (size_type k = 0; k < numCols; ++k) {
        m_r(i,k) = m_a*m_x(i,k) - m_y(i,k);
      }
    }
    if (scalar_x == 2 && scalar_y == 1) {
#ifdef KOKKOS_HAVE_PRAGMA_IVDEP
#pragma ivdep
#endif
#ifdef KOKKOS_HAVE_PRAGMA_VECTOR
#pragma vector always
#endif
      for (size_type k = 0; k < numCols; ++k) {
        m_r(i,k) = m_a*m_x(i,k) + m_y(i,k);
      }
    }
    if (scalar_x == 2 && scalar_y == 2) {
#ifdef KOKKOS_HAVE_PRAGMA_IVDEP
#pragma ivdep
#endif
#ifdef KOKKOS_HAVE_PRAGMA_VECTOR
#pragma vector always
#endif
      for (size_type k = 0; k < numCols; ++k) {
        m_r(i,k) = m_a*m_x(i,k) + m_b*m_y(i,k);
      }
    }
  }
};


// Column-unrolled variant of MV_Axpby_Functor.  The number of columns
// in X and Y, UNROLL, is a compile-time constant.
template<class RV, class aVector, class XMV, class bVector, class YMV,
         int scalar_x, int scalar_y, int UNROLL, class SizeType>
struct MV_Axpby_Unroll_Functor
{
  typedef typename RV::execution_space execution_space;
  typedef SizeType                           size_type;
  typedef Kokkos::Details::ArithTraits<typename RV::non_const_value_type> ATS;

  RV m_r;
  XMV m_x;
  YMV m_y;
  aVector m_a;
  bVector m_b;

  MV_Axpby_Unroll_Functor (const RV& r, const XMV& x, const YMV& y,
                           const aVector& a, const bVector& b) :
    m_r (r), m_x (x), m_y (y), m_a (a), m_b (b)
  {}

  KOKKOS_INLINE_FUNCTION
  void operator() (const size_type& i) const
  {
    // scalar_x and scalar_y are compile-time constants (since they
    // are template parameters), so the compiler should evaluate these
    // branches at compile time.
    if (scalar_x == 0 && scalar_y == 0) {
#ifdef KOKKOS_HAVE_PRAGMA_UNROLL
#pragma unroll
#endif
      for (int k = 0; k < UNROLL; ++k) {
        m_r(i,k) = ATS::zero ();
      }
    }
    if (scalar_x == 0 && scalar_y == -1) {
#ifdef KOKKOS_HAVE_PRAGMA_UNROLL
#pragma unroll
#endif
      for (int k = 0; k < UNROLL; ++k) {
        m_r(i,k) = -m_y(i,k);
      }
    }
    if (scalar_x == 0 && scalar_y == 1) {
#ifdef KOKKOS_HAVE_PRAGMA_UNROLL
#pragma unroll
#endif
      for (int k = 0; k < UNROLL; ++k) {
        m_r(i,k) = m_y(i,k);
      }
    }
    if (scalar_x == 0 && scalar_y == 2) {
#ifdef KOKKOS_HAVE_PRAGMA_UNROLL
#pragma unroll
#endif
      for (int k = 0; k < UNROLL; ++k) {
        m_r(i,k) = m_b(k)*m_y(i,k);
      }
    }
    if (scalar_x == -1 && scalar_y == 0) {
#ifdef KOKKOS_HAVE_PRAGMA_UNROLL
#pragma unroll
#endif
      for (int k = 0; k < UNROLL; ++k) {
        m_r(i,k) = -m_x(i,k);
      }
    }
    if (scalar_x == -1 && scalar_y == -1) {
#ifdef KOKKOS_HAVE_PRAGMA_UNROLL
#pragma unroll
#endif
      for (int k = 0; k < UNROLL; ++k) {
        m_r(i,k) = -m_x(i,k) - m_y(i,k);
      }
    }
    if (scalar_x == -1 && scalar_y == 1) {
#ifdef KOKKOS_HAVE_PRAGMA_UNROLL
#pragma unroll
#endif
      for (int k = 0; k < UNROLL; ++k) {
        m_r(i,k) = -m_x(i,k) + m_y(i,k);
      }
    }
    if (scalar_x == -1 && scalar_y == 2) {
#ifdef KOKKOS_HAVE_PRAGMA_UNROLL
#pragma unroll
#endif
      for (int k = 0; k < UNROLL; ++k) {
        m_r(i,k) = -m_x(i,k) + m_b(k)*m_y(i,k);
      }
    }
    if (scalar_x == 1 && scalar_y == 0) {
#ifdef KOKKOS_HAVE_PRAGMA_UNROLL
#pragma unroll
#endif
      for (int k = 0; k < UNROLL; ++k) {
        m_r(i,k) = m_x(i,k);
      }
    }
    if (scalar_x == 1 && scalar_y == -1) {
#ifdef KOKKOS_HAVE_PRAGMA_UNROLL
#pragma unroll
#endif
      for (int k = 0; k < UNROLL; ++k) {
        m_r(i,k) = m_x(i,k) - m_y(i,k);
      }
    }
    if (scalar_x == 1 && scalar_y == 1) {
#ifdef KOKKOS_HAVE_PRAGMA_UNROLL
#pragma unroll
#endif
      for (int k = 0; k < UNROLL; ++k) {
        m_r(i,k) = m_x(i,k) + m_y(i,k);
      }
    }
    if (scalar_x == 1 && scalar_y == 2) {
#ifdef KOKKOS_HAVE_PRAGMA_UNROLL
#pragma unroll
#endif
      for (int k = 0; k < UNROLL; ++k) {
        m_r(i,k) = m_x(i,k) + m_b(k)*m_y(i,k);
      }
    }
    if (scalar_x == 2 && scalar_y == 0) {
#ifdef KOKKOS_HAVE_PRAGMA_UNROLL
#pragma unroll
#endif
      for (int k = 0; k < UNROLL; ++k) {
        m_r(i,k) = m_a(k)*m_x(i,k);
      }
    }
    if (scalar_x == 2 && scalar_y == -1) {
#ifdef KOKKOS_HAVE_PRAGMA_UNROLL
#pragma unroll
#endif
      for (int k = 0; k < UNROLL; ++k) {
        m_r(i,k) = m_a(k)*m_x(i,k) - m_y(i,k);
      }
    }
    if (scalar_x == 2 && scalar_y == 1) {
#ifdef KOKKOS_HAVE_PRAGMA_UNROLL
#pragma unroll
#endif
      for (int k = 0; k < UNROLL; ++k) {
        m_r(i,k) = m_a(k)*m_x(i,k) + m_y(i,k);
      }
    }
    if (scalar_x == 2 && scalar_y == 2) {
#ifdef KOKKOS_HAVE_PRAGMA_UNROLL
#pragma unroll
#endif
      for (int k = 0; k < UNROLL; ++k) {
        m_r(i,k) = m_a(k)*m_x(i,k) + m_b(k)*m_y(i,k);
      }
    }
  }
};


// Variant of MV_Axpby_Unroll_Functor for single coefficients (rather
// than vectors of coefficients) a and b.  The number of columns in X
// and Y, UNROLL, is a compile-time constant.
template<class RV, class XMV, class YMV,
         int scalar_x, int scalar_y, int UNROLL, class SizeType>
struct MV_Axpby_Unroll_Functor<RV, typename XMV::non_const_value_type, XMV,
                               typename YMV::non_const_value_type, YMV,
                               scalar_x, scalar_y, UNROLL, SizeType>
{
  typedef typename RV::execution_space execution_space;
  typedef SizeType                           size_type;
  typedef Kokkos::Details::ArithTraits<typename RV::non_const_value_type> ATS;

  RV m_r;
  XMV m_x;
  YMV m_y;
  const typename XMV::non_const_value_type m_a;
  const typename YMV::non_const_value_type m_b;

  MV_Axpby_Unroll_Functor (const RV& r, const XMV& x, const YMV& y,
                           const typename XMV::non_const_value_type& a,
                           const typename YMV::non_const_value_type& b) :
    m_r (r), m_x (x), m_y (y), m_a (a), m_b (b)
  {}

  KOKKOS_INLINE_FUNCTION
  void operator() (const size_type& i) const
  {
    // scalar_x and scalar_y are compile-time constants (since they
    // are template parameters), so the compiler should evaluate these
    // branches at compile time.
    if (scalar_x == 0 && scalar_y == 0) {
#ifdef KOKKOS_HAVE_PRAGMA_UNROLL
#pragma unroll
#endif
      for (int k = 0; k < UNROLL; ++k) {
        m_r(i,k) = ATS::zero ();
      }
    }
    if (scalar_x == 0 && scalar_y == -1) {
#ifdef KOKKOS_HAVE_PRAGMA_UNROLL
#pragma unroll
#endif
      for (int k = 0; k < UNROLL; ++k) {
        m_r(i,k) = -m_y(i,k);
      }
    }
    if (scalar_x == 0 && scalar_y == 1) {
#ifdef KOKKOS_HAVE_PRAGMA_UNROLL
#pragma unroll
#endif
      for (int k = 0; k < UNROLL; ++k) {
        m_r(i,k) = m_y(i,k);
      }
    }
    if (scalar_x == 0 && scalar_y == 2) {
#ifdef KOKKOS_HAVE_PRAGMA_UNROLL
#pragma unroll
#endif
      for (int k = 0; k < UNROLL; ++k) {
        m_r(i,k) = m_b*m_y(i,k);
      }
    }
    if (scalar_x == -1 && scalar_y == 0) {
#ifdef KOKKOS_HAVE_PRAGMA_UNROLL
#pragma unroll
#endif
      for (int k = 0; k < UNROLL; ++k) {
        m_r(i,k) = -m_x(i,k);
      }
    }
    if (scalar_x == -1 && scalar_y == -1) {
#ifdef KOKKOS_HAVE_PRAGMA_UNROLL
#pragma unroll
#endif
      for (int k = 0; k < UNROLL; ++k) {
        m_r(i,k) = -m_x(i,k) - m_y(i,k);
      }
    }
    if (scalar_x == -1 && scalar_y == 1) {
#ifdef KOKKOS_HAVE_PRAGMA_UNROLL
#pragma unroll
#endif
      for (int k = 0; k < UNROLL; ++k) {
        m_r(i,k) = -m_x(i,k) + m_y(i,k);
      }
    }
    if (scalar_x == -1 && scalar_y == 2) {
#ifdef KOKKOS_HAVE_PRAGMA_UNROLL
#pragma unroll
#endif
      for (int k = 0; k < UNROLL; ++k) {
        m_r(i,k) = -m_x(i,k) + m_b*m_y(i,k);
      }
    }
    if (scalar_x == 1 && scalar_y == 0) {
#ifdef KOKKOS_HAVE_PRAGMA_UNROLL
#pragma unroll
#endif
      for (int k = 0; k < UNROLL; ++k) {
        m_r(i,k) = m_x(i,k);
      }
    }
    if (scalar_x == 1 && scalar_y == -1) {
#ifdef KOKKOS_HAVE_PRAGMA_UNROLL
#pragma unroll
#endif
      for (int k = 0; k < UNROLL; ++k) {
        m_r(i,k) = m_x(i,k) - m_y(i,k);
      }
    }
    if (scalar_x == 1 && scalar_y == 1) {
#ifdef KOKKOS_HAVE_PRAGMA_UNROLL
#pragma unroll
#endif
      for (int k = 0; k < UNROLL; ++k) {
        m_r(i,k) = m_x(i,k) + m_y(i,k);
      }
    }
    if (scalar_x == 1 && scalar_y == 2) {
#ifdef KOKKOS_HAVE_PRAGMA_UNROLL
#pragma unroll
#endif
      for (int k = 0; k < UNROLL; ++k) {
        m_r(i,k) = m_x(i,k) + m_b*m_y(i,k);
      }
    }
    if (scalar_x == 2 && scalar_y == 0) {
#ifdef KOKKOS_HAVE_PRAGMA_UNROLL
#pragma unroll
#endif
      for (int k = 0; k < UNROLL; ++k) {
        m_r(i,k) = m_a*m_x(i,k);
      }
    }
    if (scalar_x == 2 && scalar_y == -1) {
#ifdef KOKKOS_HAVE_PRAGMA_UNROLL
#pragma unroll
#endif
      for (int k = 0; k < UNROLL; ++k) {
        m_r(i,k) = m_a*m_x(i,k) - m_y(i,k);
      }
    }
    if (scalar_x == 2 && scalar_y == 1) {
#ifdef KOKKOS_HAVE_PRAGMA_UNROLL
#pragma unroll
#endif
      for (int k = 0; k < UNROLL; ++k) {
        m_r(i,k) = m_a*m_x(i,k) + m_y(i,k);
      }
    }
    if (scalar_x == 2 && scalar_y == 2) {
#ifdef KOKKOS_HAVE_PRAGMA_UNROLL
#pragma unroll
#endif
      for (int k = 0; k < UNROLL; ++k) {
        m_r(i,k) = m_a*m_x(i,k) + m_b*m_y(i,k);
      }
    }
  }
};

// Single-vector version of MV_Axpby_Functor.  By default, a and b are
// still 1-D Views.  Below is a partial specialization that lets both
// of them be scalars.  This functor computes any of the following:
//
// 1. Y(i) = alpha*X(i) + beta*Y(i) for alpha,beta in -1,0,1
// 2. Y(i) = a(0)*X(i) + beta*Y(i) for beta in -1,0,1
// 3. Y(i) = alpha*X(i) + b(0)*Y(i) for alpha in -1,0,1
// 4. Y(i) = a(0)*X(i) + b(0)*Y(i)
//
// The template parameters scalar_x and scalar_y correspond to alpha
// resp. beta in the operation y = alpha*x + beta*y.  The values -1,
// 0, and -1 correspond to literal values of those coefficients.  The
// value 2 tells the functor to use the corresponding vector of
// coefficients.  Any literal coefficient of zero has BLAS semantics
// of ignoring the corresponding (multi)vector entry.  This does not
// apply to coefficients in the a and b vectors, if they are used.
template<class RV, class AV, class XV, class BV, class YV,
         int scalar_x, int scalar_y, class SizeType>
struct V_Axpby_Functor {
  typedef typename RV::execution_space execution_space;
  typedef SizeType                           size_type;
  typedef Kokkos::Details::ArithTraits<typename RV::non_const_value_type> ATS;

  RV m_r;
  XV m_x;
  YV m_y;
  AV m_a;
  BV m_b;

  V_Axpby_Functor (const RV& r, const XV& x, const YV& y,
                   const AV& a, const BV& b) :
    m_r (r), m_x (x), m_y (y), m_a (a), m_b (b)
  {}

  KOKKOS_INLINE_FUNCTION
  void operator() (const size_type& i) const
  {
    // scalar_x and scalar_y are compile-time constants (since they
    // are template parameters), so the compiler should evaluate these
    // branches at compile time.
    if (scalar_x == 0 && scalar_y == 0) {
      m_r(i) = ATS::zero ();
    }
    if (scalar_x == 0 && scalar_y == -1) {
      m_r(i) = -m_y(i);
    }
    if (scalar_x == 0 && scalar_y == 1) {
      m_r(i) = m_y(i);
    }
    if (scalar_x == 0 && scalar_y == 2) {
      m_r(i) = m_b(0)*m_y(i);
    }
    if (scalar_x == -1 && scalar_y == 0) {
      m_r(i) = -m_x(i);
    }
    if (scalar_x == -1 && scalar_y == -1) {
      m_r(i) = -m_x(i) - m_y(i);
    }
    if (scalar_x == -1 && scalar_y == 1) {
      m_r(i) = -m_x(i) + m_y(i);
    }
    if (scalar_x == -1 && scalar_y == 2) {
      m_r(i) = -m_x(i) + m_b(0)*m_y(i);
    }
    if (scalar_x == 1 && scalar_y == 0) {
      m_r(i) = m_x(i);
    }
    if (scalar_x == 1 && scalar_y == -1) {
      m_r(i) = m_x(i) - m_y(i);
    }
    if (scalar_x == 1 && scalar_y == 1) {
      m_r(i) = m_x(i) + m_y(i);
    }
    if (scalar_x == 1 && scalar_y == 2) {
      m_r(i) = m_x(i) + m_b(0)*m_y(i);
    }
    if (scalar_x == 2 && scalar_y == 0) {
      m_r(i) = m_a(0)*m_x(i);
    }
    if (scalar_x == 2 && scalar_y == -1) {
      m_r(i) = m_a(0)*m_x(i) - m_y(i);
    }
    if (scalar_x == 2 && scalar_y == 1) {
      m_r(i) = m_a(0)*m_x(i) + m_y(i);
    }
    if (scalar_x == 2 && scalar_y == 2) {
      m_r(i) = m_a(0)*m_x(i) + m_b(0)*m_y(i);
    }
  }
};


// Partial specialization of V_Axpby_Functor that lets a and b be
// scalars (rather than 1-D Views, as in the most general version
// above).  This functor computes any of the following:
//
// 1. Y(i) = alpha*X(i) + beta*Y(i) for alpha,beta in -1,0,1
// 2. Y(i) = a*X(i) + beta*Y(i) for beta in -1,0,1
// 3. Y(i) = alpha*X(i) + b*Y(i) for alpha in -1,0,1
// 4. Y(i) = a*X(i) + b*Y(i)
//
// The template parameters scalar_x and scalar_y correspond to alpha
// resp. beta in the operation y = alpha*x + beta*y.  The values -1,
// 0, and -1 correspond to literal values of those coefficients.  The
// value 2 tells the functor to use the corresponding vector of
// coefficients.  Any literal coefficient of zero has BLAS semantics
// of ignoring the corresponding (multi)vector entry.  This does not
// apply to coefficients in the a and b vectors, if they are used.
template<class RV, class XV, class YV,
         int scalar_x, int scalar_y, class SizeType>
struct V_Axpby_Functor<RV, typename XV::non_const_value_type, XV,
                       typename YV::non_const_value_type, YV,
                       scalar_x, scalar_y, SizeType> {
  typedef typename RV::execution_space execution_space;
  typedef SizeType                           size_type;
  typedef Kokkos::Details::ArithTraits<typename RV::non_const_value_type> ATS;

  RV m_r;
  XV m_x;
  YV m_y;
  const typename XV::non_const_value_type m_a;
  const typename YV::non_const_value_type m_b;

  V_Axpby_Functor (const RV& r, const XV& x, const YV& y,
                   const typename XV::non_const_value_type& a,
                   const typename YV::non_const_value_type& b) :
    m_r (r), m_x (x), m_y (y), m_a (a), m_b (b)
  {}

  KOKKOS_INLINE_FUNCTION
  void operator() (const size_type& i) const
  {
    // scalar_x and scalar_y are compile-time constants (since they
    // are template parameters), so the compiler should evaluate these
    // branches at compile time.
    if (scalar_x == 0 && scalar_y == 0) {
      m_r(i) = ATS::zero ();
    }
    if (scalar_x == 0 && scalar_y == -1) {
      m_r(i) = -m_y(i);
    }
    if (scalar_x == 0 && scalar_y == 1) {
      m_r(i) = m_y(i);
    }
    if (scalar_x == 0 && scalar_y == 2) {
      m_r(i) = m_b*m_y(i);
    }
    if (scalar_x == -1 && scalar_y == 0) {
      m_r(i) = -m_x(i);
    }
    if (scalar_x == -1 && scalar_y == -1) {
      m_r(i) = -m_x(i) - m_y(i);
    }
    if (scalar_x == -1 && scalar_y == 1) {
      m_r(i) = -m_x(i) + m_y(i);
    }
    if (scalar_x == -1 && scalar_y == 2) {
      m_r(i) = -m_x(i) + m_b*m_y(i);
    }
    if (scalar_x == 1 && scalar_y == 0) {
      m_r(i) = m_x(i);
    }
    if (scalar_x == 1 && scalar_y == -1) {
      m_r(i) = m_x(i) - m_y(i);
    }
    if (scalar_x == 1 && scalar_y == 1) {
      m_r(i) = m_x(i) + m_y(i);
    }
    if (scalar_x == 1 && scalar_y == 2) {
      m_r(i) = m_x(i) + m_b*m_y(i);
    }
    if (scalar_x == 2 && scalar_y == 0) {
      m_r(i) = m_a*m_x(i);
    }
    if (scalar_x == 2 && scalar_y == -1) {
      m_r(i) = m_a*m_x(i) - m_y(i);
    }
    if (scalar_x == 2 && scalar_y == 1) {
      m_r(i) = m_a*m_x(i) + m_y(i);
    }
    if (scalar_x == 2 && scalar_y == 2) {
      m_r(i) = m_a*m_x(i) + m_b*m_y(i);
    }
  }
};

// Invoke the unrolled multivector functor that computes any of the
// following:
//
// 1. R(i,j) = a*X(i,j) + b*Y(i,j) for a,b in -1,0,1
// 2. R(i,j) = av(j)*X(i,j) + b*Y(i,j) for b in -1,0,1
// 3. R(i,j) = a*X(i,j) + b*Y(i,j) for a in -1,0,1
// 4. R(i,j) = av(j)*X(i,j) + bv(j)*Y(i,j)
//
// a and b come in as integers.  The values -1, 0, and 1 correspond to
// the literal values of the coefficients.  The value 2 tells the
// functor to use the corresponding vector of coefficients: a == 2
// means use av, and b == 2 means use bv.  Otherwise, av resp. vb are
// ignored.
//
// Any literal coefficient of zero has BLAS semantics of ignoring the
// corresponding (multi)vector entry.  This does NOT apply to
// coefficients in av and bv vectors, if they are used.
//
// Either av and bv are both 1-D Views, or av and bv are both scalars.
template<class RMV, class aVector, class XMV,
         class bVector, class YMV, int UNROLL, class SizeType>
void
MV_Axpby_Unrolled (const RMV& r, const aVector& av, const XMV& x,
                   const bVector& bv, const YMV& y,
                   int a = 2, int b = 2)
{
  typedef typename XMV::execution_space execution_space;

  if (a == 0 && b == 0) {
    MV_Axpby_Unroll_Functor<RMV, aVector, XMV, bVector, YMV, 0, 0, UNROLL, SizeType> op (r, x, y, av, bv);
    const SizeType numRows = x.dimension_0 ();
    Kokkos::RangePolicy<execution_space, SizeType> policy (0, numRows);
    Kokkos::parallel_for (policy, op);
    return;
  }
  if (a == 0 && b == -1) {
    MV_Axpby_Unroll_Functor<RMV, aVector, XMV, bVector, YMV, 0, -1, UNROLL, SizeType> op (r, x, y, av, bv);
    const SizeType numRows = x.dimension_0 ();
    Kokkos::RangePolicy<execution_space, SizeType> policy (0, numRows);
    Kokkos::parallel_for (policy, op);
    return;
  }
  if (a == 0 && b == 1) {
    MV_Axpby_Unroll_Functor<RMV, aVector, XMV, bVector, YMV, 0, 1, UNROLL, SizeType> op (r, x, y, av, bv);
    const SizeType numRows = x.dimension_0 ();
    Kokkos::RangePolicy<execution_space, SizeType> policy (0, numRows);
    Kokkos::parallel_for (policy, op);
    return;
  }
  if (a == 0 && b == 2) {
    MV_Axpby_Unroll_Functor<RMV, aVector, XMV, bVector, YMV, 0, 2, UNROLL, SizeType> op (r, x, y, av, bv);
    const SizeType numRows = x.dimension_0 ();
    Kokkos::RangePolicy<execution_space, SizeType> policy (0, numRows);
    Kokkos::parallel_for (policy, op);
    return;
  }
  // a == -1
  if (a == -1 && b == 0) {
    MV_Axpby_Unroll_Functor<RMV, aVector, XMV, bVector, YMV, -1, 0, UNROLL, SizeType> op (r, x, y, av, bv);
    const SizeType numRows = x.dimension_0 ();
    Kokkos::RangePolicy<execution_space, SizeType> policy (0, numRows);
    Kokkos::parallel_for (policy, op);
    return;
  }
  if (a == -1 && b == -1) {
    MV_Axpby_Unroll_Functor<RMV, aVector, XMV, bVector, YMV, -1, -1, UNROLL, SizeType> op (r, x, y, av, bv);
    const SizeType numRows = x.dimension_0 ();
    Kokkos::RangePolicy<execution_space, SizeType> policy (0, numRows);
    Kokkos::parallel_for (policy, op);
    return;
  }
  if (a == -1 && b == 1) {
    MV_Axpby_Unroll_Functor<RMV, aVector, XMV, bVector, YMV, -1, 1, UNROLL, SizeType> op (r, x, y, av, bv);
    const SizeType numRows = x.dimension_0 ();
    Kokkos::RangePolicy<execution_space, SizeType> policy (0, numRows);
    Kokkos::parallel_for (policy, op);
    return;
  }
  if (a == -1 && b == 2) {
    MV_Axpby_Unroll_Functor<RMV, aVector, XMV, bVector, YMV, -1, 2, UNROLL, SizeType> op (r, x, y, av, bv);
    const SizeType numRows = x.dimension_0 ();
    Kokkos::RangePolicy<execution_space, SizeType> policy (0, numRows);
    Kokkos::parallel_for (policy, op);
    return;
  }
  // a == 1
  if (a == 1 && b == 0) {
    MV_Axpby_Unroll_Functor<RMV, aVector, XMV, bVector, YMV, 1, 0, UNROLL, SizeType> op (r, x, y, av, bv);
    const SizeType numRows = x.dimension_0 ();
    Kokkos::RangePolicy<execution_space, SizeType> policy (0, numRows);
    Kokkos::parallel_for (policy, op);
    return;
  }
  if (a == 1 && b == -1) {
    MV_Axpby_Unroll_Functor<RMV, aVector, XMV, bVector, YMV, 1, -1, UNROLL, SizeType> op (r, x, y, av, bv);
    const SizeType numRows = x.dimension_0 ();
    Kokkos::RangePolicy<execution_space, SizeType> policy (0, numRows);
    Kokkos::parallel_for (policy, op);
    return;
  }
  if (a == 1 && b == 1) {
    MV_Axpby_Unroll_Functor<RMV, aVector, XMV, bVector, YMV, 1, 1, UNROLL, SizeType> op (r, x, y, av, bv);
    const SizeType numRows = x.dimension_0 ();
    Kokkos::RangePolicy<execution_space, SizeType> policy (0, numRows);
    Kokkos::parallel_for (policy, op);
    return;
  }
  if (a == 1 && b == 2) {
    MV_Axpby_Unroll_Functor<RMV, aVector, XMV, bVector, YMV, 1, 2, UNROLL, SizeType> op (r, x, y, av, bv);
    const SizeType numRows = x.dimension_0 ();
    Kokkos::RangePolicy<execution_space, SizeType> policy (0, numRows);
    Kokkos::parallel_for (policy, op);
    return;
  }
  // a and b arbitrary (not -1, 0, or 1)

  MV_Axpby_Unroll_Functor<RMV, aVector, XMV, bVector, YMV, 2, 2, UNROLL, SizeType> op (r, x, y, av, bv);
  const SizeType numRows = x.dimension_0 ();
  Kokkos::RangePolicy<execution_space, SizeType> policy (0, numRows);
  Kokkos::parallel_for (policy, op);
}


// Invoke the "generic" (not unrolled) multivector functor that
// computes any of the following:
//
// 1. R(i,j) = a*X(i,j) + b*Y(i,j) for a,b in -1,0,1
// 2. R(i,j) = av(j)*X(i,j) + b*Y(i,j) for b in -1,0,1
// 3. R(i,j) = a*X(i,j) + b*Y(i,j) for a in -1,0,1
// 4. R(i,j) = av(j)*X(i,j) + bv(j)*Y(i,j)
//
// a and b come in as integers.  The values -1, 0, and 1 correspond to
// the literal values of the coefficients.  The value 2 tells the
// functor to use the corresponding vector of coefficients: a == 2
// means use av, and b == 2 means use bv.  Otherwise, av resp. vb are
// ignored.
//
// Any literal coefficient of zero has BLAS semantics of ignoring the
// corresponding (multi)vector entry.  This does NOT apply to
// coefficients in av and bv vectors, if they are used.
//
// Either av and bv are both 1-D Views, or av and bv are both scalars.
template<class RVector, class aVector, class XVector,
         class bVector, class YVector, class SizeType>
void
MV_Axpby_Generic (const RVector& r, const aVector& av, const XVector& x,
                  const bVector& bv, const YVector& y,
                  int a = 2, int b = 2)
{
  typedef typename XVector::execution_space execution_space;

  if (a == 0 && b == 0) {
    MV_Axpby_Functor<RVector, aVector, XVector, bVector, YVector, 0, 0, SizeType> op (r, x, y, av, bv);
    const SizeType numRows = x.dimension_0 ();
    Kokkos::RangePolicy<execution_space, SizeType> policy (0, numRows);
    Kokkos::parallel_for (policy, op);
    return;
  }
  if (a == 0 && b == -1) {
    MV_Axpby_Functor<RVector, aVector, XVector, bVector, YVector, 0, -1, SizeType> op (r, x, y, av, bv);
    const SizeType numRows = x.dimension_0 ();
    Kokkos::RangePolicy<execution_space, SizeType> policy (0, numRows);
    Kokkos::parallel_for (policy, op);
    return;
  }
  if (a == 0 && b == 1) {
    MV_Axpby_Functor<RVector, aVector, XVector, bVector, YVector, 0, 1, SizeType> op (r, x, y, av, bv);
    const SizeType numRows = x.dimension_0 ();
    Kokkos::RangePolicy<execution_space, SizeType> policy (0, numRows);
    Kokkos::parallel_for (policy, op);
    return;
  }
  if (a == 0 && b == 2) {
    MV_Axpby_Functor<RVector, aVector, XVector, bVector, YVector, 0, 2, SizeType> op (r, x, y, av, bv);
    const SizeType numRows = x.dimension_0 ();
    Kokkos::RangePolicy<execution_space, SizeType> policy (0, numRows);
    Kokkos::parallel_for (policy, op);
    return;
  }
  // a == -1
  if (a == -1 && b == 0) {
    MV_Axpby_Functor<RVector, aVector, XVector, bVector, YVector, -1, 0, SizeType> op (r, x, y, av, bv);
    const SizeType numRows = x.dimension_0 ();
    Kokkos::RangePolicy<execution_space, SizeType> policy (0, numRows);
    Kokkos::parallel_for (policy, op);
    return;
  }
  if (a == -1 && b == -1) {
    MV_Axpby_Functor<RVector, aVector, XVector, bVector, YVector, -1, -1, SizeType> op (r, x, y, av, bv);
    const SizeType numRows = x.dimension_0 ();
    Kokkos::RangePolicy<execution_space, SizeType> policy (0, numRows);
    Kokkos::parallel_for (policy, op);
    return;
  }
  if (a == -1 && b == 1) {
    MV_Axpby_Functor<RVector, aVector, XVector, bVector, YVector, -1, 1, SizeType> op (r, x, y, av, bv);
    const SizeType numRows = x.dimension_0 ();
    Kokkos::RangePolicy<execution_space, SizeType> policy (0, numRows);
    Kokkos::parallel_for (policy, op);
    return;
  }
  if (a == -1 && b == 2) {
    MV_Axpby_Functor<RVector, aVector, XVector, bVector, YVector, -1, 2, SizeType> op (r, x, y, av, bv);
    const SizeType numRows = x.dimension_0 ();
    Kokkos::RangePolicy<execution_space, SizeType> policy (0, numRows);
    Kokkos::parallel_for (policy, op);
    return;
  }
  // a == 1
  if (a == 1 && b == 0) {
    MV_Axpby_Functor<RVector, aVector, XVector, bVector, YVector, 1, 0, SizeType> op (r, x, y, av, bv);
    const SizeType numRows = x.dimension_0 ();
    Kokkos::RangePolicy<execution_space, SizeType> policy (0, numRows);
    Kokkos::parallel_for (policy, op);
    return;
  }
  if (a == 1 && b == -1) {
    MV_Axpby_Functor<RVector, aVector, XVector, bVector, YVector, 1, -1, SizeType> op (r, x, y, av, bv);
    const SizeType numRows = x.dimension_0 ();
    Kokkos::RangePolicy<execution_space, SizeType> policy (0, numRows);
    Kokkos::parallel_for (policy, op);
    return;
  }
  if (a == 1 && b == 1) {
    MV_Axpby_Functor<RVector, aVector, XVector, bVector, YVector, 1, 1, SizeType> op (r, x, y, av, bv);
    const SizeType numRows = x.dimension_0 ();
    Kokkos::RangePolicy<execution_space, SizeType> policy (0, numRows);
    Kokkos::parallel_for (policy, op);
    return;
  }
  if (a == 1 && b == 2) {
    MV_Axpby_Functor<RVector, aVector, XVector, bVector, YVector, 1, 2, SizeType> op (r, x, y, av, bv);
    const SizeType numRows = x.dimension_0 ();
    Kokkos::RangePolicy<execution_space, SizeType> policy (0, numRows);
    Kokkos::parallel_for (policy, op);
    return;
  }
  // a and b arbitrary (not -1, 0, or 1)

  MV_Axpby_Functor<RVector, aVector, XVector, bVector, YVector, 2, 2, SizeType> op (r, x, y, av, bv);
  const SizeType numRows = x.dimension_0 ();
  Kokkos::RangePolicy<execution_space, SizeType> policy (0, numRows);
  Kokkos::parallel_for (policy, op);
}

// Variant of MV_Axpby_Generic for single vectors (1-D Views) r, x,
// and y.  As above, either av and bv are both 1-D Views (and only the
// first entry of each will be read), or both av and bv are scalars.
template<class RV, class aVector, class XV,
         class bVector, class YV, class SizeType>
void
V_Axpby_Generic (const RV& r, const aVector& av, const XV& x,
                 const bVector& bv, const YV& y,
                 int a = 2, int b = 2)
{
  typedef typename RV::execution_space execution_space;

  if (a == 0 && b == 0) {
    V_Axpby_Functor<RV, aVector, XV, bVector, YV, 0, 0, SizeType> op (r, x, y, av, bv);
    const SizeType numRows = x.dimension_0 ();
    Kokkos::RangePolicy<execution_space, SizeType> policy (0, numRows);
    Kokkos::parallel_for (policy, op);
    return;
  }
  if (a == 0 && b == -1) {
    V_Axpby_Functor<RV, aVector, XV, bVector, YV, 0, -1, SizeType> op (r, x, y, av, bv);
    const SizeType numRows = x.dimension_0 ();
    Kokkos::RangePolicy<execution_space, SizeType> policy (0, numRows);
    Kokkos::parallel_for (policy, op);
    return;
  }
  if (a == 0 && b == 1) {
    V_Axpby_Functor<RV, aVector, XV, bVector, YV, 0, 1, SizeType> op (r, x, y, av, bv);
    const SizeType numRows = x.dimension_0 ();
    Kokkos::RangePolicy<execution_space, SizeType> policy (0, numRows);
    Kokkos::parallel_for (policy, op);
    return;
  }
  if (a == 0 && b == 2) {
    V_Axpby_Functor<RV, aVector, XV, bVector, YV, 0, 2, SizeType> op (r, x, y, av, bv);
    const SizeType numRows = x.dimension_0 ();
    Kokkos::RangePolicy<execution_space, SizeType> policy (0, numRows);
    Kokkos::parallel_for (policy, op);
    return;
  }
  // a == -1
  if (a == -1 && b == 0) {
    V_Axpby_Functor<RV, aVector, XV, bVector, YV, -1, 0, SizeType> op (r, x, y, av, bv);
    const SizeType numRows = x.dimension_0 ();
    Kokkos::RangePolicy<execution_space, SizeType> policy (0, numRows);
    Kokkos::parallel_for (policy, op);
    return;
  }
  if (a == -1 && b == -1) {
    V_Axpby_Functor<RV, aVector, XV, bVector, YV, -1, -1, SizeType> op (r, x, y, av, bv);
    const SizeType numRows = x.dimension_0 ();
    Kokkos::RangePolicy<execution_space, SizeType> policy (0, numRows);
    Kokkos::parallel_for (policy, op);
    return;
  }
  if (a == -1 && b == 1) {
    V_Axpby_Functor<RV, aVector, XV, bVector, YV, -1, 1, SizeType> op (r, x, y, av, bv);
    const SizeType numRows = x.dimension_0 ();
    Kokkos::RangePolicy<execution_space, SizeType> policy (0, numRows);
    Kokkos::parallel_for (policy, op);
    return;
  }
  if (a == -1 && b == 2) {
    V_Axpby_Functor<RV, aVector, XV, bVector, YV, -1, 2, SizeType> op (r, x, y, av, bv);
    const SizeType numRows = x.dimension_0 ();
    Kokkos::RangePolicy<execution_space, SizeType> policy (0, numRows);
    Kokkos::parallel_for (policy, op);
    return;
  }
  // a == 1
  if (a == 1 && b == 0) {
    V_Axpby_Functor<RV, aVector, XV, bVector, YV, 1, 0, SizeType> op (r, x, y, av, bv);
    const SizeType numRows = x.dimension_0 ();
    Kokkos::RangePolicy<execution_space, SizeType> policy (0, numRows);
    Kokkos::parallel_for (policy, op);
    return;
  }
  if (a == 1 && b == -1) {
    V_Axpby_Functor<RV, aVector, XV, bVector, YV, 1, -1, SizeType> op (r, x, y, av, bv);
    const SizeType numRows = x.dimension_0 ();
    Kokkos::RangePolicy<execution_space, SizeType> policy (0, numRows);
    Kokkos::parallel_for (policy, op);
    return;
  }
  if (a == 1 && b == 1) {
    V_Axpby_Functor<RV, aVector, XV, bVector, YV, 1, 1, SizeType> op (r, x, y, av, bv);
    const SizeType numRows = x.dimension_0 ();
    Kokkos::RangePolicy<execution_space, SizeType> policy (0, numRows);
    Kokkos::parallel_for (policy, op);
    return;
  }
  if (a == 1 && b == 2) {
    V_Axpby_Functor<RV, aVector, XV, bVector, YV, 1, 2, SizeType> op (r, x, y, av, bv);
    const SizeType numRows = x.dimension_0 ();
    Kokkos::RangePolicy<execution_space, SizeType> policy (0, numRows);
    Kokkos::parallel_for (policy, op);
    return;
  }

  // a and b arbitrary (not -1, 0, or 1)
  V_Axpby_Functor<RV, aVector, XV, bVector, YV, 2, 2, SizeType> op (r, x, y, av, bv);
  const SizeType numRows = x.dimension_0 ();
  Kokkos::RangePolicy<execution_space, SizeType> policy (0, numRows);
  Kokkos::parallel_for (policy, op);
}

// Compute any of the following, in a way optimized for X, Y, and R
// being LayoutLeft:
//
// 1. R(i,j) = a*X(i,j) + b*Y(i,j) for a,b in -1,0,1
// 2. R(i,j) = av(j)*X(i,j) + b*Y(i,j) for b in -1,0,1
// 3. R(i,j) = a*X(i,j) + b*Y(i,j) for a in -1,0,1
// 4. R(i,j) = av(j)*X(i,j) + bv(j)*Y(i,j)
//
// a and b come in as integers.  The values -1, 0, and 1 correspond to
// the literal values of the coefficients.  The value 2 tells the
// functor to use the corresponding vector of coefficients: a == 2
// means use av, and b == 2 means use bv.  Otherwise, av resp. vb are
// ignored.
//
// Any literal coefficient of zero has BLAS semantics of ignoring the
// corresponding (multi)vector entry.  This does NOT apply to
// coefficients in av and bv vectors, if they are used.
//
// Either av and bv are both 1-D Views, or av and bv are both scalars.
template<class RMV, class aVector, class XMV,
         class bVector, class YMV, class SizeType>
void
MV_Axpby_Invoke_Left (const RMV& r, const aVector& av, const XMV& x,
                      const bVector& bv, const YMV& y,
                      int a = 2, int b = 2)
{
  const SizeType numCols = x.dimension_1 ();

  switch (numCols) {
  case 1: {
    typedef Kokkos::View<typename RMV::value_type*, typename RMV::array_layout,
      typename RMV::device_type, typename RMV::memory_traits,
      typename RMV::specialize> RV;
    typedef Kokkos::View<typename XMV::value_type*, typename XMV::array_layout,
      typename XMV::device_type, typename XMV::memory_traits,
      typename XMV::specialize> XV;
    typedef Kokkos::View<typename YMV::value_type*, typename YMV::array_layout,
      typename YMV::device_type, typename YMV::memory_traits,
      typename YMV::specialize> YV;

    RV r_0 = Kokkos::subview (r, Kokkos::ALL (), 0);
    XV x_0 = Kokkos::subview (x, Kokkos::ALL (), 0);
    YV y_0 = Kokkos::subview (y, Kokkos::ALL (), 0);
    V_Axpby_Generic<RV, aVector, XV, bVector, YV, SizeType> (r_0, av, x_0, bv, y_0, a, b);
    break;
  }
  case 2:
    MV_Axpby_Unrolled<RMV, aVector, XMV, bVector, YMV, 2, SizeType> (r, av, x, bv, y, a, b);
    break;
  case 3:
    MV_Axpby_Unrolled<RMV, aVector, XMV, bVector, YMV, 3, SizeType> (r, av, x, bv, y, a, b);
    break;
  case 4:
    MV_Axpby_Unrolled<RMV, aVector, XMV, bVector, YMV, 4, SizeType> (r, av, x, bv, y, a, b);
    break;
  case 5:
    MV_Axpby_Unrolled<RMV, aVector, XMV, bVector, YMV, 5, SizeType> (r, av, x, bv, y, a, b);
    break;
  case 6:
    MV_Axpby_Unrolled<RMV, aVector, XMV, bVector, YMV, 6, SizeType> (r, av, x, bv, y, a, b);
    break;
  case 7:
    MV_Axpby_Unrolled<RMV, aVector, XMV, bVector, YMV, 7, SizeType> (r, av, x, bv, y, a, b);
    break;
  case 8:
    MV_Axpby_Unrolled<RMV, aVector, XMV, bVector, YMV, 8, SizeType> (r, av, x, bv, y, a, b);
    break;
  case 9:
    MV_Axpby_Unrolled<RMV, aVector, XMV, bVector, YMV, 9, SizeType> (r, av, x, bv, y, a, b);
    break;
  case 10:
    MV_Axpby_Unrolled<RMV, aVector, XMV, bVector, YMV, 10, SizeType> (r, av, x, bv, y, a, b);
    break;
  case 11:
    MV_Axpby_Unrolled<RMV, aVector, XMV, bVector, YMV, 11, SizeType> (r, av, x, bv, y, a, b);
    break;
  case 12:
    MV_Axpby_Unrolled<RMV, aVector, XMV, bVector, YMV, 12, SizeType> (r, av, x, bv, y, a, b);
    break;
  case 13:
    MV_Axpby_Unrolled<RMV, aVector, XMV, bVector, YMV, 13, SizeType> (r, av, x, bv, y, a, b);
    break;
  case 14:
    MV_Axpby_Unrolled<RMV, aVector, XMV, bVector, YMV, 14, SizeType> (r, av, x, bv, y, a, b);
    break;
  case 15:
    MV_Axpby_Unrolled<RMV, aVector, XMV, bVector, YMV, 15, SizeType> (r, av, x, bv, y, a, b);
    break;
  case 16:
    MV_Axpby_Unrolled<RMV, aVector, XMV, bVector, YMV, 16, SizeType> (r, av, x, bv, y, a, b);
    break;
  default:
    MV_Axpby_Generic<RMV, aVector, XMV, bVector, YMV, SizeType> (r, av, x, bv, y, a, b);
  }
}

// Compute any of the following, in a way optimized for X, Y, and R
// being LayoutRight:
//
// 1. R(i,j) = a*X(i,j) + b*Y(i,j) for a,b in -1,0,1
// 2. R(i,j) = av(j)*X(i,j) + b*Y(i,j) for b in -1,0,1
// 3. R(i,j) = a*X(i,j) + b*Y(i,j) for a in -1,0,1
// 4. R(i,j) = av(j)*X(i,j) + bv(j)*Y(i,j)
//
// a and b come in as integers.  The values -1, 0, and 1 correspond to
// the literal values of the coefficients.  The value 2 tells the
// functor to use the corresponding vector of coefficients: a == 2
// means use av, and b == 2 means use bv.  Otherwise, av resp. vb are
// ignored.
//
// Any literal coefficient of zero has BLAS semantics of ignoring the
// corresponding (multi)vector entry.  This does NOT apply to
// coefficients in av and bv vectors, if they are used.
//
// Either av and bv are both 1-D Views, or av and bv are both scalars.
template<class RMV, class aVector, class XMV,
         class bVector, class YMV, class SizeType>
void
MV_Axpby_Invoke_Right (const RMV& r, const aVector& av, const XMV& x,
                       const bVector& bv, const YMV& y,
                       int a = 2, int b = 2)
{
  const SizeType numCols = x.dimension_1 ();

  if (numCols == 1) {
    typedef Kokkos::View<typename RMV::value_type*, typename RMV::array_layout,
      typename RMV::device_type, typename RMV::memory_traits,
      typename RMV::specialize> RV;
    typedef Kokkos::View<typename XMV::value_type*, typename XMV::array_layout,
      typename XMV::device_type, typename XMV::memory_traits,
      typename XMV::specialize> XV;
    typedef Kokkos::View<typename YMV::value_type*, typename YMV::array_layout,
      typename YMV::device_type, typename YMV::memory_traits,
      typename YMV::specialize> YV;

    RV r_0 = Kokkos::subview (r, Kokkos::ALL (), 0);
    XV x_0 = Kokkos::subview (x, Kokkos::ALL (), 0);
    YV y_0 = Kokkos::subview (y, Kokkos::ALL (), 0);
    V_Axpby_Generic<RMV, aVector, XMV, bVector, YMV, 1, SizeType> (r_0, av, x_0, bv, y_0, a, b);
  }
  else {
    MV_Axpby_Generic<RMV, aVector, XMV, bVector, YMV, SizeType> (r, av, x, bv, y, a, b);
  }
}

//! Implementation of KokkosBlas::axpby for multivectors.
template<class RT, class RL, class RD, class RM, class RS,
         class XT, class XL, class XD, class XM, class XS,
         class YT, class YL, class YD, class YM, class YS>
struct Axpby_MV {
  typedef Kokkos::View<RT,RL,RD,RM,RS> RMV;
  typedef Kokkos::View<XT,XL,XD,XM,XS> XMV;
  typedef Kokkos::View<YT,YL,YD,YM,YS> YMV;
  typedef typename XMV::size_type size_type;
  typedef Kokkos::Details::ArithTraits<typename XMV::non_const_value_type> ATA;
  typedef Kokkos::Details::ArithTraits<typename YMV::non_const_value_type> ATB;

  static void
  axpby (const RMV& R, const typename XMV::non_const_value_type& alpha,
         const XMV& X, const typename YMV::non_const_value_type& beta,
         const YMV& Y)
  {
    const size_type numRows = X.dimension_0 ();
    const size_type numCols = X.dimension_1 ();
    int a, b;
    if (alpha == ATA::zero ()) {
      a = 0;
    }
    else if (alpha == -ATA::one ()) {
      a = -1;
    }
    else if (alpha == ATA::one ()) {
      a = 1;
    }
    else {
      a = 2;
    }
    if (beta == ATB::zero ()) {
      b = 0;
    }
    else if (beta == -ATB::one ()) {
      b = -1;
    }
    else if (beta == ATB::one ()) {
      b = 1;
    }
    else {
      b = 2;
    }

    if (numRows < static_cast<size_type> (INT_MAX) &&
        numRows * numCols < static_cast<size_type> (INT_MAX)) {
      typedef int index_type;
      MV_Axpby_Invoke_Left<RMV, typename XMV::non_const_value_type, XMV,
        typename YMV::non_const_value_type, YMV, index_type> (R, alpha, X,
                                                              beta, Y, a, b);
    }
    else {
      typedef typename XMV::size_type index_type;
      MV_Axpby_Invoke_Left<RMV, typename XMV::non_const_value_type, XMV,
        typename YMV::non_const_value_type, YMV, index_type> (R, alpha, X,
                                                              beta, Y, a, b);
    }
  }
};

// Compute any of the following:
//
// 1. R(i,j) = a*X(i,j) + b*Y(i,j) for a,b in -1,0,1
// 2. R(i,j) = av(j)*X(i,j) + b*Y(i,j) for b in -1,0,1
// 3. R(i,j) = a*X(i,j) + b*Y(i,j) for a in -1,0,1
// 4. R(i,j) = av(j)*X(i,j) + bv(j)*Y(i,j)
//
// a and b come in as integers.  The values -1, 0, and 1 correspond to
// the literal values of the coefficients.  The value 2 tells the
// functor to use the corresponding vector of coefficients: a == 2
// means use av, and b == 2 means use bv.  Otherwise, av resp. vb are
// ignored.
//
// Any literal coefficient of zero has BLAS semantics of ignoring the
// corresponding (multi)vector entry.  This does NOT apply to
// coefficients in av and bv vectors, if they are used.
//
// Either av and bv are both 1-D Views, or av and bv are both scalars.
template<class RT, class RL, class RD, class RM, class RS,
         class AT, class AL, class AD, class AM, class AS,
         class XT, class XL, class XD, class XM, class XS,
         class BT, class BL, class BD, class BM, class BS,
         class YT, class YL, class YD, class YM, class YS>
struct Axpby_MV_V {
  typedef Kokkos::View<RT,RL,RD,RM,RS> RMV;
  typedef Kokkos::View<AT,AL,AD,AM,AS> AV;
  typedef Kokkos::View<XT,XL,XD,XM,XS> XMV;
  typedef Kokkos::View<BT,BL,BD,BM,BS> BV;
  typedef Kokkos::View<YT,YL,YD,YM,YS> YMV;
  typedef typename XMV::size_type size_type;
  typedef Kokkos::Details::ArithTraits<typename XMV::non_const_value_type> ATA;
  typedef Kokkos::Details::ArithTraits<typename YMV::non_const_value_type> ATB;

  static void
  axpby (const RMV& R, const AV& av, const XMV& X, const BV& bv, const YMV& Y)
  {
    const size_type numRows = X.dimension_0 ();
    const size_type numCols = X.dimension_1 ();
    int a = 2, b = 2;
    if (av.dimension_0 () == 0) {
      a = 0;
    }
    if (bv.dimension_0 () == 0) {
      b = 0;
    }
    if (numRows < static_cast<size_type> (INT_MAX) &&
        numRows * numCols < static_cast<size_type> (INT_MAX)) {
      typedef int index_type;
      MV_Axpby_Invoke_Left<RMV, AV, XMV, BV, YMV, index_type> (R, av, X,
                                                               bv, Y, a, b);
    }
    else {
      typedef typename XMV::size_type index_type;
      MV_Axpby_Invoke_Left<RMV, AV, XMV, BV, YMV, index_type> (R, av, X,
                                                               bv, Y, a, b);
    }
  }
};

#ifdef KOKKOS_HAVE_SERIAL

template<>
struct Axpby_MV<double**,
                Kokkos::LayoutLeft,
                Kokkos::Device<Kokkos::Serial, Kokkos::HostSpace>,
                Kokkos::MemoryTraits<Kokkos::Unmanaged>,
                Kokkos::Impl::ViewDefault,
                const double**,
                Kokkos::LayoutLeft,
                Kokkos::Device<Kokkos::Serial, Kokkos::HostSpace>,
                Kokkos::MemoryTraits<Kokkos::Unmanaged>,
                Kokkos::Impl::ViewDefault,
                const double**,
                Kokkos::LayoutLeft,
                Kokkos::Device<Kokkos::Serial, Kokkos::HostSpace>,
                Kokkos::MemoryTraits<Kokkos::Unmanaged>,
                Kokkos::Impl::ViewDefault> {
  typedef double** RT;
  typedef Kokkos::LayoutLeft RL;
  typedef Kokkos::Device<Kokkos::Serial, Kokkos::HostSpace> RD;
  typedef Kokkos::MemoryTraits<Kokkos::Unmanaged> RM;
  typedef Kokkos::Impl::ViewDefault RS;
  typedef Kokkos::View<RT,RL,RD,RM,RS> RMV;

  typedef const double** XT;
  typedef Kokkos::LayoutLeft XL;
  typedef Kokkos::Device<Kokkos::Serial, Kokkos::HostSpace> XD;
  typedef Kokkos::MemoryTraits<Kokkos::Unmanaged> XM;
  typedef Kokkos::Impl::ViewDefault XS;
  typedef Kokkos::View<XT,XL,XD,XM,XS> XMV;

  typedef const double** YT;
  typedef Kokkos::LayoutLeft YL;
  typedef Kokkos::Device<Kokkos::Serial, Kokkos::HostSpace> YD;
  typedef Kokkos::MemoryTraits<Kokkos::Unmanaged> YM;
  typedef Kokkos::Impl::ViewDefault YS;
  typedef Kokkos::View<YT,YL,YD,YM,YS> YMV;

  typedef XMV::size_type size_type;
  typedef Kokkos::Details::ArithTraits<XMV::non_const_value_type> ATA;
  typedef Kokkos::Details::ArithTraits<YMV::non_const_value_type> ATB;

  static void
  axpby (const RMV& R, const XMV::non_const_value_type& alpha, const XMV& X,
         const YMV::non_const_value_type& beta, const YMV& Y);
};

#endif // KOKKOS_HAVE_SERIAL

#ifdef KOKKOS_HAVE_OPENMP

template<>
struct Axpby_MV<double**,
                Kokkos::LayoutLeft,
                Kokkos::Device<Kokkos::OpenMP, Kokkos::HostSpace>,
                Kokkos::MemoryTraits<Kokkos::Unmanaged>,
                Kokkos::Impl::ViewDefault,
                const double**,
                Kokkos::LayoutLeft,
                Kokkos::Device<Kokkos::OpenMP, Kokkos::HostSpace>,
                Kokkos::MemoryTraits<Kokkos::Unmanaged>,
                Kokkos::Impl::ViewDefault,
                const double**,
                Kokkos::LayoutLeft,
                Kokkos::Device<Kokkos::OpenMP, Kokkos::HostSpace>,
                Kokkos::MemoryTraits<Kokkos::Unmanaged>,
                Kokkos::Impl::ViewDefault> {
  typedef double** RT;
  typedef Kokkos::LayoutLeft RL;
  typedef Kokkos::Device<Kokkos::OpenMP, Kokkos::HostSpace> RD;
  typedef Kokkos::MemoryTraits<Kokkos::Unmanaged> RM;
  typedef Kokkos::Impl::ViewDefault RS;
  typedef Kokkos::View<RT,RL,RD,RM,RS> RMV;

  typedef const double** XT;
  typedef Kokkos::LayoutLeft XL;
  typedef Kokkos::Device<Kokkos::OpenMP, Kokkos::HostSpace> XD;
  typedef Kokkos::MemoryTraits<Kokkos::Unmanaged> XM;
  typedef Kokkos::Impl::ViewDefault XS;
  typedef Kokkos::View<XT,XL,XD,XM,XS> XMV;

  typedef const double** YT;
  typedef Kokkos::LayoutLeft YL;
  typedef Kokkos::Device<Kokkos::OpenMP, Kokkos::HostSpace> YD;
  typedef Kokkos::MemoryTraits<Kokkos::Unmanaged> YM;
  typedef Kokkos::Impl::ViewDefault YS;
  typedef Kokkos::View<YT,YL,YD,YM,YS> YMV;

  typedef XMV::size_type size_type;
  typedef Kokkos::Details::ArithTraits<XMV::non_const_value_type> ATA;
  typedef Kokkos::Details::ArithTraits<YMV::non_const_value_type> ATB;

  static void
  axpby (const RMV& R, const XMV::non_const_value_type& alpha, const XMV& X,
         const YMV::non_const_value_type& beta, const YMV& Y);
};

#endif // KOKKOS_HAVE_OPENMP

#ifdef KOKKOS_HAVE_PTHREAD

template<>
struct Axpby_MV<double**,
                Kokkos::LayoutLeft,
                Kokkos::Device<Kokkos::Threads, Kokkos::HostSpace>,
                Kokkos::MemoryTraits<Kokkos::Unmanaged>,
                Kokkos::Impl::ViewDefault,
                const double**,
                Kokkos::LayoutLeft,
                Kokkos::Device<Kokkos::Threads, Kokkos::HostSpace>,
                Kokkos::MemoryTraits<Kokkos::Unmanaged>,
                Kokkos::Impl::ViewDefault,
                const double**,
                Kokkos::LayoutLeft,
                Kokkos::Device<Kokkos::Threads, Kokkos::HostSpace>,
                Kokkos::MemoryTraits<Kokkos::Unmanaged>,
                Kokkos::Impl::ViewDefault> {
  typedef double** RT;
  typedef Kokkos::LayoutLeft RL;
  typedef Kokkos::Device<Kokkos::Threads, Kokkos::HostSpace> RD;
  typedef Kokkos::MemoryTraits<Kokkos::Unmanaged> RM;
  typedef Kokkos::Impl::ViewDefault RS;
  typedef Kokkos::View<RT,RL,RD,RM,RS> RMV;

  typedef const double** XT;
  typedef Kokkos::LayoutLeft XL;
  typedef Kokkos::Device<Kokkos::Threads, Kokkos::HostSpace> XD;
  typedef Kokkos::MemoryTraits<Kokkos::Unmanaged> XM;
  typedef Kokkos::Impl::ViewDefault XS;
  typedef Kokkos::View<XT,XL,XD,XM,XS> XMV;

  typedef const double** YT;
  typedef Kokkos::LayoutLeft YL;
  typedef Kokkos::Device<Kokkos::Threads, Kokkos::HostSpace> YD;
  typedef Kokkos::MemoryTraits<Kokkos::Unmanaged> YM;
  typedef Kokkos::Impl::ViewDefault YS;
  typedef Kokkos::View<YT,YL,YD,YM,YS> YMV;

  typedef XMV::size_type size_type;
  typedef Kokkos::Details::ArithTraits<XMV::non_const_value_type> ATA;
  typedef Kokkos::Details::ArithTraits<YMV::non_const_value_type> ATB;

  static void
  axpby (const RMV& R, const XMV::non_const_value_type& alpha, const XMV& X,
         const YMV::non_const_value_type& beta, const YMV& Y);
};

#endif // KOKKOS_HAVE_PTHREAD

#ifdef KOKKOS_HAVE_CUDA

template<>
struct Axpby_MV<double**,
                Kokkos::LayoutLeft,
                Kokkos::Device<Kokkos::Cuda, Kokkos::CudaSpace>,
                Kokkos::MemoryTraits<Kokkos::Unmanaged>,
                Kokkos::Impl::ViewDefault,
                const double**,
                Kokkos::LayoutLeft,
                Kokkos::Device<Kokkos::Cuda, Kokkos::CudaSpace>,
                Kokkos::MemoryTraits<Kokkos::Unmanaged>,
                Kokkos::Impl::ViewDefault,
                const double**,
                Kokkos::LayoutLeft,
                Kokkos::Device<Kokkos::Cuda, Kokkos::CudaSpace>,
                Kokkos::MemoryTraits<Kokkos::Unmanaged>,
                Kokkos::Impl::ViewDefault> {
  typedef double** RT;
  typedef Kokkos::LayoutLeft RL;
  typedef Kokkos::Device<Kokkos::Cuda, Kokkos::CudaSpace> RD;
  typedef Kokkos::MemoryTraits<Kokkos::Unmanaged> RM;
  typedef Kokkos::Impl::ViewDefault RS;
  typedef Kokkos::View<RT,RL,RD,RM,RS> RMV;

  typedef const double** XT;
  typedef Kokkos::LayoutLeft XL;
  typedef Kokkos::Device<Kokkos::Cuda, Kokkos::CudaSpace> XD;
  typedef Kokkos::MemoryTraits<Kokkos::Unmanaged> XM;
  typedef Kokkos::Impl::ViewDefault XS;
  typedef Kokkos::View<XT,XL,XD,XM,XS> XMV;

  typedef const double** YT;
  typedef Kokkos::LayoutLeft YL;
  typedef Kokkos::Device<Kokkos::Cuda, Kokkos::CudaSpace> YD;
  typedef Kokkos::MemoryTraits<Kokkos::Unmanaged> YM;
  typedef Kokkos::Impl::ViewDefault YS;
  typedef Kokkos::View<YT,YL,YD,YM,YS> YMV;

  typedef XMV::size_type size_type;
  typedef Kokkos::Details::ArithTraits<XMV::non_const_value_type> ATA;
  typedef Kokkos::Details::ArithTraits<YMV::non_const_value_type> ATB;

  static void
  axpby (const RMV& R, const XMV::non_const_value_type& alpha, const XMV& X,
         const YMV::non_const_value_type& beta, const YMV& Y);
};

#endif // KOKKOS_HAVE_CUDA

#ifdef KOKKOS_HAVE_CUDA

template<>
struct Axpby_MV<double**,
                Kokkos::LayoutLeft,
                Kokkos::Device<Kokkos::Cuda, Kokkos::CudaUVMSpace>,
                Kokkos::MemoryTraits<Kokkos::Unmanaged>,
                Kokkos::Impl::ViewDefault,
                const double**,
                Kokkos::LayoutLeft,
                Kokkos::Device<Kokkos::Cuda, Kokkos::CudaUVMSpace>,
                Kokkos::MemoryTraits<Kokkos::Unmanaged>,
                Kokkos::Impl::ViewDefault,
                const double**,
                Kokkos::LayoutLeft,
                Kokkos::Device<Kokkos::Cuda, Kokkos::CudaUVMSpace>,
                Kokkos::MemoryTraits<Kokkos::Unmanaged>,
                Kokkos::Impl::ViewDefault> {
  typedef double** RT;
  typedef Kokkos::LayoutLeft RL;
  typedef Kokkos::Device<Kokkos::Cuda, Kokkos::CudaUVMSpace> RD;
  typedef Kokkos::MemoryTraits<Kokkos::Unmanaged> RM;
  typedef Kokkos::Impl::ViewDefault RS;
  typedef Kokkos::View<RT,RL,RD,RM,RS> RMV;

  typedef const double** XT;
  typedef Kokkos::LayoutLeft XL;
  typedef Kokkos::Device<Kokkos::Cuda, Kokkos::CudaUVMSpace> XD;
  typedef Kokkos::MemoryTraits<Kokkos::Unmanaged> XM;
  typedef Kokkos::Impl::ViewDefault XS;
  typedef Kokkos::View<XT,XL,XD,XM,XS> XMV;

  typedef const double** YT;
  typedef Kokkos::LayoutLeft YL;
  typedef Kokkos::Device<Kokkos::Cuda, Kokkos::CudaUVMSpace> YD;
  typedef Kokkos::MemoryTraits<Kokkos::Unmanaged> YM;
  typedef Kokkos::Impl::ViewDefault YS;
  typedef Kokkos::View<YT,YL,YD,YM,YS> YMV;

  typedef XMV::size_type size_type;
  typedef Kokkos::Details::ArithTraits<XMV::non_const_value_type> ATA;
  typedef Kokkos::Details::ArithTraits<YMV::non_const_value_type> ATB;

  static void
  axpby (const RMV& R, const XMV::non_const_value_type& alpha, const XMV& X,
         const YMV::non_const_value_type& beta, const YMV& Y);
};

#endif // KOKKOS_HAVE_CUDA

//
// scal
//

// Functor for multivectors R and X and 1-D View a, that computes any
// of the following:
//
// 1. R(i,j) = alpha*X(i,j) for alpha in -1,0,1
// 2. R(i,j) = a(j)*X(i,j)
//
// The template parameter scalar_x corresponds to alpha in the
// operation y = alpha*x.  The values -1, 0, and -1 correspond to
// literal values of this coefficient.  The value 2 tells the functor
// to use the corresponding vector of coefficients.  Any literal
// coefficient of zero has BLAS semantics of ignoring the
// corresponding (multi)vector entry.  This does not apply to
// coefficients in the a vector, if they are used.
template<class RMV, class aVector, class XMV, int scalar_x,
         class SizeType = typename RMV::size_type>
struct MV_Scal_Functor
{
  typedef typename RMV::execution_space execution_space;
  typedef SizeType                            size_type;
  typedef Kokkos::Details::ArithTraits<typename RMV::non_const_value_type> ATS;

  const size_type numCols;
  RMV R_;
  XMV X_;
  aVector a_;

  MV_Scal_Functor (const RMV& R, const XMV& X, const aVector& a) :
    numCols (X.dimension_1 ()), R_ (R), X_ (X), a_ (a)
  {}

  KOKKOS_INLINE_FUNCTION
  void operator() (const size_type& i) const
  {
    // scalar_x is a compile-time constant (since it is a template
    // parameter), so the compiler should evaluate these branches at
    // compile time.
    if (scalar_x == 0) {
#ifdef KOKKOS_HAVE_PRAGMA_IVDEP
#pragma ivdep
#endif
#ifdef KOKKOS_HAVE_PRAGMA_VECTOR
#pragma vector always
#endif
      for (size_type k = 0; k < numCols; ++k) {
        R_(i,k) = ATS::zero ();
      }
    }
    if (scalar_x == -1) {
#ifdef KOKKOS_HAVE_PRAGMA_IVDEP
#pragma ivdep
#endif
#ifdef KOKKOS_HAVE_PRAGMA_VECTOR
#pragma vector always
#endif
      for (size_type k = 0; k < numCols; ++k) {
        R_(i,k) = -X_(i,k);
      }
    }
    if (scalar_x == 1) {
#ifdef KOKKOS_HAVE_PRAGMA_IVDEP
#pragma ivdep
#endif
#ifdef KOKKOS_HAVE_PRAGMA_VECTOR
#pragma vector always
#endif
      for (size_type k = 0; k < numCols; ++k) {
        R_(i,k) = X_(i,k);
      }
    }
    if (scalar_x == 2) {
#ifdef KOKKOS_HAVE_PRAGMA_IVDEP
#pragma ivdep
#endif
#ifdef KOKKOS_HAVE_PRAGMA_VECTOR
#pragma vector always
#endif
      for (size_type k = 0; k < numCols; ++k) {
        R_(i,k) = a_(k)*X_(i,k);
      }
    }
  }
};

// Variant of MV_Scal_Functor, where a is a scalar.
// This functor computes any of the following:
//
// 1. R(i,j) = alpha*X(i,j) for alpha,beta in -1,0,1
// 2. R(i,j) = a*X(i,j)
//
// This version works by partial specialization on aVector.
// In this partial specialization, aVector is a scalar.
template<class RMV, class XMV, int scalar_x, class SizeType>
struct MV_Scal_Functor<RMV, typename XMV::non_const_value_type,
                       XMV, scalar_x, SizeType>
{
  typedef typename RMV::execution_space execution_space;
  typedef SizeType                            size_type;
  typedef Kokkos::Details::ArithTraits<typename RMV::non_const_value_type> ATS;

  const size_type numCols;
  RMV m_r;
  XMV m_x;
  const typename XMV::non_const_value_type m_a;

  MV_Scal_Functor (const RMV& r, const XMV& x,
                   const typename XMV::non_const_value_type& a) :
    numCols (x.dimension_1 ()), m_r (r), m_x (x), m_a (a)
  {}

  KOKKOS_INLINE_FUNCTION
  void operator() (const size_type& i) const
  {
    // scalar_x and scalar_y are compile-time constants (since they
    // are template parameters), so the compiler should evaluate these
    // branches at compile time.
    if (scalar_x == 0) {
#ifdef KOKKOS_HAVE_PRAGMA_IVDEP
#pragma ivdep
#endif
#ifdef KOKKOS_HAVE_PRAGMA_VECTOR
#pragma vector always
#endif
      for (size_type k = 0; k < numCols; ++k) {
        m_r(i,k) = ATS::zero ();
      }
    }
    if (scalar_x == -1) {
#ifdef KOKKOS_HAVE_PRAGMA_IVDEP
#pragma ivdep
#endif
#ifdef KOKKOS_HAVE_PRAGMA_VECTOR
#pragma vector always
#endif
      for (size_type k = 0; k < numCols; ++k) {
        m_r(i,k) = -m_x(i,k);
      }
    }
    if (scalar_x == 1) {
#ifdef KOKKOS_HAVE_PRAGMA_IVDEP
#pragma ivdep
#endif
#ifdef KOKKOS_HAVE_PRAGMA_VECTOR
#pragma vector always
#endif
      for (size_type k = 0; k < numCols; ++k) {
        m_r(i,k) = m_x(i,k);
      }
    }
    if (scalar_x == 2) {
#ifdef KOKKOS_HAVE_PRAGMA_IVDEP
#pragma ivdep
#endif
#ifdef KOKKOS_HAVE_PRAGMA_VECTOR
#pragma vector always
#endif
      for (size_type k = 0; k < numCols; ++k) {
        m_r(i,k) = m_a*m_x(i,k);
      }
    }
  }
};


// Column-unrolled variant of MV_Scal_Functor.  The number of columns
// in X and Y, UNROLL, is a compile-time constant.
template<class RMV, class aVector, class XMV,
         int scalar_x, int UNROLL, class SizeType>
struct MV_Scal_Unroll_Functor
{
  typedef typename RMV::execution_space execution_space;
  typedef SizeType                            size_type;
  typedef Kokkos::Details::ArithTraits<typename RMV::non_const_value_type> ATS;

  RMV m_r;
  XMV m_x;
  aVector m_a;

  MV_Scal_Unroll_Functor (const RMV& r, const XMV& x, const aVector& a) :
    m_r (r), m_x (x), m_a (a)
  {}

  KOKKOS_INLINE_FUNCTION
  void operator() (const size_type& i) const
  {
    if (scalar_x == 0) {
#ifdef KOKKOS_HAVE_PRAGMA_UNROLL
#pragma unroll
#endif
      for (int k = 0; k < UNROLL; ++k) {
        m_r(i,k) = ATS::zero ();
      }
    }
    if (scalar_x == -1) {
#ifdef KOKKOS_HAVE_PRAGMA_UNROLL
#pragma unroll
#endif
      for (int k = 0; k < UNROLL; ++k) {
        m_r(i,k) = -m_x(i,k);
      }
    }
    if (scalar_x == 1) {
#ifdef KOKKOS_HAVE_PRAGMA_UNROLL
#pragma unroll
#endif
      for (int k = 0; k < UNROLL; ++k) {
        m_r(i,k) = m_x(i,k);
      }
    }
    if (scalar_x == 2) {
#ifdef KOKKOS_HAVE_PRAGMA_UNROLL
#pragma unroll
#endif
      for (int k = 0; k < UNROLL; ++k) {
        m_r(i,k) = m_a(k)*m_x(i,k);
      }
    }
  }
};


// Variant of MV_Scal_Unroll_Functor for a single coefficient (rather
// than a vector of coefficients) a.  The number of columns in X,
// UNROLL, is a compile-time constant.
template<class RMV, class XMV, int scalar_x, int UNROLL, class SizeType>
struct MV_Scal_Unroll_Functor<RMV, typename XMV::non_const_value_type,
                              XMV, scalar_x, UNROLL, SizeType>
{
  typedef typename RMV::execution_space execution_space;
  typedef SizeType                            size_type;
  typedef Kokkos::Details::ArithTraits<typename RMV::non_const_value_type> ATS;

  RMV m_r;
  XMV m_x;
  const typename XMV::non_const_value_type m_a;

  MV_Scal_Unroll_Functor (const RMV& r, const XMV& x,
                          const typename XMV::non_const_value_type& a) :
    m_r (r), m_x (x), m_a (a)
  {}

  KOKKOS_INLINE_FUNCTION
  void operator() (const size_type& i) const
  {
    if (scalar_x == 0) {
#ifdef KOKKOS_HAVE_PRAGMA_UNROLL
#pragma unroll
#endif
      for (int k = 0; k < UNROLL; ++k) {
        m_r(i,k) = ATS::zero ();
      }
    }
    if (scalar_x == -1) {
#ifdef KOKKOS_HAVE_PRAGMA_UNROLL
#pragma unroll
#endif
      for (int k = 0; k < UNROLL; ++k) {
        m_r(i,k) = -m_x(i,k);
      }
    }
    if (scalar_x == 1) {
#ifdef KOKKOS_HAVE_PRAGMA_UNROLL
#pragma unroll
#endif
      for (int k = 0; k < UNROLL; ++k) {
        m_r(i,k) = m_x(i,k);
      }
    }
    if (scalar_x == 2) {
#ifdef KOKKOS_HAVE_PRAGMA_UNROLL
#pragma unroll
#endif
      for (int k = 0; k < UNROLL; ++k) {
        m_r(i,k) = m_a*m_x(i,k);
      }
    }
  }
};

// Single-vector version of MV_Scal_Functor.  By default, a is still a
// 1-D View.  Below is a partial specialization that lets a be a
// scalar.  This functor computes any of the following:
//
// 1. Y(i) = alpha*X(i) for alpha in -1,0,1
// 2. Y(i) = a(0)*X(i)
//
// The template parameter scalar_x corresponds to alpha in the
// operation y = alpha*x + beta*y.  The values -1, 0, and -1
// correspond to literal values of this coefficient.  The value 2
// tells the functor to use the corresponding vector of coefficients.
// Any literal coefficient of zero has BLAS semantics of ignoring the
// corresponding (multi)vector entry.  This does not apply to
// coefficients in the a vector, if used.
template<class RV, class AV, class XV, int scalar_x, class SizeType>
struct V_Scal_Functor {
  typedef typename RV::execution_space execution_space;
  typedef SizeType                           size_type;
  typedef Kokkos::Details::ArithTraits<typename RV::non_const_value_type> ATS;

  RV m_r;
  XV m_x;
  AV m_a;

  V_Scal_Functor (const RV& r, const XV& x, const AV& a) :
    m_r (r), m_x (x), m_a (a)
  {}

  KOKKOS_INLINE_FUNCTION
  void operator() (const size_type& i) const
  {
    // scalar_x is a compile-time constant (since it is a template
    // parameter), so the compiler should evaluate these branches at
    // compile time.
    if (scalar_x == 0) {
      m_r(i) = ATS::zero ();
    }
    if (scalar_x == -1) {
      m_r(i) = -m_x(i);
    }
    if (scalar_x == 1) {
      m_r(i) = m_x(i);
    }
    if (scalar_x == 2) {
      m_r(i) = m_a(0)*m_x(i);
    }
  }
};


// Partial specialization of V_Scal_Functor that lets a be a scalar
// (rather than a 1-D View, as in the most general version above).
// This functor computes any of the following:
//
// 1. Y(i) = alpha*X(i) for alpha in -1,0,1
// 2. Y(i) = a*X(i)
template<class RV, class XV, int scalar_x, class SizeType>
struct V_Scal_Functor<RV, typename XV::non_const_value_type,
                      XV, scalar_x, SizeType> {
  typedef typename RV::execution_space execution_space;
  typedef SizeType                           size_type;
  typedef Kokkos::Details::ArithTraits<typename RV::non_const_value_type> ATS;

  RV m_r;
  XV m_x;
  const typename XV::non_const_value_type m_a;

  V_Scal_Functor (const RV& r, const XV& x,
                  const typename XV::non_const_value_type& a) :
    m_r (r), m_x (x), m_a (a)
  {}

  KOKKOS_INLINE_FUNCTION
  void operator() (const size_type& i) const
  {
    if (scalar_x == 0) {
      m_r(i) = ATS::zero ();
    }
    if (scalar_x == -1) {
      m_r(i) = -m_x(i);
    }
    if (scalar_x == -1) {
      m_r(i) = -m_x(i);
    }
    if (scalar_x == 2) {
      m_r(i) = m_a*m_x(i);
    }
  }
};

// Invoke the unrolled multivector functor that computes any of the
// following:
//
// 1. R(i,j) = a*X(i,j) for a in -1,0,1
// 2. R(i,j) = av(j)*X(i,j)
//
// a comes in as an int.  The values -1, 0, and 1 correspond to the
// literal values of this coefficient.  The value 2 tells the functor
// to use av, which may be either a 1-D View or a scalar.  Otherwise,
// av is ignored.
//
// Any literal coefficient of zero has BLAS semantics of ignoring the
// corresponding (multi)vector entry.  This does NOT apply to
// coefficient(s) in av, if used.
template<class RMV, class aVector, class XMV, int UNROLL, class SizeType>
void
MV_Scal_Unrolled (const RMV& r, const aVector& av, const XMV& x, int a = 2)
{
  typedef typename XMV::execution_space execution_space;

  if (a == 0) {
    MV_Scal_Unroll_Functor<RMV, aVector, XMV, 0, UNROLL, SizeType> op (r, x, av);
    const SizeType numRows = x.dimension_0 ();
    Kokkos::RangePolicy<execution_space, SizeType> policy (0, numRows);
    Kokkos::parallel_for (policy, op);
    return;
  }
  if (a == -1) {
    MV_Scal_Unroll_Functor<RMV, aVector, XMV, -1, UNROLL, SizeType> op (r, x, av);
    const SizeType numRows = x.dimension_0 ();
    Kokkos::RangePolicy<execution_space, SizeType> policy (0, numRows);
    Kokkos::parallel_for (policy, op);
    return;
  }
  if (a == 1) {
    MV_Scal_Unroll_Functor<RMV, aVector, XMV, 1, UNROLL, SizeType> op (r, x, av);
    const SizeType numRows = x.dimension_0 ();
    Kokkos::RangePolicy<execution_space, SizeType> policy (0, numRows);
    Kokkos::parallel_for (policy, op);
    return;
  }

  // a arbitrary (not -1, 0, or 1)
  MV_Scal_Unroll_Functor<RMV, aVector, XMV, 2, UNROLL, SizeType> op (r, x, av);
  const SizeType numRows = x.dimension_0 ();
  Kokkos::RangePolicy<execution_space, SizeType> policy (0, numRows);
  Kokkos::parallel_for (policy, op);
}


// Invoke the "generic" (not unrolled) multivector functor that
// computes any of the following:
//
// 1. R(i,j) = a*X(i,j) for a in -1,0,1
// 2. R(i,j) = av(j)*X(i,j)
//
// a comes in as an int.  The values -1, 0, and 1 correspond to the
// literal values of this coefficient.  The value 2 tells the functor
// to use av, which may be either a 1-D View or a scalar.  Otherwise,
// av is ignored.
//
// Any literal coefficient of zero has BLAS semantics of ignoring the
// corresponding (multi)vector entry.  This does NOT apply to
// coefficient(s) in av, if used.
template<class RVector, class aVector, class XVector, class SizeType>
void
MV_Scal_Generic (const RVector& r, const aVector& av,
                 const XVector& x, int a = 2)
{
  typedef typename XVector::execution_space execution_space;
  const SizeType numRows = x.dimension_0 ();
  Kokkos::RangePolicy<execution_space, SizeType> policy (0, numRows);

  if (a == 0) {
    MV_Scal_Functor<RVector, aVector, XVector, 0, SizeType> op (r, x, av);
    Kokkos::parallel_for (policy, op);
    return;
  }
  if (a == -1) {
    MV_Scal_Functor<RVector, aVector, XVector, -1, SizeType> op (r, x, av);
    Kokkos::parallel_for (policy, op);
    return;
  }
  if (a == 1) {
    MV_Scal_Functor<RVector, aVector, XVector, 1, SizeType> op (r, x, av);
    Kokkos::parallel_for (policy, op);
    return;
  }

  // a arbitrary (not -1, 0, or 1)
  MV_Scal_Functor<RVector, aVector, XVector, 2, SizeType> op (r, x, av);
  Kokkos::parallel_for (policy, op);
}

// Variant of MV_Scal_Generic for single vectors (1-D Views) r and x.
// As above, av is either a 1-D View (and only its first entry will be
// read), or a scalar.
template<class RV, class AV, class XV, class SizeType>
void
V_Scal_Generic (const RV& r, const AV& av, const XV& x, int a = 2)
{
#ifdef KOKKOS_HAVE_CXX11
  static_assert (Kokkos::Impl::is_view<RV>::value,
                 "V_Scal_Generic: RV is not a Kokkos::View.");
  static_assert (Kokkos::Impl::is_view<XV>::value,
                 "V_Scal_Generic: XV is not a Kokkos::View.");
  static_assert (RV::rank == 1,
                 "V_Scal_Generic: RV is not rank 1.");
  static_assert (XV::rank == 1,
                 "V_Scal_Generic: XV is not rank 1.");
#endif // KOKKOS_HAVE_CXX11

  typedef typename RV::execution_space execution_space;
  const SizeType numRows = x.dimension_0 ();
  Kokkos::RangePolicy<execution_space, SizeType> policy (0, numRows);

  if (a == 0) {
    V_Scal_Functor<RV, AV, XV, 0, SizeType> op (r, x, av);
    Kokkos::parallel_for (policy, op);
    return;
  }
  if (a == -1) {
    V_Scal_Functor<RV, AV, XV, -1, SizeType> op (r, x, av);
    Kokkos::parallel_for (policy, op);
    return;
  }
  if (a == 1) {
    V_Scal_Functor<RV, AV, XV, 1, SizeType> op (r, x, av);
    Kokkos::parallel_for (policy, op);
    return;
  }

  // a arbitrary (not -1, 0, or 1)
  V_Scal_Functor<RV, AV, XV, 2, SizeType> op (r, x, av);
  Kokkos::parallel_for (policy, op);
}

// Compute any of the following, in a way optimized for X, Y, and R
// being LayoutLeft:
//
// 1. R(i,j) = a*X(i,j) for a in -1,0,1
// 2. R(i,j) = av(j)*X(i,j)
//
// a comes in as an int.  The values -1, 0, and 1 correspond to the
// literal values of this coefficient.  The value 2 tells the functor
// to use av, which may be either a 1-D View or a scalar.  Otherwise,
// av is ignored.
//
// Any literal coefficient of zero has BLAS semantics of ignoring the
// corresponding (multi)vector entry.  This does NOT apply to
// coefficient(s) in av, if used.
template<class RMV, class AV, class XMV, class SizeType>
void
MV_Scal_Invoke_Left (const RMV& r, const AV& av, const XMV& x, int a = 2)
{
  const SizeType numCols = x.dimension_1 ();

  switch (numCols) {
  case 1: {
    typedef Kokkos::View<typename RMV::value_type*, typename RMV::array_layout,
      typename RMV::device_type, typename RMV::memory_traits,
      typename RMV::specialize> RV;
    typedef Kokkos::View<typename XMV::value_type*, typename XMV::array_layout,
      typename XMV::device_type, typename XMV::memory_traits,
      typename XMV::specialize> XV;

    RV r_0 = Kokkos::subview (r, Kokkos::ALL (), 0);
    XV x_0 = Kokkos::subview (x, Kokkos::ALL (), 0);
    V_Scal_Generic<RV, AV, XV, SizeType> (r_0, av, x_0, a);
    break;
  }
  case 2:
    MV_Scal_Unrolled<RMV, AV, XMV, 2, SizeType> (r, av, x, a);
    break;
  case 3:
    MV_Scal_Unrolled<RMV, AV, XMV, 3, SizeType> (r, av, x, a);
    break;
  case 4:
    MV_Scal_Unrolled<RMV, AV, XMV, 4, SizeType> (r, av, x, a);
    break;
  case 5:
    MV_Scal_Unrolled<RMV, AV, XMV, 5, SizeType> (r, av, x, a);
    break;
  case 6:
    MV_Scal_Unrolled<RMV, AV, XMV, 6, SizeType> (r, av, x, a);
    break;
  case 7:
    MV_Scal_Unrolled<RMV, AV, XMV, 7, SizeType> (r, av, x, a);
    break;
  case 8:
    MV_Scal_Unrolled<RMV, AV, XMV, 8, SizeType> (r, av, x, a);
    break;
  case 9:
    MV_Scal_Unrolled<RMV, AV, XMV, 9, SizeType> (r, av, x, a);
    break;
  case 10:
    MV_Scal_Unrolled<RMV, AV, XMV, 10, SizeType> (r, av, x, a);
    break;
  case 11:
    MV_Scal_Unrolled<RMV, AV, XMV, 11, SizeType> (r, av, x, a);
    break;
  case 12:
    MV_Scal_Unrolled<RMV, AV, XMV, 12, SizeType> (r, av, x, a);
    break;
  case 13:
    MV_Scal_Unrolled<RMV, AV, XMV, 13, SizeType> (r, av, x, a);
    break;
  case 14:
    MV_Scal_Unrolled<RMV, AV, XMV, 14, SizeType> (r, av, x, a);
    break;
  case 15:
    MV_Scal_Unrolled<RMV, AV, XMV, 15, SizeType> (r, av, x, a);
    break;
  case 16:
    MV_Scal_Unrolled<RMV, AV, XMV, 16, SizeType> (r, av, x, a);
    break;
  default:
    MV_Scal_Generic<RMV, AV, XMV, SizeType> (r, av, x, a);
  }
}

// Compute any of the following, in a way optimized for X, Y, and R
// being LayoutRight:
//
// 1. R(i,j) = a*X(i,j) for a in -1,0,1
// 2. R(i,j) = av(j)*X(i,j)
//
// a comes in as an int.  The values -1, 0, and 1 correspond to the
// literal values of this coefficient.  The value 2 tells the functor
// to use av, which may be either a 1-D View or a scalar.  Otherwise,
// av is ignored.
//
// Any literal coefficient of zero has BLAS semantics of ignoring the
// corresponding (multi)vector entry.  This does NOT apply to
// coefficient(s) in av, if used.
template<class RMV, class aVector, class XMV, class SizeType>
void
MV_Scal_Invoke_Right (const RMV& r, const aVector& av, const XMV& x, int a = 2)
{
  const SizeType numCols = x.dimension_1 ();

  if (numCols == 1) {
    typedef Kokkos::View<typename RMV::value_type*, typename RMV::array_layout,
      typename RMV::device_type, typename RMV::memory_traits,
      typename RMV::specialize> RV;
    typedef Kokkos::View<typename XMV::value_type*, typename XMV::array_layout,
      typename XMV::device_type, typename XMV::memory_traits,
      typename XMV::specialize> XV;

    RV r_0 = Kokkos::subview (r, Kokkos::ALL (), 0);
    XV x_0 = Kokkos::subview (x, Kokkos::ALL (), 0);
    V_Scal_Generic<RMV, aVector, XMV, 1, SizeType> (r_0, av, x_0, a);
  }
  else {
    MV_Scal_Generic<RMV, aVector, XMV, SizeType> (r, av, x, a);
  }
}

/// \brief Implementation of KokkosBlas::scal for multivectors.
///
/// Compute any of the following:
///
/// 1. R(i,j) = a*X(i,j) for a in -1,0,1
/// 2. R(i,j) = av(j)*X(i,j)
template<class RMV, class AV, class XMV,
         int rank = RMV::rank>
struct Scal {};

template<class RMV, class AV, class XMV>
struct Scal<RMV, AV, XMV, 2> {
  typedef typename XMV::size_type size_type;
  typedef Kokkos::Details::ArithTraits<typename XMV::non_const_value_type> ATA;

  static void
  scal (const RMV& R, const AV& av, const XMV& X)
  {
#ifdef KOKKOS_HAVE_CXX11
    static_assert (Kokkos::Impl::is_view<RMV>::value, "KokkosBlas::Impl::"
                   "Scal<2-D>: RMV is not a Kokkos::View.");
    static_assert (Kokkos::Impl::is_view<AV>::value, "KokkosBlas::Impl::"
                   "Scal<2-D>: AV is not a Kokkos::View.");
    static_assert (Kokkos::Impl::is_view<XMV>::value, "KokkosBlas::Impl::"
                   "Scal<2-D>: XMV is not a Kokkos::View.");
    static_assert (RMV::rank == 2, "KokkosBlas::Impl::Scal<2-D>: "
                   "RMV is not rank 2.");
    static_assert (AV::rank == 1, "KokkosBlas::Impl::Scal<2-D>: "
                   "AV is not rank 1.");
    static_assert (XMV::rank == 2, "KokkosBlas::Impl::Scal<2-D>: "
                   "XMV is not rank 2.");
#endif // KOKKOS_HAVE_CXX11
    const size_type numRows = X.dimension_0 ();
    const size_type numCols = X.dimension_1 ();
    const int a = (av.dimension_0 () == 0) ? 0 : 2;
    if (numRows < static_cast<size_type> (INT_MAX) &&
        numRows * numCols < static_cast<size_type> (INT_MAX)) {
      typedef int index_type;
      MV_Scal_Invoke_Left<RMV, AV, XMV, index_type> (R, av, X, a);
    }
    else {
      typedef typename XMV::size_type index_type;
      MV_Scal_Invoke_Left<RMV, AV, XMV, index_type> (R, av, X, a);
    }
  }
};

/// \brief Partial specialization of Scal for scalar AV (instead of 1-D View).
///
/// Compute any of the following:
///
/// 1. R(i,j) = a*X(i,j) for a in -1,0,1
/// 2. R(i,j) = alpha*X(i,j)
template<class RMV, class XMV>
struct Scal<RMV, typename XMV::non_const_value_type, XMV, 2> {
  typedef typename XMV::non_const_value_type AV;
  typedef typename XMV::size_type size_type;
  typedef Kokkos::Details::ArithTraits<typename XMV::non_const_value_type> ATA;

  static void
  scal (const RMV& R, const AV& alpha, const XMV& X)
  {
#ifdef KOKKOS_HAVE_CXX11
    static_assert (Kokkos::Impl::is_view<RMV>::value, "KokkosBlas::Impl::"
                   "Scal<2-D, AV=scalar>: RMV is not a Kokkos::View.");
    static_assert (Kokkos::Impl::is_view<XMV>::value, "KokkosBlas::Impl::"
                   "Scal<2-D, AV=scalar>: XMV is not a Kokkos::View.");
    static_assert (RMV::rank == 2, "KokkosBlas::Impl::Scal<2-D, AV=scalar>: "
                   "RMV is not rank 2.");
    static_assert (XMV::rank == 2, "KokkosBlas::Impl::Scal<2-D, AV=scalar>: "
                   "XMV is not rank 2.");
#endif // KOKKOS_HAVE_CXX11

    const size_type numRows = X.dimension_0 ();
    const size_type numCols = X.dimension_1 ();
    int a = 2;
    if (alpha == ATA::zero ()) {
      a = 0;
    }
    else if (alpha == -ATA::one ()) {
      a = -1;
    }
    else if (alpha == ATA::one ()) {
      a = 1;
    }

    if (numRows < static_cast<size_type> (INT_MAX) &&
        numRows * numCols < static_cast<size_type> (INT_MAX)) {
      typedef int index_type;
      MV_Scal_Invoke_Left<RMV, typename XMV::non_const_value_type, XMV,
        index_type> (R, alpha, X, a);
    }
    else {
      typedef typename XMV::size_type index_type;
      MV_Scal_Invoke_Left<RMV, typename XMV::non_const_value_type, XMV,
        index_type> (R, alpha, X, a);
    }
  }
};


/// \brief Partial specialization of Scal for scalar AV (instead of
///   1-D View) and 1-D RMV and XMV.
///
/// Compute any of the following:
///
/// 1. R(i) = a*X(i) for a in -1,0,1
/// 2. R(i) = alpha*X(i)
template<class RMV, class XMV>
struct Scal<RMV, typename RMV::non_const_value_type, XMV, 1>
{
  typedef typename XMV::non_const_value_type AV;
  typedef typename XMV::size_type size_type;
  typedef Kokkos::Details::ArithTraits<typename XMV::non_const_value_type> ATA;

  static void
  scal (const RMV& R, const AV& alpha, const XMV& X)
  {
#ifdef KOKKOS_HAVE_CXX11
    static_assert (Kokkos::Impl::is_view<RMV>::value, "KokkosBlas::Impl::"
                   "Scal<1-D>: RMV is not a Kokkos::View.");
    static_assert (Kokkos::Impl::is_view<XMV>::value, "KokkosBlas::Impl::"
                   "Scal<1-D>: XMV is not a Kokkos::View.");
    static_assert (RMV::rank == 1, "KokkosBlas::Impl::Scal<1-D>: "
                   "RMV is not rank 1.");
    static_assert (XMV::rank == 1, "KokkosBlas::Impl::Scal<1-D>: "
                   "XMV is not rank 1.");
#endif // KOKKOS_HAVE_CXX11

    const size_type numRows = X.dimension_0 ();
    const size_type numCols = X.dimension_1 ();
    int a = 2;
    if (alpha == ATA::zero ()) {
      a = 0;
    }
    else if (alpha == -ATA::one ()) {
      a = -1;
    }
    else if (alpha == ATA::one ()) {
      a = 1;
    }

    if (numRows < static_cast<size_type> (INT_MAX) &&
        numRows * numCols < static_cast<size_type> (INT_MAX)) {
      typedef int index_type;
      V_Scal_Generic<RMV, typename XMV::non_const_value_type, XMV,
        index_type> (R, alpha, X, a);
    }
    else {
      typedef typename XMV::size_type index_type;
      V_Scal_Generic<RMV, typename XMV::non_const_value_type, XMV,
        index_type> (R, alpha, X, a);
    }
  }
};


//
// mfh 08 Apr 2015: For now, we only provide full specializations for
// the AV=scalar case of Scal<RMV, AV, XMV>.  The AV = 1-D View case
// is less commonly used, and the generic kernel should work fine
// there.
//

#ifdef KOKKOS_HAVE_SERIAL

template<>
struct Scal<Kokkos::View<double**,
                         Kokkos::LayoutLeft,
                         Kokkos::Device<Kokkos::Serial, Kokkos::HostSpace>,
                         Kokkos::MemoryTraits<Kokkos::Unmanaged>,
                         Kokkos::Impl::ViewDefault>,
            double,
            Kokkos::View<const double**,
                         Kokkos::LayoutLeft,
                         Kokkos::Device<Kokkos::Serial, Kokkos::HostSpace>,
                         Kokkos::MemoryTraits<Kokkos::Unmanaged>,
                         Kokkos::Impl::ViewDefault>,
            2>
{
  typedef Kokkos::View<double**,
                       Kokkos::LayoutLeft,
                       Kokkos::Device<Kokkos::Serial, Kokkos::HostSpace>,
                       Kokkos::MemoryTraits<Kokkos::Unmanaged>,
                       Kokkos::Impl::ViewDefault> RMV;
  typedef Kokkos::View<double**,
                       Kokkos::LayoutLeft,
                       Kokkos::Device<Kokkos::Serial, Kokkos::HostSpace>,
                       Kokkos::MemoryTraits<Kokkos::Unmanaged>,
                       Kokkos::Impl::ViewDefault> XMV;
  typedef XMV::non_const_value_type AV;
  typedef XMV::size_type size_type;
  typedef Kokkos::Details::ArithTraits<XMV::non_const_value_type> ATA;

  static void
  scal (const RMV& R, XMV::non_const_value_type& alpha, const XMV& X);
};

#endif // KOKKOS_HAVE_SERIAL

#ifdef KOKKOS_HAVE_OPENMP

template<>
struct Scal<Kokkos::View<double**,
                         Kokkos::LayoutLeft,
                         Kokkos::Device<Kokkos::OpenMP, Kokkos::HostSpace>,
                         Kokkos::MemoryTraits<Kokkos::Unmanaged>,
                         Kokkos::Impl::ViewDefault>,
            double,
            Kokkos::View<const double**,
                         Kokkos::LayoutLeft,
                         Kokkos::Device<Kokkos::OpenMP, Kokkos::HostSpace>,
                         Kokkos::MemoryTraits<Kokkos::Unmanaged>,
                         Kokkos::Impl::ViewDefault>, 2>
{
  typedef Kokkos::View<double**,
                       Kokkos::LayoutLeft,
                       Kokkos::Device<Kokkos::OpenMP, Kokkos::HostSpace>,
                       Kokkos::MemoryTraits<Kokkos::Unmanaged>,
                       Kokkos::Impl::ViewDefault> RMV;
  typedef Kokkos::View<double**,
                       Kokkos::LayoutLeft,
                       Kokkos::Device<Kokkos::OpenMP, Kokkos::HostSpace>,
                       Kokkos::MemoryTraits<Kokkos::Unmanaged>,
                       Kokkos::Impl::ViewDefault> XMV;
  typedef XMV::non_const_value_type AV;
  typedef XMV::size_type size_type;
  typedef Kokkos::Details::ArithTraits<XMV::non_const_value_type> ATA;

  static void
  scal (const RMV& R, XMV::non_const_value_type& alpha, const XMV& X);
};

#endif // KOKKOS_HAVE_OPENMP

#ifdef KOKKOS_HAVE_PTHREAD

template<>
struct Scal<Kokkos::View<double**,
                         Kokkos::LayoutLeft,
                         Kokkos::Device<Kokkos::Threads, Kokkos::HostSpace>,
                         Kokkos::MemoryTraits<Kokkos::Unmanaged>,
                         Kokkos::Impl::ViewDefault>,
            double,
            Kokkos::View<const double**,
                         Kokkos::LayoutLeft,
                         Kokkos::Device<Kokkos::Threads, Kokkos::HostSpace>,
                         Kokkos::MemoryTraits<Kokkos::Unmanaged>,
                         Kokkos::Impl::ViewDefault>, 2>
{
  typedef Kokkos::View<double**,
                       Kokkos::LayoutLeft,
                       Kokkos::Device<Kokkos::Threads, Kokkos::HostSpace>,
                       Kokkos::MemoryTraits<Kokkos::Unmanaged>,
                       Kokkos::Impl::ViewDefault> RMV;
  typedef Kokkos::View<double**,
                       Kokkos::LayoutLeft,
                       Kokkos::Device<Kokkos::Threads, Kokkos::HostSpace>,
                       Kokkos::MemoryTraits<Kokkos::Unmanaged>,
                       Kokkos::Impl::ViewDefault> XMV;
  typedef XMV::non_const_value_type AV;
  typedef XMV::size_type size_type;
  typedef Kokkos::Details::ArithTraits<XMV::non_const_value_type> ATA;

  static void
  scal (const RMV& R, XMV::non_const_value_type& alpha, const XMV& X);
};

#endif // KOKKOS_HAVE_PTHREAD

#ifdef KOKKOS_HAVE_CUDA

template<>
struct Scal<Kokkos::View<double**,
                         Kokkos::LayoutLeft,
                         Kokkos::Device<Kokkos::Cuda, Kokkos::CudaSpace>,
                         Kokkos::MemoryTraits<Kokkos::Unmanaged>,
                         Kokkos::Impl::ViewDefault>,
            double,
            Kokkos::View<const double**,
                         Kokkos::LayoutLeft,
                         Kokkos::Device<Kokkos::Cuda, Kokkos::CudaSpace>,
                         Kokkos::MemoryTraits<Kokkos::Unmanaged>,
                         Kokkos::Impl::ViewDefault>, 2>
{
  typedef Kokkos::View<double**,
                       Kokkos::LayoutLeft,
                       Kokkos::Device<Kokkos::Cuda, Kokkos::CudaSpace>,
                       Kokkos::MemoryTraits<Kokkos::Unmanaged>,
                       Kokkos::Impl::ViewDefault> RMV;
  typedef Kokkos::View<double**,
                       Kokkos::LayoutLeft,
                       Kokkos::Device<Kokkos::Cuda, Kokkos::CudaSpace>,
                       Kokkos::MemoryTraits<Kokkos::Unmanaged>,
                       Kokkos::Impl::ViewDefault> XMV;
  typedef XMV::non_const_value_type AV;
  typedef XMV::size_type size_type;
  typedef Kokkos::Details::ArithTraits<XMV::non_const_value_type> ATA;

  static void
  scal (const RMV& R, XMV::non_const_value_type& alpha, const XMV& X);
};

#endif // KOKKOS_HAVE_CUDA

#ifdef KOKKOS_HAVE_CUDA

template<>
struct Scal<Kokkos::View<double**,
                         Kokkos::LayoutLeft,
                         Kokkos::Device<Kokkos::Cuda, Kokkos::CudaUVMSpace>,
                         Kokkos::MemoryTraits<Kokkos::Unmanaged>,
                         Kokkos::Impl::ViewDefault>,
            double,
            Kokkos::View<const double**,
                         Kokkos::LayoutLeft,
                         Kokkos::Device<Kokkos::Cuda, Kokkos::CudaUVMSpace>,
                         Kokkos::MemoryTraits<Kokkos::Unmanaged>,
                         Kokkos::Impl::ViewDefault>, 2>
{
  typedef Kokkos::View<double**,
                       Kokkos::LayoutLeft,
                       Kokkos::Device<Kokkos::Cuda, Kokkos::CudaUVMSpace>,
                       Kokkos::MemoryTraits<Kokkos::Unmanaged>,
                       Kokkos::Impl::ViewDefault> RMV;
  typedef Kokkos::View<double**,
                       Kokkos::LayoutLeft,
                       Kokkos::Device<Kokkos::Cuda, Kokkos::CudaUVMSpace>,
                       Kokkos::MemoryTraits<Kokkos::Unmanaged>,
                       Kokkos::Impl::ViewDefault> XMV;
  typedef XMV::non_const_value_type AV;
  typedef XMV::size_type size_type;
  typedef Kokkos::Details::ArithTraits<XMV::non_const_value_type> ATA;

  static void
  scal (const RMV& R, XMV::non_const_value_type& alpha, const XMV& X);
};

#endif // KOKKOS_HAVE_CUDA

} // namespace Impl
} // namespace KokkosBlas

#endif // KOKKOS_BLAS1_MV_IMPL_HPP_
