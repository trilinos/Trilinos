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
#ifndef KOKKOS_BLAS1_MV_IMPL_DOT_HPP_
#define KOKKOS_BLAS1_MV_IMPL_DOT_HPP_

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
  {
#ifdef KOKKOS_HAVE_CXX11
    static_assert (Kokkos::Impl::is_view<RV>::value, "KokkosBlas::Impl::V_Dot_Functor: "
                   "R is not a Kokkos::View.");
    static_assert (Kokkos::Impl::is_view<XV>::value, "KokkosBlas::Impl::V_Dot_Functor: "
                   "X is not a Kokkos::View.");
    static_assert (Kokkos::Impl::is_view<YV>::value, "KokkosBlas::Impl::V_Dot_Functor: "
                   "Y is not a Kokkos::View.");
    static_assert (Kokkos::Impl::is_same<typename RV::value_type,
                   typename RV::non_const_value_type>::value,
                   "KokkosBlas::Impl::V_Dot_Functor: R is const.  "
                   "It must be nonconst, because it is an output argument "
                   "(we have to be able to write to its entries).");
    static_assert (XV::rank == YV::rank,
                   "KokkosBlas::Impl::V_Dot_Functor: "
                   "X and Y must have the same rank.");
    static_assert (RV::rank == 0 && XV::rank == 1,
                   "KokkosBlas::Impl::V_Dot_Functor: "
                   "RV must have rank 0 and XV and YV must have rank 1.");
#endif // KOKKOS_HAVE_CXX11
  }

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
  {
#ifdef KOKKOS_HAVE_CXX11
    static_assert (Kokkos::Impl::is_view<RV>::value, "KokkosBlas::Impl::MV_Dot_Right_FunctorVector: "
                   "R is not a Kokkos::View.");
    static_assert (Kokkos::Impl::is_view<XMV>::value, "KokkosBlas::Impl::MV_Dot_Right_FunctorVector: "
                   "X is not a Kokkos::View.");
    static_assert (Kokkos::Impl::is_view<YMV>::value, "KokkosBlas::Impl::MV_Dot_Right_FunctorVector: "
                   "Y is not a Kokkos::View.");
    static_assert (Kokkos::Impl::is_same<typename RV::value_type,
                   typename RV::non_const_value_type>::value,
                   "KokkosBlas::Impl::MV_Dot_Right_FunctorVector: R is const.  "
                   "It must be nonconst, because it is an output argument "
                   "(we have to be able to write to its entries).");
    static_assert (XMV::rank == YMV::rank,
                   "KokkosBlas::Impl::MV_Dot_Right_FunctorVector: "
                   "X and Y must have the same rank.");
    static_assert (RV::rank == 1 && XMV::rank == 2,
                   "KokkosBlas::Impl::MV_Dot_Right_FunctorVector: "
                   "RV must have rank 1 and XMV and YMV must have rank 2.");
#endif // KOKKOS_HAVE_CXX11
  }

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
  {
#ifdef KOKKOS_HAVE_CXX11
    static_assert (Kokkos::Impl::is_view<RV>::value, "KokkosBlas::Impl::MV_Dot_Right_FunctorUnroll: "
                   "R is not a Kokkos::View.");
    static_assert (Kokkos::Impl::is_view<XMV>::value, "KokkosBlas::Impl::MV_Dot_Right_FunctorUnroll: "
                   "X is not a Kokkos::View.");
    static_assert (Kokkos::Impl::is_view<YMV>::value, "KokkosBlas::Impl::MV_Dot_Right_FunctorUnroll: "
                   "Y is not a Kokkos::View.");
    static_assert (Kokkos::Impl::is_same<typename RV::value_type,
                   typename RV::non_const_value_type>::value,
                   "KokkosBlas::Impl::MV_Dot_Right_FunctorUnroll: R is const.  "
                   "It must be nonconst, because it is an output argument "
                   "(we have to be able to write to its entries).");
    static_assert (XMV::rank == YMV::rank,
                   "KokkosBlas::Impl::MV_Dot_Right_FunctorUnroll: "
                   "X and Y must have the same rank.");
    static_assert (RV::rank == 1 && XMV::rank == 2,
                   "KokkosBlas::Impl::MV_Dot_Right_FunctorUnroll: "
                   "RV must have rank 1 and XMV and YMV must have rank 2.");
#endif // KOKKOS_HAVE_CXX11
  }

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


template<class RV, class XMV, class YMV, class SizeType>
void
MV_Dot_Invoke (const RV& r, const XMV& X, const YMV& Y)
{
  const SizeType numRows = static_cast<SizeType> (X.dimension_0 ());
  const SizeType numVecs = static_cast<SizeType> (X.dimension_1 ());
  Kokkos::RangePolicy<typename XMV::execution_space, SizeType> policy (0, numRows);

  if (numVecs > 16) {
    MV_Dot_Right_FunctorVector<RV, XMV, YMV, SizeType> op (r, X, Y);
    Kokkos::parallel_reduce (policy, op);
  }
  else {
    switch (numVecs) {
    case 16: {
      MV_Dot_Right_FunctorUnroll<RV, XMV, YMV, 16, SizeType> op (r, X, Y);
      Kokkos::parallel_reduce (policy, op);
      break;
    }
    case 15: {
      MV_Dot_Right_FunctorUnroll<RV, XMV, YMV, 15, SizeType> op (r, X, Y);
      Kokkos::parallel_reduce (policy, op);
      break;
    }
    case 14: {
      MV_Dot_Right_FunctorUnroll<RV, XMV, YMV, 14, SizeType> op (r, X, Y);
      Kokkos::parallel_reduce (policy, op);
      break;
    }
    case 13: {
      MV_Dot_Right_FunctorUnroll<RV, XMV, YMV, 13, SizeType> op (r, X, Y);
      Kokkos::parallel_reduce (policy, op);
      break;
    }
    case 12: {
      MV_Dot_Right_FunctorUnroll<RV, XMV, YMV, 12, SizeType> op (r, X, Y);
      Kokkos::parallel_reduce (policy, op);
      break;
    }
    case 11: {
      MV_Dot_Right_FunctorUnroll<RV, XMV, YMV, 11, SizeType> op (r, X, Y);
      Kokkos::parallel_reduce (policy, op);
      break;
    }
    case 10: {
      MV_Dot_Right_FunctorUnroll<RV, XMV, YMV, 10, SizeType> op (r, X, Y);
      Kokkos::parallel_reduce (policy, op);
      break;
    }
    case 9: {
      MV_Dot_Right_FunctorUnroll<RV, XMV, YMV, 9, SizeType> op (r, X, Y);
      Kokkos::parallel_reduce (policy, op);
      break;
    }
    case 8: {
      MV_Dot_Right_FunctorUnroll<RV, XMV, YMV, 8, SizeType> op (r, X, Y);
      Kokkos::parallel_reduce (policy, op);
      break;
    }
    case 7: {
      MV_Dot_Right_FunctorUnroll<RV, XMV, YMV, 7, SizeType> op (r, X, Y);
      Kokkos::parallel_reduce (policy, op);
      break;
    }
    case 6: {
      MV_Dot_Right_FunctorUnroll<RV, XMV, YMV, 6, SizeType> op (r, X, Y);
      Kokkos::parallel_reduce (policy, op);
      break;
    }
    case 5: {
      MV_Dot_Right_FunctorUnroll<RV, XMV, YMV, 5, SizeType> op (r, X, Y);
      Kokkos::parallel_reduce (policy, op);
      break;
    }
    case 4: {
      MV_Dot_Right_FunctorUnroll<RV, XMV, YMV, 4, SizeType> op (r, X, Y);
      Kokkos::parallel_reduce (policy, op);
      break;
    }
    case 3: {
      MV_Dot_Right_FunctorUnroll<RV, XMV, YMV, 3, SizeType> op (r, X, Y);
      Kokkos::parallel_reduce (policy, op);
      break;
    }
    case 2: {
      MV_Dot_Right_FunctorUnroll<RV, XMV, YMV, 2, SizeType> op (r, X, Y);
      Kokkos::parallel_reduce (policy, op);
      break;
    }
    case 1: {
      // RV needs to turn 0-D, and XMV and YMV need to turn 1-D.
      using Kokkos::ALL;
      using Kokkos::subview;
#ifdef KOKKOS_HAVE_CXX11
      auto r_0 = subview (r, 0);
      auto X_0 = subview (X, ALL (), 0);
      auto Y_0 = subview (Y, ALL (), 0);
      typedef decltype (r_0) RV0D;
      typedef decltype (X_0) XMV1D;
      typedef decltype (Y_0) YMV1D;
#else
      typedef Kokkos::View<typename RV::value_type,
          typename RV::array_layout,
          typename RV::device_type, typename RV::memory_traits,
          typename RV::specialize> RV0D;
      typedef Kokkos::View<typename XMV::const_value_type*,
          typename XMV::device_type, typename XMV::memory_traits,
          typename XMV::specialize> XMV1D;
      typedef Kokkos::View<typename YMV::const_value_type*,
          typename YMV::device_type, typename YMV::memory_traits,
          typename YMV::specialize> YMV1D;
      RV0D r_0 = subview (r, 0);
      XMV1D X_0 = subview (X, ALL (), 0);
      YMV1D Y_0 = subview (Y, ALL (), 0);
#endif // KOKKOS_HAVE_CXX11

      typedef V_Dot_Functor<RV0D, XMV1D, YMV1D, SizeType> op_type;
      op_type op (r_0, X_0, Y_0);
      Kokkos::parallel_reduce (policy, op);
      break;
    }
    } // switch
  } // if-else
}

/// \brief Implementation of KokkosBlas::dot for multivectors or
///   single vectors.
///
/// The fourth template parameter \c rank is the rank of XMV (and
/// YMV).  If 2, they are multivectors; if 1, they are single vectors.
template<class RV, class XMV, class YMV, int rank = XMV::rank>
struct Dot_MV {};

//! Partial specialization for rank = 2 (MultiVectors).
template<class RV, class XMV, class YMV>
struct Dot_MV<RV, XMV, YMV, 2> {
  /// \brief Compute the dot product(s) of the column(s) of the
  ///   multivectors (2-D views) x and y, and store result(s) in the
  ///   1-D View r.
  static void dot (const RV& r, const XMV& X, const YMV& Y)
  {
    typedef typename XMV::size_type size_type;

    const size_type numRows = X.dimension_0 ();
    const size_type numCols = X.dimension_1 ();
    if (numRows < static_cast<size_type> (INT_MAX) &&
        numRows * numCols < static_cast<size_type> (INT_MAX)) {
      MV_Dot_Invoke<RV, XMV, YMV, int> (r, X, Y);
    }
    else {
      MV_Dot_Invoke<RV, XMV, YMV, size_type> (r, X, Y);
    }
  }
};

//! Partial specialization for rank = 1 (single vectors).
template<class RV, class XV, class YV>
struct Dot_MV<RV, XV, YV, 1> {
  /// \brief Compute the dot product of the single vectors X and Y,
  ///   and store result in the 0-D View r.
  static void
  dot (const RV& r, const XV& X, const YV& Y)
  {
    typedef typename XV::size_type size_type;

    const size_type numRows = X.dimension_0 ();
    if (numRows < static_cast<size_type> (INT_MAX)) {
      typedef V_Dot_Functor<RV, XV, YV, int> op_type;
      op_type op (r, X, Y);
      Kokkos::parallel_reduce (numRows, op);
    }
    else {
      typedef V_Dot_Functor<RV, XV, YV, size_type> op_type;
      op_type op (r, X, Y);
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
#define KOKKOSBLAS_IMPL_MV_EXEC_SPACE Kokkos::Serial
#define KOKKOSBLAS_IMPL_MV_MEM_SPACE Kokkos::HostSpace
#define KOKKOSBLAS_IMPL_MV_SCALAR double

template<>
struct Dot_MV<Kokkos::View<Kokkos::Details::InnerProductSpaceTraits<KOKKOSBLAS_IMPL_MV_SCALAR>::dot_type*,
                           KOKKOSBLAS_IMPL_MV_EXEC_SPACE::array_layout,
                           Kokkos::Device<KOKKOSBLAS_IMPL_MV_EXEC_SPACE, KOKKOSBLAS_IMPL_MV_MEM_SPACE>,
                           Kokkos::MemoryTraits<Kokkos::Unmanaged>,
                           Kokkos::Impl::ViewDefault>,
              Kokkos::View<const KOKKOSBLAS_IMPL_MV_SCALAR**,
                           Kokkos::LayoutLeft,
                           Kokkos::Device<KOKKOSBLAS_IMPL_MV_EXEC_SPACE, KOKKOSBLAS_IMPL_MV_MEM_SPACE>,
                           Kokkos::MemoryTraits<Kokkos::Unmanaged>,
                           Kokkos::Impl::ViewDefault>,
              Kokkos::View<const KOKKOSBLAS_IMPL_MV_SCALAR**,
                           Kokkos::LayoutLeft,
                           Kokkos::Device<KOKKOSBLAS_IMPL_MV_EXEC_SPACE, KOKKOSBLAS_IMPL_MV_MEM_SPACE>,
                           Kokkos::MemoryTraits<Kokkos::Unmanaged>,
                           Kokkos::Impl::ViewDefault>,
              2>
{
  typedef Kokkos::View<Kokkos::Details::InnerProductSpaceTraits<KOKKOSBLAS_IMPL_MV_SCALAR>::dot_type*,
    KOKKOSBLAS_IMPL_MV_EXEC_SPACE::array_layout,
    Kokkos::Device<KOKKOSBLAS_IMPL_MV_EXEC_SPACE, KOKKOSBLAS_IMPL_MV_MEM_SPACE>,
    Kokkos::MemoryTraits<Kokkos::Unmanaged>,
    Kokkos::Impl::ViewDefault> RV;
  typedef Kokkos::View<const KOKKOSBLAS_IMPL_MV_SCALAR**,
    Kokkos::LayoutLeft,
    Kokkos::Device<KOKKOSBLAS_IMPL_MV_EXEC_SPACE, KOKKOSBLAS_IMPL_MV_MEM_SPACE>,
    Kokkos::MemoryTraits<Kokkos::Unmanaged>,
    Kokkos::Impl::ViewDefault> XMV;
  typedef Kokkos::View<const KOKKOSBLAS_IMPL_MV_SCALAR**,
    Kokkos::LayoutLeft,
    Kokkos::Device<KOKKOSBLAS_IMPL_MV_EXEC_SPACE, KOKKOSBLAS_IMPL_MV_MEM_SPACE>,
    Kokkos::MemoryTraits<Kokkos::Unmanaged>,
    Kokkos::Impl::ViewDefault> YMV;

  static void dot (const RV& r, const XMV& X, const YMV& Y);
};

template<>
struct Dot_MV<Kokkos::View<Kokkos::Details::InnerProductSpaceTraits<KOKKOSBLAS_IMPL_MV_SCALAR>::dot_type,
                           KOKKOSBLAS_IMPL_MV_EXEC_SPACE::array_layout,
                           Kokkos::Device<KOKKOSBLAS_IMPL_MV_EXEC_SPACE, KOKKOSBLAS_IMPL_MV_MEM_SPACE>,
                           Kokkos::MemoryTraits<Kokkos::Unmanaged>,
                           Kokkos::Impl::ViewDefault>,
              Kokkos::View<const KOKKOSBLAS_IMPL_MV_SCALAR*,
                           Kokkos::LayoutLeft,
                           Kokkos::Device<KOKKOSBLAS_IMPL_MV_EXEC_SPACE, KOKKOSBLAS_IMPL_MV_MEM_SPACE>,
                           Kokkos::MemoryTraits<Kokkos::Unmanaged>,
                           Kokkos::Impl::ViewDefault>,
              Kokkos::View<const KOKKOSBLAS_IMPL_MV_SCALAR*,
                           Kokkos::LayoutLeft,
                           Kokkos::Device<KOKKOSBLAS_IMPL_MV_EXEC_SPACE, KOKKOSBLAS_IMPL_MV_MEM_SPACE>,
                           Kokkos::MemoryTraits<Kokkos::Unmanaged>,
                           Kokkos::Impl::ViewDefault>,
              1>
{
  typedef Kokkos::View<Kokkos::Details::InnerProductSpaceTraits<KOKKOSBLAS_IMPL_MV_SCALAR>::dot_type,
    KOKKOSBLAS_IMPL_MV_EXEC_SPACE::array_layout,
    Kokkos::Device<KOKKOSBLAS_IMPL_MV_EXEC_SPACE, KOKKOSBLAS_IMPL_MV_MEM_SPACE>,
    Kokkos::MemoryTraits<Kokkos::Unmanaged>,
    Kokkos::Impl::ViewDefault> RV;
  typedef Kokkos::View<const KOKKOSBLAS_IMPL_MV_SCALAR*,
    Kokkos::LayoutLeft,
    Kokkos::Device<KOKKOSBLAS_IMPL_MV_EXEC_SPACE, KOKKOSBLAS_IMPL_MV_MEM_SPACE>,
    Kokkos::MemoryTraits<Kokkos::Unmanaged>,
    Kokkos::Impl::ViewDefault> XV;
  typedef Kokkos::View<const KOKKOSBLAS_IMPL_MV_SCALAR*,
    Kokkos::LayoutLeft,
    Kokkos::Device<KOKKOSBLAS_IMPL_MV_EXEC_SPACE, KOKKOSBLAS_IMPL_MV_MEM_SPACE>,
    Kokkos::MemoryTraits<Kokkos::Unmanaged>,
    Kokkos::Impl::ViewDefault> YV;

  static void dot (const RV& r, const XV& X, const YV& Y);
};

#undef KOKKOSBLAS_IMPL_MV_EXEC_SPACE
#undef KOKKOSBLAS_IMPL_MV_MEM_SPACE
#undef KOKKOSBLAS_IMPL_MV_SCALAR
#endif // KOKKOS_HAVE_SERIAL


#ifdef KOKKOS_HAVE_OPENMP
#define KOKKOSBLAS_IMPL_MV_EXEC_SPACE Kokkos::OpenMP
#define KOKKOSBLAS_IMPL_MV_MEM_SPACE Kokkos::HostSpace
#define KOKKOSBLAS_IMPL_MV_SCALAR double
template<>
struct Dot_MV<Kokkos::View<Kokkos::Details::InnerProductSpaceTraits<KOKKOSBLAS_IMPL_MV_SCALAR>::dot_type*,
                           KOKKOSBLAS_IMPL_MV_EXEC_SPACE::array_layout,
                           Kokkos::Device<KOKKOSBLAS_IMPL_MV_EXEC_SPACE, KOKKOSBLAS_IMPL_MV_MEM_SPACE>,
                           Kokkos::MemoryTraits<Kokkos::Unmanaged>,
                           Kokkos::Impl::ViewDefault>,
              Kokkos::View<const KOKKOSBLAS_IMPL_MV_SCALAR**,
                           Kokkos::LayoutLeft,
                           Kokkos::Device<KOKKOSBLAS_IMPL_MV_EXEC_SPACE, KOKKOSBLAS_IMPL_MV_MEM_SPACE>,
                           Kokkos::MemoryTraits<Kokkos::Unmanaged>,
                           Kokkos::Impl::ViewDefault>,
              Kokkos::View<const KOKKOSBLAS_IMPL_MV_SCALAR**,
                           Kokkos::LayoutLeft,
                           Kokkos::Device<KOKKOSBLAS_IMPL_MV_EXEC_SPACE, KOKKOSBLAS_IMPL_MV_MEM_SPACE>,
                           Kokkos::MemoryTraits<Kokkos::Unmanaged>,
                           Kokkos::Impl::ViewDefault>,
              2>
{
  typedef Kokkos::View<Kokkos::Details::InnerProductSpaceTraits<KOKKOSBLAS_IMPL_MV_SCALAR>::dot_type*,
    KOKKOSBLAS_IMPL_MV_EXEC_SPACE::array_layout,
    Kokkos::Device<KOKKOSBLAS_IMPL_MV_EXEC_SPACE, KOKKOSBLAS_IMPL_MV_MEM_SPACE>,
    Kokkos::MemoryTraits<Kokkos::Unmanaged>,
    Kokkos::Impl::ViewDefault> RV;
  typedef Kokkos::View<const KOKKOSBLAS_IMPL_MV_SCALAR**,
    Kokkos::LayoutLeft,
    Kokkos::Device<KOKKOSBLAS_IMPL_MV_EXEC_SPACE, KOKKOSBLAS_IMPL_MV_MEM_SPACE>,
    Kokkos::MemoryTraits<Kokkos::Unmanaged>,
    Kokkos::Impl::ViewDefault> XMV;
  typedef Kokkos::View<const KOKKOSBLAS_IMPL_MV_SCALAR**,
    Kokkos::LayoutLeft,
    Kokkos::Device<KOKKOSBLAS_IMPL_MV_EXEC_SPACE, KOKKOSBLAS_IMPL_MV_MEM_SPACE>,
    Kokkos::MemoryTraits<Kokkos::Unmanaged>,
    Kokkos::Impl::ViewDefault> YMV;

  static void dot (const RV& r, const XMV& X, const YMV& Y);
};

template<>
struct Dot_MV<Kokkos::View<Kokkos::Details::InnerProductSpaceTraits<KOKKOSBLAS_IMPL_MV_SCALAR>::dot_type,
                           KOKKOSBLAS_IMPL_MV_EXEC_SPACE::array_layout,
                           Kokkos::Device<KOKKOSBLAS_IMPL_MV_EXEC_SPACE, KOKKOSBLAS_IMPL_MV_MEM_SPACE>,
                           Kokkos::MemoryTraits<Kokkos::Unmanaged>,
                           Kokkos::Impl::ViewDefault>,
              Kokkos::View<const KOKKOSBLAS_IMPL_MV_SCALAR*,
                           Kokkos::LayoutLeft,
                           Kokkos::Device<KOKKOSBLAS_IMPL_MV_EXEC_SPACE, KOKKOSBLAS_IMPL_MV_MEM_SPACE>,
                           Kokkos::MemoryTraits<Kokkos::Unmanaged>,
                           Kokkos::Impl::ViewDefault>,
              Kokkos::View<const KOKKOSBLAS_IMPL_MV_SCALAR*,
                           Kokkos::LayoutLeft,
                           Kokkos::Device<KOKKOSBLAS_IMPL_MV_EXEC_SPACE, KOKKOSBLAS_IMPL_MV_MEM_SPACE>,
                           Kokkos::MemoryTraits<Kokkos::Unmanaged>,
                           Kokkos::Impl::ViewDefault>,
              1>
{
  typedef Kokkos::View<Kokkos::Details::InnerProductSpaceTraits<KOKKOSBLAS_IMPL_MV_SCALAR>::dot_type,
    KOKKOSBLAS_IMPL_MV_EXEC_SPACE::array_layout,
    Kokkos::Device<KOKKOSBLAS_IMPL_MV_EXEC_SPACE, KOKKOSBLAS_IMPL_MV_MEM_SPACE>,
    Kokkos::MemoryTraits<Kokkos::Unmanaged>,
    Kokkos::Impl::ViewDefault> RV;
  typedef Kokkos::View<const KOKKOSBLAS_IMPL_MV_SCALAR*,
    Kokkos::LayoutLeft,
    Kokkos::Device<KOKKOSBLAS_IMPL_MV_EXEC_SPACE, KOKKOSBLAS_IMPL_MV_MEM_SPACE>,
    Kokkos::MemoryTraits<Kokkos::Unmanaged>,
    Kokkos::Impl::ViewDefault> XV;
  typedef Kokkos::View<const KOKKOSBLAS_IMPL_MV_SCALAR*,
    Kokkos::LayoutLeft,
    Kokkos::Device<KOKKOSBLAS_IMPL_MV_EXEC_SPACE, KOKKOSBLAS_IMPL_MV_MEM_SPACE>,
    Kokkos::MemoryTraits<Kokkos::Unmanaged>,
    Kokkos::Impl::ViewDefault> YV;

  static void dot (const RV& r, const XV& X, const YV& Y);
};
#undef KOKKOSBLAS_IMPL_MV_EXEC_SPACE
#undef KOKKOSBLAS_IMPL_MV_MEM_SPACE
#undef KOKKOSBLAS_IMPL_MV_SCALAR
#endif // KOKKOS_HAVE_OPENMP


#ifdef KOKKOS_HAVE_PTHREAD
#define KOKKOSBLAS_IMPL_MV_EXEC_SPACE Kokkos::Threads
#define KOKKOSBLAS_IMPL_MV_MEM_SPACE Kokkos::HostSpace
#define KOKKOSBLAS_IMPL_MV_SCALAR double
template<>
struct Dot_MV<Kokkos::View<Kokkos::Details::InnerProductSpaceTraits<KOKKOSBLAS_IMPL_MV_SCALAR>::dot_type*,
                           KOKKOSBLAS_IMPL_MV_EXEC_SPACE::array_layout,
                           Kokkos::Device<KOKKOSBLAS_IMPL_MV_EXEC_SPACE, KOKKOSBLAS_IMPL_MV_MEM_SPACE>,
                           Kokkos::MemoryTraits<Kokkos::Unmanaged>,
                           Kokkos::Impl::ViewDefault>,
              Kokkos::View<const KOKKOSBLAS_IMPL_MV_SCALAR**,
                           Kokkos::LayoutLeft,
                           Kokkos::Device<KOKKOSBLAS_IMPL_MV_EXEC_SPACE, KOKKOSBLAS_IMPL_MV_MEM_SPACE>,
                           Kokkos::MemoryTraits<Kokkos::Unmanaged>,
                           Kokkos::Impl::ViewDefault>,
              Kokkos::View<const KOKKOSBLAS_IMPL_MV_SCALAR**,
                           Kokkos::LayoutLeft,
                           Kokkos::Device<KOKKOSBLAS_IMPL_MV_EXEC_SPACE, KOKKOSBLAS_IMPL_MV_MEM_SPACE>,
                           Kokkos::MemoryTraits<Kokkos::Unmanaged>,
                           Kokkos::Impl::ViewDefault>,
              2>
{
  typedef Kokkos::View<Kokkos::Details::InnerProductSpaceTraits<KOKKOSBLAS_IMPL_MV_SCALAR>::dot_type*,
    KOKKOSBLAS_IMPL_MV_EXEC_SPACE::array_layout,
    Kokkos::Device<KOKKOSBLAS_IMPL_MV_EXEC_SPACE, KOKKOSBLAS_IMPL_MV_MEM_SPACE>,
    Kokkos::MemoryTraits<Kokkos::Unmanaged>,
    Kokkos::Impl::ViewDefault> RV;
  typedef Kokkos::View<const KOKKOSBLAS_IMPL_MV_SCALAR**,
    Kokkos::LayoutLeft,
    Kokkos::Device<KOKKOSBLAS_IMPL_MV_EXEC_SPACE, KOKKOSBLAS_IMPL_MV_MEM_SPACE>,
    Kokkos::MemoryTraits<Kokkos::Unmanaged>,
    Kokkos::Impl::ViewDefault> XMV;
  typedef Kokkos::View<const KOKKOSBLAS_IMPL_MV_SCALAR**,
    Kokkos::LayoutLeft,
    Kokkos::Device<KOKKOSBLAS_IMPL_MV_EXEC_SPACE, KOKKOSBLAS_IMPL_MV_MEM_SPACE>,
    Kokkos::MemoryTraits<Kokkos::Unmanaged>,
    Kokkos::Impl::ViewDefault> YMV;

  static void dot (const RV& r, const XMV& X, const YMV& Y);
};

template<>
struct Dot_MV<Kokkos::View<Kokkos::Details::InnerProductSpaceTraits<KOKKOSBLAS_IMPL_MV_SCALAR>::dot_type,
                           KOKKOSBLAS_IMPL_MV_EXEC_SPACE::array_layout,
                           Kokkos::Device<KOKKOSBLAS_IMPL_MV_EXEC_SPACE, KOKKOSBLAS_IMPL_MV_MEM_SPACE>,
                           Kokkos::MemoryTraits<Kokkos::Unmanaged>,
                           Kokkos::Impl::ViewDefault>,
              Kokkos::View<const KOKKOSBLAS_IMPL_MV_SCALAR*,
                           Kokkos::LayoutLeft,
                           Kokkos::Device<KOKKOSBLAS_IMPL_MV_EXEC_SPACE, KOKKOSBLAS_IMPL_MV_MEM_SPACE>,
                           Kokkos::MemoryTraits<Kokkos::Unmanaged>,
                           Kokkos::Impl::ViewDefault>,
              Kokkos::View<const KOKKOSBLAS_IMPL_MV_SCALAR*,
                           Kokkos::LayoutLeft,
                           Kokkos::Device<KOKKOSBLAS_IMPL_MV_EXEC_SPACE, KOKKOSBLAS_IMPL_MV_MEM_SPACE>,
                           Kokkos::MemoryTraits<Kokkos::Unmanaged>,
                           Kokkos::Impl::ViewDefault>,
              1>
{
  typedef Kokkos::View<Kokkos::Details::InnerProductSpaceTraits<KOKKOSBLAS_IMPL_MV_SCALAR>::dot_type,
    KOKKOSBLAS_IMPL_MV_EXEC_SPACE::array_layout,
    Kokkos::Device<KOKKOSBLAS_IMPL_MV_EXEC_SPACE, KOKKOSBLAS_IMPL_MV_MEM_SPACE>,
    Kokkos::MemoryTraits<Kokkos::Unmanaged>,
    Kokkos::Impl::ViewDefault> RV;
  typedef Kokkos::View<const KOKKOSBLAS_IMPL_MV_SCALAR*,
    Kokkos::LayoutLeft,
    Kokkos::Device<KOKKOSBLAS_IMPL_MV_EXEC_SPACE, KOKKOSBLAS_IMPL_MV_MEM_SPACE>,
    Kokkos::MemoryTraits<Kokkos::Unmanaged>,
    Kokkos::Impl::ViewDefault> XV;
  typedef Kokkos::View<const KOKKOSBLAS_IMPL_MV_SCALAR*,
    Kokkos::LayoutLeft,
    Kokkos::Device<KOKKOSBLAS_IMPL_MV_EXEC_SPACE, KOKKOSBLAS_IMPL_MV_MEM_SPACE>,
    Kokkos::MemoryTraits<Kokkos::Unmanaged>,
    Kokkos::Impl::ViewDefault> YV;

  static void dot (const RV& r, const XV& X, const YV& Y);
};
#undef KOKKOSBLAS_IMPL_MV_EXEC_SPACE
#undef KOKKOSBLAS_IMPL_MV_MEM_SPACE
#undef KOKKOSBLAS_IMPL_MV_SCALAR
#endif // KOKKOS_HAVE_PTHREAD


#ifdef KOKKOS_HAVE_CUDA
#define KOKKOSBLAS_IMPL_MV_EXEC_SPACE Kokkos::Cuda
#define KOKKOSBLAS_IMPL_MV_MEM_SPACE Kokkos::CudaSpace
#define KOKKOSBLAS_IMPL_MV_SCALAR double
template<>
struct Dot_MV<Kokkos::View<Kokkos::Details::InnerProductSpaceTraits<KOKKOSBLAS_IMPL_MV_SCALAR>::dot_type*,
                           KOKKOSBLAS_IMPL_MV_EXEC_SPACE::array_layout,
                           Kokkos::Device<KOKKOSBLAS_IMPL_MV_EXEC_SPACE, KOKKOSBLAS_IMPL_MV_MEM_SPACE>,
                           Kokkos::MemoryTraits<Kokkos::Unmanaged>,
                           Kokkos::Impl::ViewDefault>,
              Kokkos::View<const KOKKOSBLAS_IMPL_MV_SCALAR**,
                           Kokkos::LayoutLeft,
                           Kokkos::Device<KOKKOSBLAS_IMPL_MV_EXEC_SPACE, KOKKOSBLAS_IMPL_MV_MEM_SPACE>,
                           Kokkos::MemoryTraits<Kokkos::Unmanaged>,
                           Kokkos::Impl::ViewDefault>,
              Kokkos::View<const KOKKOSBLAS_IMPL_MV_SCALAR**,
                           Kokkos::LayoutLeft,
                           Kokkos::Device<KOKKOSBLAS_IMPL_MV_EXEC_SPACE, KOKKOSBLAS_IMPL_MV_MEM_SPACE>,
                           Kokkos::MemoryTraits<Kokkos::Unmanaged>,
                           Kokkos::Impl::ViewDefault>,
              2>
{
  typedef Kokkos::View<Kokkos::Details::InnerProductSpaceTraits<KOKKOSBLAS_IMPL_MV_SCALAR>::dot_type*,
    KOKKOSBLAS_IMPL_MV_EXEC_SPACE::array_layout,
    Kokkos::Device<KOKKOSBLAS_IMPL_MV_EXEC_SPACE, KOKKOSBLAS_IMPL_MV_MEM_SPACE>,
    Kokkos::MemoryTraits<Kokkos::Unmanaged>,
    Kokkos::Impl::ViewDefault> RV;
  typedef Kokkos::View<const KOKKOSBLAS_IMPL_MV_SCALAR**,
    Kokkos::LayoutLeft,
    Kokkos::Device<KOKKOSBLAS_IMPL_MV_EXEC_SPACE, KOKKOSBLAS_IMPL_MV_MEM_SPACE>,
    Kokkos::MemoryTraits<Kokkos::Unmanaged>,
    Kokkos::Impl::ViewDefault> XMV;
  typedef Kokkos::View<const KOKKOSBLAS_IMPL_MV_SCALAR**,
    Kokkos::LayoutLeft,
    Kokkos::Device<KOKKOSBLAS_IMPL_MV_EXEC_SPACE, KOKKOSBLAS_IMPL_MV_MEM_SPACE>,
    Kokkos::MemoryTraits<Kokkos::Unmanaged>,
    Kokkos::Impl::ViewDefault> YMV;

  static void dot (const RV& r, const XMV& X, const YMV& Y);
};

template<>
struct Dot_MV<Kokkos::View<Kokkos::Details::InnerProductSpaceTraits<KOKKOSBLAS_IMPL_MV_SCALAR>::dot_type,
                           KOKKOSBLAS_IMPL_MV_EXEC_SPACE::array_layout,
                           Kokkos::Device<KOKKOSBLAS_IMPL_MV_EXEC_SPACE, KOKKOSBLAS_IMPL_MV_MEM_SPACE>,
                           Kokkos::MemoryTraits<Kokkos::Unmanaged>,
                           Kokkos::Impl::ViewDefault>,
              Kokkos::View<const KOKKOSBLAS_IMPL_MV_SCALAR*,
                           Kokkos::LayoutLeft,
                           Kokkos::Device<KOKKOSBLAS_IMPL_MV_EXEC_SPACE, KOKKOSBLAS_IMPL_MV_MEM_SPACE>,
                           Kokkos::MemoryTraits<Kokkos::Unmanaged>,
                           Kokkos::Impl::ViewDefault>,
              Kokkos::View<const KOKKOSBLAS_IMPL_MV_SCALAR*,
                           Kokkos::LayoutLeft,
                           Kokkos::Device<KOKKOSBLAS_IMPL_MV_EXEC_SPACE, KOKKOSBLAS_IMPL_MV_MEM_SPACE>,
                           Kokkos::MemoryTraits<Kokkos::Unmanaged>,
                           Kokkos::Impl::ViewDefault>,
              1>
{
  typedef Kokkos::View<Kokkos::Details::InnerProductSpaceTraits<KOKKOSBLAS_IMPL_MV_SCALAR>::dot_type,
    KOKKOSBLAS_IMPL_MV_EXEC_SPACE::array_layout,
    Kokkos::Device<KOKKOSBLAS_IMPL_MV_EXEC_SPACE, KOKKOSBLAS_IMPL_MV_MEM_SPACE>,
    Kokkos::MemoryTraits<Kokkos::Unmanaged>,
    Kokkos::Impl::ViewDefault> RV;
  typedef Kokkos::View<const KOKKOSBLAS_IMPL_MV_SCALAR*,
    Kokkos::LayoutLeft,
    Kokkos::Device<KOKKOSBLAS_IMPL_MV_EXEC_SPACE, KOKKOSBLAS_IMPL_MV_MEM_SPACE>,
    Kokkos::MemoryTraits<Kokkos::Unmanaged>,
    Kokkos::Impl::ViewDefault> XV;
  typedef Kokkos::View<const KOKKOSBLAS_IMPL_MV_SCALAR*,
    Kokkos::LayoutLeft,
    Kokkos::Device<KOKKOSBLAS_IMPL_MV_EXEC_SPACE, KOKKOSBLAS_IMPL_MV_MEM_SPACE>,
    Kokkos::MemoryTraits<Kokkos::Unmanaged>,
    Kokkos::Impl::ViewDefault> YV;

  static void dot (const RV& r, const XV& X, const YV& Y);
};
#undef KOKKOSBLAS_IMPL_MV_EXEC_SPACE
#undef KOKKOSBLAS_IMPL_MV_MEM_SPACE
#undef KOKKOSBLAS_IMPL_MV_SCALAR
#endif // KOKKOS_HAVE_CUDA


#ifdef KOKKOS_HAVE_CUDA
#define KOKKOSBLAS_IMPL_MV_EXEC_SPACE Kokkos::Cuda
#define KOKKOSBLAS_IMPL_MV_MEM_SPACE Kokkos::CudaUVMSpace
#define KOKKOSBLAS_IMPL_MV_SCALAR double
template<>
struct Dot_MV<Kokkos::View<Kokkos::Details::InnerProductSpaceTraits<KOKKOSBLAS_IMPL_MV_SCALAR>::dot_type*,
                           KOKKOSBLAS_IMPL_MV_EXEC_SPACE::array_layout,
                           Kokkos::Device<KOKKOSBLAS_IMPL_MV_EXEC_SPACE, KOKKOSBLAS_IMPL_MV_MEM_SPACE>,
                           Kokkos::MemoryTraits<Kokkos::Unmanaged>,
                           Kokkos::Impl::ViewDefault>,
              Kokkos::View<const KOKKOSBLAS_IMPL_MV_SCALAR**,
                           Kokkos::LayoutLeft,
                           Kokkos::Device<KOKKOSBLAS_IMPL_MV_EXEC_SPACE, KOKKOSBLAS_IMPL_MV_MEM_SPACE>,
                           Kokkos::MemoryTraits<Kokkos::Unmanaged>,
                           Kokkos::Impl::ViewDefault>,
              Kokkos::View<const KOKKOSBLAS_IMPL_MV_SCALAR**,
                           Kokkos::LayoutLeft,
                           Kokkos::Device<KOKKOSBLAS_IMPL_MV_EXEC_SPACE, KOKKOSBLAS_IMPL_MV_MEM_SPACE>,
                           Kokkos::MemoryTraits<Kokkos::Unmanaged>,
                           Kokkos::Impl::ViewDefault>,
              2>
{
  typedef Kokkos::View<Kokkos::Details::InnerProductSpaceTraits<KOKKOSBLAS_IMPL_MV_SCALAR>::dot_type*,
    KOKKOSBLAS_IMPL_MV_EXEC_SPACE::array_layout,
    Kokkos::Device<KOKKOSBLAS_IMPL_MV_EXEC_SPACE, KOKKOSBLAS_IMPL_MV_MEM_SPACE>,
    Kokkos::MemoryTraits<Kokkos::Unmanaged>,
    Kokkos::Impl::ViewDefault> RV;
  typedef Kokkos::View<const KOKKOSBLAS_IMPL_MV_SCALAR**,
    Kokkos::LayoutLeft,
    Kokkos::Device<KOKKOSBLAS_IMPL_MV_EXEC_SPACE, KOKKOSBLAS_IMPL_MV_MEM_SPACE>,
    Kokkos::MemoryTraits<Kokkos::Unmanaged>,
    Kokkos::Impl::ViewDefault> XMV;
  typedef Kokkos::View<const KOKKOSBLAS_IMPL_MV_SCALAR**,
    Kokkos::LayoutLeft,
    Kokkos::Device<KOKKOSBLAS_IMPL_MV_EXEC_SPACE, KOKKOSBLAS_IMPL_MV_MEM_SPACE>,
    Kokkos::MemoryTraits<Kokkos::Unmanaged>,
    Kokkos::Impl::ViewDefault> YMV;

  static void dot (const RV& r, const XMV& X, const YMV& Y);
};

template<>
struct Dot_MV<Kokkos::View<Kokkos::Details::InnerProductSpaceTraits<KOKKOSBLAS_IMPL_MV_SCALAR>::dot_type,
                           KOKKOSBLAS_IMPL_MV_EXEC_SPACE::array_layout,
                           Kokkos::Device<KOKKOSBLAS_IMPL_MV_EXEC_SPACE, KOKKOSBLAS_IMPL_MV_MEM_SPACE>,
                           Kokkos::MemoryTraits<Kokkos::Unmanaged>,
                           Kokkos::Impl::ViewDefault>,
              Kokkos::View<const KOKKOSBLAS_IMPL_MV_SCALAR*,
                           Kokkos::LayoutLeft,
                           Kokkos::Device<KOKKOSBLAS_IMPL_MV_EXEC_SPACE, KOKKOSBLAS_IMPL_MV_MEM_SPACE>,
                           Kokkos::MemoryTraits<Kokkos::Unmanaged>,
                           Kokkos::Impl::ViewDefault>,
              Kokkos::View<const KOKKOSBLAS_IMPL_MV_SCALAR*,
                           Kokkos::LayoutLeft,
                           Kokkos::Device<KOKKOSBLAS_IMPL_MV_EXEC_SPACE, KOKKOSBLAS_IMPL_MV_MEM_SPACE>,
                           Kokkos::MemoryTraits<Kokkos::Unmanaged>,
                           Kokkos::Impl::ViewDefault>,
              1>
{
  typedef Kokkos::View<Kokkos::Details::InnerProductSpaceTraits<KOKKOSBLAS_IMPL_MV_SCALAR>::dot_type,
    KOKKOSBLAS_IMPL_MV_EXEC_SPACE::array_layout,
    Kokkos::Device<KOKKOSBLAS_IMPL_MV_EXEC_SPACE, KOKKOSBLAS_IMPL_MV_MEM_SPACE>,
    Kokkos::MemoryTraits<Kokkos::Unmanaged>,
    Kokkos::Impl::ViewDefault> RV;
  typedef Kokkos::View<const KOKKOSBLAS_IMPL_MV_SCALAR*,
    Kokkos::LayoutLeft,
    Kokkos::Device<KOKKOSBLAS_IMPL_MV_EXEC_SPACE, KOKKOSBLAS_IMPL_MV_MEM_SPACE>,
    Kokkos::MemoryTraits<Kokkos::Unmanaged>,
    Kokkos::Impl::ViewDefault> XV;
  typedef Kokkos::View<const KOKKOSBLAS_IMPL_MV_SCALAR*,
    Kokkos::LayoutLeft,
    Kokkos::Device<KOKKOSBLAS_IMPL_MV_EXEC_SPACE, KOKKOSBLAS_IMPL_MV_MEM_SPACE>,
    Kokkos::MemoryTraits<Kokkos::Unmanaged>,
    Kokkos::Impl::ViewDefault> YV;

  static void dot (const RV& r, const XV& X, const YV& Y);
};
#undef KOKKOSBLAS_IMPL_MV_EXEC_SPACE
#undef KOKKOSBLAS_IMPL_MV_MEM_SPACE
#undef KOKKOSBLAS_IMPL_MV_SCALAR
#endif // KOKKOS_HAVE_CUDA


} // namespace Impl
} // namespace KokkosBlas

#endif // KOKKOS_BLAS1_MV_IMPL_DOT_HPP_
