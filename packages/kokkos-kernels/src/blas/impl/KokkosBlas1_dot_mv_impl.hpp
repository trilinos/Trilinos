/*
//@HEADER
// ************************************************************************
//
//               KokkosKernels 0.9: Linear Algebra and Graph Kernels
//                 Copyright 2017 Sandia Corporation
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
// Questions? Contact Siva Rajamanickam (srajama@sandia.gov)
//
// ************************************************************************
//@HEADER
*/
#ifndef KOKKOSBLAS1_IMPL_DOT_MV_IMPL_HPP_
#define KOKKOSBLAS1_IMPL_DOT_MV_IMPL_HPP_

#ifndef KOKKOSBLAS_OPTIMIZATION_LEVEL_DOT
#define KOKKOSBLAS_OPTIMIZATION_LEVEL_DOT 2
#endif // KOKKOSBLAS_OPTIMIZATION_LEVEL_DOT

#include <KokkosKernels_config.h>
#include <Kokkos_Core.hpp>
#include <Kokkos_InnerProductSpaceTraits.hpp>
#include <type_traits>
#include <KokkosBlas1_dot_impl.hpp>

namespace KokkosBlas {
namespace Impl {


/// \brief Dot product functor for a multivector times a single
///   vector, or vice versa.
///
/// "Multivector dot a single vector" means that each column of the
/// multivector gets dotted with the same vector, as if the latter
/// vector were replicated.  "Single vector dot a multivector" means
/// the (conjugate) transpose of that.  We combine both (multivector
/// dot vector) and (vector dot multivector) cases into a single
/// functor, to avoid code duplication.
///
/// \tparam RV 1-D output View
/// \tparam XMV 2-D input View (the multivector)
/// \tparam YV 1-D input View (the single vector)
/// \tparam SizeType Index type.  Use int (32 bits) if possible.
template<class RV, class XMV, class YV, class SizeType = typename XMV::size_type>
struct MV_V_Dot_Functor
{
  typedef typename XMV::execution_space              execution_space;
  typedef SizeType                                         size_type;
  typedef typename XMV::non_const_value_type             xvalue_type;
  typedef Kokkos::Details::InnerProductSpaceTraits<xvalue_type>  IPT;
  typedef Kokkos::Details::ArithTraits<typename IPT::dot_type>    AT;
  typedef typename IPT::dot_type                        value_type[];

  size_type value_count;
  RV m_r;
  typename XMV::const_type m_x;
  typename YV::const_type m_y;
  //! If true, do y dot x instead of x dot y.
  bool reverseOrder_;

  MV_V_Dot_Functor (const RV& r, const XMV& x, const YV& y,
                    const bool reverseOrder) :
    value_count (x.extent(1)), m_r (r), m_x (x), m_y (y),
    reverseOrder_ (reverseOrder)
  {
    static_assert (Kokkos::Impl::is_view<RV>::value, "KokkosBlas::Impl::"
                   "MV_V_Dot_Functor: R is not a Kokkos::View.");
    static_assert (Kokkos::Impl::is_view<XMV>::value, "KokkosBlas::Impl::"
                   "MV_V_Dot_Functor: X is not a Kokkos::View.");
    static_assert (Kokkos::Impl::is_view<YV>::value, "KokkosBlas::Impl::"
                   "MV_V_Dot_Functor: Y is not a Kokkos::View.");
    static_assert (Kokkos::Impl::is_same<typename RV::value_type,
                   typename RV::non_const_value_type>::value,
                   "KokkosBlas::Impl::MV_V_Dot_Functor: R is const.  "
                   "It must be nonconst, because it is an output argument "
                   "(we have to be able to write to its entries).");
    static_assert (static_cast<int> (XMV::rank) == 2,
                   "KokkosBlas::Impl::MV_V_Dot_Functor: X must have rank 2.");
    static_assert (static_cast<int> (YV::rank) == 1,
                   "KokkosBlas::Impl::MV_V_Dot_Functor: Y must have rank 1.");
    static_assert (static_cast<int> (RV::rank) == 1,
                   "KokkosBlas::Impl::MV_V_Dot_Functor: RV must have rank 1.");
  }

  KOKKOS_INLINE_FUNCTION
  void operator() (const size_type& i, value_type sum) const
  {
    const size_type numVecs = value_count;
    if (reverseOrder_) {
#ifdef KOKKOS_ENABLE_PRAGMA_IVDEP
#pragma ivdep
#endif
#ifdef KOKKOS_ENABLE_PRAGMA_VECTOR
#pragma vector always
#endif
      for (size_type k = 0; k < numVecs; ++k) {
        sum[k] += IPT::dot (m_y(i), m_x(i,k)); // m_x(i,k) * m_y(i)
      }
    }
    else {
#ifdef KOKKOS_ENABLE_PRAGMA_IVDEP
#pragma ivdep
#endif
#ifdef KOKKOS_ENABLE_PRAGMA_VECTOR
#pragma vector always
#endif
      for (size_type k = 0; k < numVecs; ++k) {
        sum[k] += IPT::dot (m_x(i,k), m_y(i)); // m_x(i,k) * m_y(i)
      }
    }
  }

  KOKKOS_INLINE_FUNCTION void init (value_type update) const
  {
    const size_type numVecs = value_count;
#ifdef KOKKOS_ENABLE_PRAGMA_IVDEP
#pragma ivdep
#endif
#ifdef KOKKOS_ENABLE_PRAGMA_VECTOR
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
#ifdef KOKKOS_ENABLE_PRAGMA_IVDEP
#pragma ivdep
#endif
#ifdef KOKKOS_ENABLE_PRAGMA_VECTOR
#pragma vector always
#endif
    for (size_type k = 0; k < numVecs; ++k) {
      update[k] += source[k];
    }
  }

  // On device, write the reduction result to the output View.
  /*KOKKOS_INLINE_FUNCTION void
  final (const value_type dst) const
  {
    const size_type numVecs = value_count;
#ifdef KOKKOS_ENABLE_PRAGMA_IVDEP
#pragma ivdep
#endif
#ifdef KOKKOS_ENABLE_PRAGMA_VECTOR
#pragma vector always
#endif
    for (size_type k = 0; k < numVecs; ++k) {
      m_r(k) = dst[k];
    }
  }*/
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
    value_count (x.extent(1)), m_r (r), m_x (x), m_y (y)
  {
    static_assert (Kokkos::Impl::is_view<RV>::value, "KokkosBlas::Impl::"
                   "MV_Dot_Right_FunctorVector: R is not a Kokkos::View.");
    static_assert (Kokkos::Impl::is_view<XMV>::value, "KokkosBlas::Impl::"
                   "MV_Dot_Right_FunctorVector: X is not a Kokkos::View.");
    static_assert (Kokkos::Impl::is_view<YMV>::value, "KokkosBlas::Impl::"
                   "MV_Dot_Right_FunctorVector: Y is not a Kokkos::View.");
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
  }

  KOKKOS_INLINE_FUNCTION
  void operator() (const size_type& i, value_type sum) const
  {
    const size_type numVecs = value_count;
#ifdef KOKKOS_ENABLE_PRAGMA_IVDEP
#pragma ivdep
#endif
#ifdef KOKKOS_ENABLE_PRAGMA_VECTOR
#pragma vector always
#endif
    for (size_type k = 0; k < numVecs; ++k) {
      sum[k] += IPT::dot (m_x(i,k), m_y(i,k)); // m_x(i,k) * m_y(i,k)
    }
  }

  KOKKOS_INLINE_FUNCTION void init (value_type update) const
  {
    const size_type numVecs = value_count;
#ifdef KOKKOS_ENABLE_PRAGMA_IVDEP
#pragma ivdep
#endif
#ifdef KOKKOS_ENABLE_PRAGMA_VECTOR
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
#ifdef KOKKOS_ENABLE_PRAGMA_IVDEP
#pragma ivdep
#endif
#ifdef KOKKOS_ENABLE_PRAGMA_VECTOR
#pragma vector always
#endif
    for (size_type k = 0; k < numVecs; ++k) {
      update[k] += source[k];
    }
  }

  // On device, write the reduction result to the output View.
  /*KOKKOS_INLINE_FUNCTION void
  final (const value_type dst) const
  {
    const size_type numVecs = value_count;
#ifdef KOKKOS_ENABLE_PRAGMA_IVDEP
#pragma ivdep
#endif
#ifdef KOKKOS_ENABLE_PRAGMA_VECTOR
#pragma vector always
#endif
    for (size_type k = 0; k < numVecs; ++k) {
      m_r(k) = dst[k];
    }
  }*/
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
    value_count (x.extent(1)), m_r (r), m_x (x), m_y (y)
  {
    static_assert (Kokkos::Impl::is_view<RV>::value, "KokkosBlas::Impl::"
                   "MV_Dot_Right_FunctorUnroll: R is not a Kokkos::View.");
    static_assert (Kokkos::Impl::is_view<XMV>::value, "KokkosBlas::Impl::"
                   "MV_Dot_Right_FunctorUnroll: X is not a Kokkos::View.");
    static_assert (Kokkos::Impl::is_view<YMV>::value, "KokkosBlas::Impl::"
                   "MV_Dot_Right_FunctorUnroll: Y is not a Kokkos::View.");
    static_assert (Kokkos::Impl::is_same<typename RV::value_type,
                   typename RV::non_const_value_type>::value,
                   "KokkosBlas::Impl::MV_Dot_Right_FunctorUnroll: R is const.  "
                   "It must be nonconst, because it is an output argument "
                   "(we have to be able to write to its entries).");
    static_assert (int(XMV::rank) == int(YMV::rank),
                   "KokkosBlas::Impl::MV_Dot_Right_FunctorUnroll: "
                   "X and Y must have the same rank.");
    static_assert (RV::rank == 1 && XMV::rank == 2,
                   "KokkosBlas::Impl::MV_Dot_Right_FunctorUnroll: "
                   "RV must have rank 1 and XMV and YMV must have rank 2.");
  }

  KOKKOS_INLINE_FUNCTION
  void operator() (const size_type& i, value_type sum) const
  {
#ifdef KOKKOS_ENABLE_PRAGMA_UNROLL
#pragma unroll
#endif
    for (int k = 0; k < UNROLL; ++k) {
      sum[k] += IPT::dot (m_x(i,k), m_y(i,k)); // m_x(i,k) * m_y(i,k)
    }
  }

  KOKKOS_INLINE_FUNCTION void init (volatile value_type update) const
  {
#ifdef KOKKOS_ENABLE_PRAGMA_UNROLL
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
#ifdef KOKKOS_ENABLE_PRAGMA_UNROLL
#pragma unroll
#endif
    for (int k = 0; k < UNROLL; ++k) {
      update[k] += source[k];
    }
  }

  // On device, write the reduction result to the output View.
  /*KOKKOS_INLINE_FUNCTION void
  final (const value_type dst) const
  {
#ifdef KOKKOS_ENABLE_PRAGMA_UNROLL
#pragma unroll
#endif
    for (int k = 0; k < UNROLL; ++k) {
      m_r(k) = dst[k];
    }
  }*/
};


//! Implementation detail of MV_V_Dot_Invoke (see below).
template<class RV, class XMV, class YMV, class SizeType,
         const int XMV_rank = XMV::rank,
         const int YMV_rank = YMV::rank>
struct MV_V_Dot_Invoke_Impl
{
  static void
  run (const RV& r, const XMV& X, const YMV& Y, const SizeType numRows);
};

template<class RV, class XMV, class YMV, class SizeType>
struct MV_V_Dot_Invoke_Impl<RV, XMV, YMV, SizeType, 2, 1>
{
  static void
  run (const RV& r, const XMV& X, const YMV& Y, const SizeType numRows)
  {
    static_assert (static_cast<int> (XMV::rank) == 2 && static_cast<int> (YMV::rank) == 1,
                   "XMV must have rank 2, and YMV must have rank 1.");
    static_assert (static_cast<int> (RV::rank) == 1, "Rank of r must be 1.");
    static_assert (std::is_integral<SizeType>::value,
                   "SizeType must be a built-in integer type.");

    typedef typename XMV::execution_space execution_space;
    typedef SizeType size_type;
    typedef Kokkos::RangePolicy<execution_space, size_type> range_type;

    typedef MV_V_Dot_Functor<RV, XMV, YMV, size_type> op_type;
    constexpr bool reverseOrder = false;
    op_type op (r, X, Y, reverseOrder);
    Kokkos::parallel_reduce (range_type (0, numRows), op, r);
  }
};

template<class RV, class XMV, class YMV, class SizeType>
struct MV_V_Dot_Invoke_Impl<RV, XMV, YMV, SizeType, 1, 2>
{
  static void
  run (const RV& r, const XMV& X, const YMV& Y, const SizeType numRows)
  {
    static_assert (static_cast<int> (XMV::rank) == 1 && static_cast<int> (YMV::rank) == 2,
                   "XMV must have rank 1, and YMV must have rank 2.");
    static_assert (static_cast<int> (RV::rank) == 1, "Rank of r must be 1.");
    static_assert (std::is_integral<SizeType>::value,
                   "SizeType must be a built-in integer type.");

    typedef typename XMV::execution_space execution_space;
    typedef SizeType size_type;
    typedef Kokkos::RangePolicy<execution_space, size_type> range_type;

    // Use the "reverse arguments" mode of the functor.
    typedef MV_V_Dot_Functor<RV, YMV, XMV, size_type> op_type;
    constexpr bool reverseOrder = true;
    op_type op (r, Y, X, reverseOrder);
    Kokkos::parallel_reduce (range_type (0, numRows), op, r);
  }
};

//! Special case where XMV has rank 2 and YMV has rank 1, or vice versa.
template<class RV, class XMV, class YMV, class SizeType>
void
MV_V_Dot_Invoke (const RV& r, const XMV& X, const YMV& Y, const SizeType numRows)
{
  MV_V_Dot_Invoke_Impl<RV, XMV, YMV, SizeType>::run (r, X, Y, numRows);
}

//! Special case where XMV and YMV both have rank 2.
template<class RV, class XMV, class YMV, class SizeType>
void
MV_Dot_Invoke (const RV& r, const XMV& X, const YMV& Y)
{
  const SizeType numRows = static_cast<SizeType> (X.extent(0));
  const SizeType numCols = static_cast<SizeType> (X.extent(1));
  Kokkos::RangePolicy<typename XMV::execution_space, SizeType> policy (0, numRows);

  if (static_cast<int> (X.extent(1)) != 1 && static_cast<int> (Y.extent(1)) == 1) {
    // X has > 1 columns, and Y has 1 column.
    auto Y_0 = Kokkos::subview (Y, Kokkos::ALL (), 0);
    typedef typename decltype (Y_0)::const_type YV;
    MV_V_Dot_Invoke<RV, XMV, YV, SizeType> (r, X, Y_0, numRows);
    return;
  }
  else if (static_cast<int> (X.extent(1)) == 1 && static_cast<int> (Y.extent(1)) != 1) {
    // X has 1 column, and Y has > 1 columns.
    auto X_0 = Kokkos::subview (X, Kokkos::ALL (), 0);
    typedef typename decltype (X_0)::const_type XV;
    MV_V_Dot_Invoke<RV, XV, YMV, SizeType> (r, X_0, Y, numRows);
    return;
  }

#if KOKKOSBLAS_OPTIMIZATION_LEVEL_DOT <= 2

  // Strip-mine by 8, then 4.  After that, do one column at a time.
  // We limit the number of strip-mine values in order to keep down
  // the amount of code to compile.

  SizeType j = 0; // the current column of X and Y
  for ( ; j + 8 <= numCols; j += 8) {
    auto X_cur = Kokkos::subview (X, Kokkos::ALL (), std::make_pair (j, j+8));
    auto Y_cur = Kokkos::subview (Y, Kokkos::ALL (), std::make_pair (j, j+8));
    auto r_cur = Kokkos::subview (r, std::make_pair (j, j+8));

    MV_Dot_Right_FunctorUnroll<RV, XMV, YMV, 8, SizeType> op (r_cur, X_cur, Y_cur);
    Kokkos::parallel_reduce (policy, op, r_cur);
  }
  for ( ; j + 4 <= numCols; j += 4) {
    auto X_cur = Kokkos::subview (X, Kokkos::ALL (), std::make_pair (j, j+4));
    auto Y_cur = Kokkos::subview (Y, Kokkos::ALL (), std::make_pair (j, j+4));
    auto r_cur = Kokkos::subview (r, std::make_pair (j, j+4));

    MV_Dot_Right_FunctorUnroll<RV, XMV, YMV, 4, SizeType> op (r_cur, X_cur, Y_cur);
    Kokkos::parallel_reduce (policy, op, r_cur);
  }
  for ( ; j < numCols; ++j) {
    // RV needs to turn 0-D, and XMV and YMV need to turn 1-D.
    auto x_cur = Kokkos::subview (X, Kokkos::ALL (), j);
    auto y_cur = Kokkos::subview (Y, Kokkos::ALL (), j);
    auto r_cur = Kokkos::subview (r, j);
    typedef decltype (r_cur) RV0D;
    typedef decltype (x_cur) XMV1D;
    typedef decltype (y_cur) YMV1D;

    DotFunctor<RV0D, XMV1D, YMV1D, SizeType> op(x_cur, y_cur);
    Kokkos::parallel_reduce (policy, op, r_cur);
  }

#else // KOKKOSBLAS_OPTIMIZATION_LEVEL_DOT > 2

  if (numCols > 16) {
    MV_Dot_Right_FunctorVector<RV, XMV, YMV, SizeType> op (r, X, Y);
    Kokkos::parallel_reduce (policy, op, r);
  }
  else {
    switch (numCols) {
    case 16: {
      MV_Dot_Right_FunctorUnroll<RV, XMV, YMV, 16, SizeType> op (r, X, Y);
      Kokkos::parallel_reduce (policy, op, r);
      break;
    }
    case 15: {
      MV_Dot_Right_FunctorUnroll<RV, XMV, YMV, 15, SizeType> op (r, X, Y);
      Kokkos::parallel_reduce (policy, op, r);
      break;
    }
    case 14: {
      MV_Dot_Right_FunctorUnroll<RV, XMV, YMV, 14, SizeType> op (r, X, Y);
      Kokkos::parallel_reduce (policy, op, r);
      break;
    }
    case 13: {
      MV_Dot_Right_FunctorUnroll<RV, XMV, YMV, 13, SizeType> op (r, X, Y);
      Kokkos::parallel_reduce (policy, op, r);
      break;
    }
    case 12: {
      MV_Dot_Right_FunctorUnroll<RV, XMV, YMV, 12, SizeType> op (r, X, Y);
      Kokkos::parallel_reduce (policy, op, r);
      break;
    }
    case 11: {
      MV_Dot_Right_FunctorUnroll<RV, XMV, YMV, 11, SizeType> op (r, X, Y);
      Kokkos::parallel_reduce (policy, op, r);
      break;
    }
    case 10: {
      MV_Dot_Right_FunctorUnroll<RV, XMV, YMV, 10, SizeType> op (r, X, Y);
      Kokkos::parallel_reduce (policy, op, r);
      break;
    }
    case 9: {
      MV_Dot_Right_FunctorUnroll<RV, XMV, YMV, 9, SizeType> op (r, X, Y);
      Kokkos::parallel_reduce (policy, op, r);
      break;
    }
    case 8: {
      MV_Dot_Right_FunctorUnroll<RV, XMV, YMV, 8, SizeType> op (r, X, Y);
      Kokkos::parallel_reduce (policy, op, r);
      break;
    }
    case 7: {
      MV_Dot_Right_FunctorUnroll<RV, XMV, YMV, 7, SizeType> op (r, X, Y);
      Kokkos::parallel_reduce (policy, op, r);
      break;
    }
    case 6: {
      MV_Dot_Right_FunctorUnroll<RV, XMV, YMV, 6, SizeType> op (r, X, Y);
      Kokkos::parallel_reduce (policy, op, r);
      break;
    }
    case 5: {
      MV_Dot_Right_FunctorUnroll<RV, XMV, YMV, 5, SizeType> op (r, X, Y);
      Kokkos::parallel_reduce (policy, op, r);
      break;
    }
    case 4: {
      MV_Dot_Right_FunctorUnroll<RV, XMV, YMV, 4, SizeType> op (r, X, Y);
      Kokkos::parallel_reduce (policy, op, r);
      break;
    }
    case 3: {
      MV_Dot_Right_FunctorUnroll<RV, XMV, YMV, 3, SizeType> op (r, X, Y);
      Kokkos::parallel_reduce (policy, op, r);
      break;
    }
    case 2: {
      MV_Dot_Right_FunctorUnroll<RV, XMV, YMV, 2, SizeType> op (r, X, Y);
      Kokkos::parallel_reduce (policy, op, r);
      break;
    }
    case 1: {
      // RV needs to turn 0-D, and XMV and YMV need to turn 1-D.
      auto r_0 = Kokkos::subview (r, 0);
      auto X_0 = Kokkos::subview (X, Kokkos::ALL (), 0);
      auto Y_0 = Kokkos::subview (Y, Kokkos::ALL (), 0);
      typedef decltype (r_0) RV0D;
      typedef decltype (X_0) XMV1D;
      typedef decltype (Y_0) YMV1D;

      typedef V_Dot_Functor<RV0D, XMV1D, YMV1D, SizeType> op_type;
      op_type op (r_0, X_0, Y_0);
      Kokkos::parallel_reduce (policy, op, r_0);
      break;
    }
    } // switch
  } // if-else

#endif // KOKKOSBLAS_OPTIMIZATION_LEVEL_DOT
}

/// \brief Implementation of KokkosBlas::dot for multivectors or
///   single vectors.
///
/// \tparam RV Type of return View of dot product results.
/// \tparam XMV Type of first input (multi)vector X View.
/// \tparam YMV Type of second input (multi)vector Y View.
/// \tparam XMV_rank The rank of XMY.  If 2, it is a multivector; if
///   1, it is a single vector.
/// \tparam YMV_rank The rank of YMV.  If 2, it is a multivector; if
///   1, it is a single vector.
template<class RV, class XMV, class YMV, class SizeType,
         const int XMV_rank = XMV::rank,
         const int YMV_rnak = YMV::rank>
struct Dot_MV;

template<class RV, class XMV, class YMV, class SizeType>
struct Dot_MV<RV, XMV, YMV, SizeType, 2, 2>
{
  /// \brief Compute the dot product(s) of the column(s) of the
  ///   multivectors (2-D views) X and Y, and store result(s) in the
  ///   1-D View r.
  static void dot (const RV& r, const XMV& X, const YMV& Y)
  {
      MV_Dot_Invoke<RV, XMV, YMV, SizeType> (r, X, Y);
  }
};

/// \brief Partial specialization for XMV_rank == 2 and YMV_rank == 1
///   (X is a multivector, and Y is a single column).
template<class RV, class XMV, class YV, class SizeType>
struct Dot_MV<RV, XMV, YV, SizeType, 2, 1> {
  /// \brief Compute the dot product(s) of each column of X with the
  ///   single vector Y, and store result(s) in the 1-D View r.
  static void dot (const RV& r, const XMV& X, const YV& Y)
  {
    MV_V_Dot_Invoke<RV, XMV, YV, SizeType> (r, X, Y, static_cast<int> (X.extent(0)));
  }
};

/// \brief Partial specialization for XMV_rank == 1 and YMV_rank == 2
///   (X is a single column, and Y is a multivector).
template<class RV, class XV, class YMV, class SizeType>
struct Dot_MV<RV, XV, YMV, SizeType, 1, 2> {
  /// \brief Compute the dot product(s) of the single vector X with
  ///   each column of Y, and store result(s) in the 1-D View r.
  static void dot (const RV& r, const XV& X, const YMV& Y)
  {
    const SizeType numRows = X.extent(0);
    MV_V_Dot_Invoke<RV, XV, YMV, SizeType> (r, X, Y, numRows);
  }
};

} // namespace Impl
} // namespace KokkosBlas

#endif // KOKKOSBLAS1_IMPL_DOT_MV_IMPL_HPP_
