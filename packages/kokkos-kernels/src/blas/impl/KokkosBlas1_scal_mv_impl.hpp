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
#ifndef KOKKOSBLAS1_SCAL_MV_IMPL_HPP_
#define KOKKOSBLAS1_SCAL_MV_IMPL_HPP_

#include <KokkosKernels_config.h>
#include <Kokkos_Core.hpp>
#include <Kokkos_InnerProductSpaceTraits.hpp>
#include <KokkosBlas1_scal_spec.hpp>
#include <KokkosBlas1_scal_impl.hpp>

#ifndef KOKKOSBLAS_OPTIMIZATION_LEVEL_SCAL
#define KOKKOSBLAS_OPTIMIZATION_LEVEL_SCAL 2
#endif // KOKKOSBLAS_OPTIMIZATION_LEVEL_SCAL

namespace KokkosBlas {
namespace Impl {

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
  typedef SizeType size_type;
  typedef Kokkos::Details::ArithTraits<typename RMV::non_const_value_type> ATS;

  const size_type numCols;
  RMV R_;
  XMV X_;
  aVector a_;

  MV_Scal_Functor (const RMV& R, const XMV& X, const aVector& a,
                   const SizeType startingColumn) :
    numCols (X.extent(1)), R_ (R), X_ (X), a_ (a)
  {
    if (startingColumn != 0) {
      auto rng = std::make_pair (startingColumn, static_cast<SizeType> (a.extent(0)));
      a_ = Kokkos::subview (a, rng);
    }
  }

  KOKKOS_INLINE_FUNCTION
  void operator() (const size_type& i) const
  {
    // scalar_x is a compile-time constant (since it is a template
    // parameter), so the compiler should evaluate these branches at
    // compile time.
    if (scalar_x == 0) {
#ifdef KOKKOS_ENABLE_PRAGMA_IVDEP
#pragma ivdep
#endif
#ifdef KOKKOS_ENABLE_PRAGMA_VECTOR
#pragma vector always
#endif
      for (size_type k = 0; k < numCols; ++k) {
        R_(i,k) = ATS::zero ();
      }
    }
    if (scalar_x == -1) {
#ifdef KOKKOS_ENABLE_PRAGMA_IVDEP
#pragma ivdep
#endif
#ifdef KOKKOS_ENABLE_PRAGMA_VECTOR
#pragma vector always
#endif
      for (size_type k = 0; k < numCols; ++k) {
        R_(i,k) = -X_(i,k);
      }
    }
    if (scalar_x == 1) {
#ifdef KOKKOS_ENABLE_PRAGMA_IVDEP
#pragma ivdep
#endif
#ifdef KOKKOS_ENABLE_PRAGMA_VECTOR
#pragma vector always
#endif
      for (size_type k = 0; k < numCols; ++k) {
        R_(i,k) = X_(i,k);
      }
    }
    if (scalar_x == 2) {
#ifdef KOKKOS_ENABLE_PRAGMA_IVDEP
#pragma ivdep
#endif
#ifdef KOKKOS_ENABLE_PRAGMA_VECTOR
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
  typedef SizeType size_type;
  typedef Kokkos::Details::ArithTraits<typename RMV::non_const_value_type> ATS;

  const size_type numCols;
  RMV m_r;
  XMV m_x;
  const typename XMV::non_const_value_type m_a;

  MV_Scal_Functor (const RMV& r, const XMV& x,
                   const typename XMV::non_const_value_type& a,
                   const SizeType /* startingColumn */) :
    numCols (x.extent(1)), m_r (r), m_x (x), m_a (a)
  {}

  KOKKOS_INLINE_FUNCTION
  void operator() (const size_type& i) const
  {
    // scalar_x and scalar_y are compile-time constants (since they
    // are template parameters), so the compiler should evaluate these
    // branches at compile time.
    if (scalar_x == 0) {
#ifdef KOKKOS_ENABLE_PRAGMA_IVDEP
#pragma ivdep
#endif
#ifdef KOKKOS_ENABLE_PRAGMA_VECTOR
#pragma vector always
#endif
      for (size_type k = 0; k < numCols; ++k) {
        m_r(i,k) = ATS::zero ();
      }
    }
    if (scalar_x == -1) {
#ifdef KOKKOS_ENABLE_PRAGMA_IVDEP
#pragma ivdep
#endif
#ifdef KOKKOS_ENABLE_PRAGMA_VECTOR
#pragma vector always
#endif
      for (size_type k = 0; k < numCols; ++k) {
        m_r(i,k) = -m_x(i,k);
      }
    }
    if (scalar_x == 1) {
#ifdef KOKKOS_ENABLE_PRAGMA_IVDEP
#pragma ivdep
#endif
#ifdef KOKKOS_ENABLE_PRAGMA_VECTOR
#pragma vector always
#endif
      for (size_type k = 0; k < numCols; ++k) {
        m_r(i,k) = m_x(i,k);
      }
    }
    if (scalar_x == 2) {
#ifdef KOKKOS_ENABLE_PRAGMA_IVDEP
#pragma ivdep
#endif
#ifdef KOKKOS_ENABLE_PRAGMA_VECTOR
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
  typedef SizeType size_type;
  typedef Kokkos::Details::ArithTraits<typename RMV::non_const_value_type> ATS;

  RMV m_r;
  XMV m_x;
  aVector m_a;

  MV_Scal_Unroll_Functor (const RMV& r, const XMV& x, const aVector& a,
                          const SizeType startingColumn) :
    m_r (r), m_x (x), m_a (a)
  {
    if (startingColumn != 0) {
      auto rng = std::make_pair (startingColumn, static_cast<SizeType> (a.extent(0)));
      m_a = Kokkos::subview (a, rng);
    }
  }

  KOKKOS_INLINE_FUNCTION
  void operator() (const size_type& i) const
  {
    if (scalar_x == 0) {
#ifdef KOKKOS_ENABLE_PRAGMA_UNROLL
#pragma unroll
#endif
      for (int k = 0; k < UNROLL; ++k) {
        m_r(i,k) = ATS::zero ();
      }
    }
    if (scalar_x == -1) {
#ifdef KOKKOS_ENABLE_PRAGMA_UNROLL
#pragma unroll
#endif
      for (int k = 0; k < UNROLL; ++k) {
        m_r(i,k) = -m_x(i,k);
      }
    }
    if (scalar_x == 1) {
#ifdef KOKKOS_ENABLE_PRAGMA_UNROLL
#pragma unroll
#endif
      for (int k = 0; k < UNROLL; ++k) {
        m_r(i,k) = m_x(i,k);
      }
    }
    if (scalar_x == 2) {
#ifdef KOKKOS_ENABLE_PRAGMA_UNROLL
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
  typedef SizeType size_type;
  typedef Kokkos::Details::ArithTraits<typename RMV::non_const_value_type> ATS;

  RMV m_r;
  XMV m_x;
  const typename XMV::non_const_value_type m_a;

  MV_Scal_Unroll_Functor (const RMV& r, const XMV& x,
                          const typename XMV::non_const_value_type& a,
                          const SizeType /* startingColumn */ ) :
    m_r (r), m_x (x), m_a (a)
  {}

  KOKKOS_INLINE_FUNCTION
  void operator() (const size_type& i) const
  {
    if (scalar_x == 0) {
#ifdef KOKKOS_ENABLE_PRAGMA_UNROLL
#pragma unroll
#endif
      for (int k = 0; k < UNROLL; ++k) {
        m_r(i,k) = ATS::zero ();
      }
    }
    if (scalar_x == -1) {
#ifdef KOKKOS_ENABLE_PRAGMA_UNROLL
#pragma unroll
#endif
      for (int k = 0; k < UNROLL; ++k) {
        m_r(i,k) = -m_x(i,k);
      }
    }
    if (scalar_x == 1) {
#ifdef KOKKOS_ENABLE_PRAGMA_UNROLL
#pragma unroll
#endif
      for (int k = 0; k < UNROLL; ++k) {
        m_r(i,k) = m_x(i,k);
      }
    }
    if (scalar_x == 2) {
#ifdef KOKKOS_ENABLE_PRAGMA_UNROLL
#pragma unroll
#endif
      for (int k = 0; k < UNROLL; ++k) {
        m_r(i,k) = m_a*m_x(i,k);
      }
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
MV_Scal_Unrolled (const RMV& r, const aVector& av, const XMV& x,
                  const SizeType startingColumn, int a = 2)
{
  typedef typename XMV::execution_space execution_space;

  if (a == 0) {
    MV_Scal_Unroll_Functor<RMV, aVector, XMV, 0, UNROLL, SizeType> op (r, x, av, startingColumn);
    const SizeType numRows = x.extent(0);
    Kokkos::RangePolicy<execution_space, SizeType> policy (0, numRows);
    Kokkos::parallel_for ("KokkosBlas::Scal::MV::S0", policy, op);
    return;
  }
  if (a == -1) {
    MV_Scal_Unroll_Functor<RMV, aVector, XMV, -1, UNROLL, SizeType> op (r, x, av, startingColumn);
    const SizeType numRows = x.extent(0);
    Kokkos::RangePolicy<execution_space, SizeType> policy (0, numRows);
    Kokkos::parallel_for ("KokkosBlas::Scal::MV::S1", policy, op);
    return;
  }
  if (a == 1) {
    MV_Scal_Unroll_Functor<RMV, aVector, XMV, 1, UNROLL, SizeType> op (r, x, av, startingColumn);
    const SizeType numRows = x.extent(0);
    Kokkos::RangePolicy<execution_space, SizeType> policy (0, numRows);
    Kokkos::parallel_for ("KokkosBlas::Scal::MV::S2", policy, op);
    return;
  }

  // a arbitrary (not -1, 0, or 1)
  MV_Scal_Unroll_Functor<RMV, aVector, XMV, 2, UNROLL, SizeType> op (r, x, av, startingColumn);
  const SizeType numRows = x.extent(0);
  Kokkos::RangePolicy<execution_space, SizeType> policy (0, numRows);
  Kokkos::parallel_for ("KokkosBlas::Scal::MV::S3", policy, op);
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
MV_Scal_Generic (const RVector& r,
                 const aVector& av,
                 const XVector& x,
                 const SizeType startingColumn,
                 int a = 2)
{
  typedef typename XVector::execution_space execution_space;
  const SizeType numRows = x.extent(0);
  Kokkos::RangePolicy<execution_space, SizeType> policy (0, numRows);

  if (a == 0) {
    MV_Scal_Functor<RVector, aVector, XVector, 0, SizeType> op (r, x, av, startingColumn);
    Kokkos::parallel_for ("KokkosBlas::Scal::MV::S4", policy, op);
    return;
  }
  if (a == -1) {
    MV_Scal_Functor<RVector, aVector, XVector, -1, SizeType> op (r, x, av, startingColumn);
    Kokkos::parallel_for ("KokkosBlas::Scal::MV::S5", policy, op);
    return;
  }
  if (a == 1) {
    MV_Scal_Functor<RVector, aVector, XVector, 1, SizeType> op (r, x, av, startingColumn);
    Kokkos::parallel_for ("KokkosBlas::Scal::MV::S6", policy, op);
    return;
  }

  // a arbitrary (not -1, 0, or 1)
  MV_Scal_Functor<RVector, aVector, XVector, 2, SizeType> op (r, x, av, startingColumn);
  Kokkos::parallel_for ("KokkosBlas::Scal::MV::S7", policy, op);
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
  const SizeType numCols = x.extent(1);

#if KOKKOSBLAS_OPTIMIZATION_LEVEL_SCAL <= 2

  // Strip-mine by 8, then 4.  After that, do one column at a time.
  // We limit the number of strip-mine values in order to keep down
  // the amount of code to compile.

  SizeType j = 0; // the current column of X and Y
  for ( ; j + 8 <= numCols; j += 8) {
    const std::pair<SizeType, SizeType> rng (j, j+8);
    auto X_cur = Kokkos::subview (x, Kokkos::ALL (), rng);
    auto R_cur = Kokkos::subview (r, Kokkos::ALL (), rng);
    typedef decltype (X_cur) XMV2D;
    typedef decltype (R_cur) RMV2D;

    MV_Scal_Unrolled<RMV2D, AV, XMV2D, 8, SizeType> (R_cur, av, X_cur, j, a);
  }
  for ( ; j + 4 <= numCols; j += 4) {
    const std::pair<SizeType, SizeType> rng (j, j+4);
    auto X_cur = Kokkos::subview (x, Kokkos::ALL (), rng);
    auto R_cur = Kokkos::subview (r, Kokkos::ALL (), rng);
    typedef decltype (X_cur) XMV2D;
    typedef decltype (R_cur) RMV2D;

    MV_Scal_Unrolled<RMV2D, AV, XMV2D, 4, SizeType> (R_cur, av, X_cur, j, a);
  }
  for ( ; j < numCols; ++j) {
    // RMV and XMV need to turn 1-D.
    auto x_cur = Kokkos::subview (x, Kokkos::ALL (), j);
    auto r_cur = Kokkos::subview (r, Kokkos::ALL (), j);
    typedef decltype (r_cur) RV;
    typedef decltype (x_cur) XV;

    V_Scal_Generic<RV, AV, XV, SizeType> (r_cur, av, x_cur, j, a);
  }

#else // KOKKOSBLAS_OPTIMIZATION_LEVEL_SCAL > 2

  switch (numCols) {
  case 1: {
    auto r_0 = Kokkos::subview (r, Kokkos::ALL (), 0);
    auto x_0 = Kokkos::subview (x, Kokkos::ALL (), 0);
    typedef decltype (r_0) RV;
    typedef decltype (x_0) XV;

    V_Scal_Generic<RV, AV, XV, SizeType> (r_0, av, x_0, 0, a);
    break;
  }
  case 2:
    MV_Scal_Unrolled<RMV, AV, XMV, 2, SizeType> (r, av, x, 0, a);
    break;
  case 3:
    MV_Scal_Unrolled<RMV, AV, XMV, 3, SizeType> (r, av, x, 0, a);
    break;
  case 4:
    MV_Scal_Unrolled<RMV, AV, XMV, 4, SizeType> (r, av, x, 0, a);
    break;
  case 5:
    MV_Scal_Unrolled<RMV, AV, XMV, 5, SizeType> (r, av, x, 0, a);
    break;
  case 6:
    MV_Scal_Unrolled<RMV, AV, XMV, 6, SizeType> (r, av, x, 0, a);
    break;
  case 7:
    MV_Scal_Unrolled<RMV, AV, XMV, 7, SizeType> (r, av, x, 0, a);
    break;
  case 8:
    MV_Scal_Unrolled<RMV, AV, XMV, 8, SizeType> (r, av, x, 0, a);
    break;
  case 9:
    MV_Scal_Unrolled<RMV, AV, XMV, 9, SizeType> (r, av, x, 0, a);
    break;
  case 10:
    MV_Scal_Unrolled<RMV, AV, XMV, 10, SizeType> (r, av, x, 0, a);
    break;
  case 11:
    MV_Scal_Unrolled<RMV, AV, XMV, 11, SizeType> (r, av, x, 0, a);
    break;
  case 12:
    MV_Scal_Unrolled<RMV, AV, XMV, 12, SizeType> (r, av, x, 0, a);
    break;
  case 13:
    MV_Scal_Unrolled<RMV, AV, XMV, 13, SizeType> (r, av, x, 0, a);
    break;
  case 14:
    MV_Scal_Unrolled<RMV, AV, XMV, 14, SizeType> (r, av, x, 0, a);
    break;
  case 15:
    MV_Scal_Unrolled<RMV, AV, XMV, 15, SizeType> (r, av, x, 0, a);
    break;
  case 16:
    MV_Scal_Unrolled<RMV, AV, XMV, 16, SizeType> (r, av, x, 0, a);
    break;
  default:
    MV_Scal_Generic<RMV, AV, XMV, SizeType> (r, av, x, 0, a);
  }

#endif // KOKKOSBLAS_OPTIMIZATION_LEVEL_SCAL
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
  const SizeType numCols = x.extent(1);

  if (numCols == 1) {
    typedef Kokkos::View<typename RMV::value_type*, typename RMV::array_layout,
      typename RMV::device_type, typename RMV::memory_traits> RV;
    typedef Kokkos::View<typename XMV::value_type*, typename XMV::array_layout,
      typename XMV::device_type, typename XMV::memory_traits> XV;

    RV r_0 = Kokkos::subview (r, Kokkos::ALL (), 0);
    XV x_0 = Kokkos::subview (x, Kokkos::ALL (), 0);
    V_Scal_Generic<RMV, aVector, XMV, 1, SizeType> (r_0, av, x_0, a);
  }
  else {
    MV_Scal_Generic<RMV, aVector, XMV, SizeType> (r, av, x, a);
  }
}



} // namespace Impl
} // namespace KokkosBlas

#endif // KOKKOSBLAS1_SCAL_MV_IMPL_HPP_
