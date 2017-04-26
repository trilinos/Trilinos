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
#ifndef KOKKOS_BLAS1_MV_IMPL_SCAL_HPP_
#define KOKKOS_BLAS1_MV_IMPL_SCAL_HPP_

#include <KokkosKernels_config.h>
#include <Kokkos_Core.hpp>
#include <Kokkos_InnerProductSpaceTraits.hpp>

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
    numCols (X.dimension_1 ()), R_ (R), X_ (X), a_ (a)
  {
    if (startingColumn != 0) {
      auto rng = std::make_pair (startingColumn, static_cast<SizeType> (a.dimension_0 ()));
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
  typedef SizeType size_type;
  typedef Kokkos::Details::ArithTraits<typename RMV::non_const_value_type> ATS;

  const size_type numCols;
  RMV m_r;
  XMV m_x;
  const typename XMV::non_const_value_type m_a;

  MV_Scal_Functor (const RMV& r, const XMV& x,
                   const typename XMV::non_const_value_type& a,
                   const SizeType /* startingColumn */) :
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
      auto rng = std::make_pair (startingColumn, static_cast<SizeType> (a.dimension_0 ()));
      m_a = Kokkos::subview (a, rng);
    }
  }

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
  typedef SizeType size_type;
  typedef Kokkos::Details::ArithTraits<typename RV::non_const_value_type> ATS;

  RV m_r;
  XV m_x;
  AV m_a;

  V_Scal_Functor (const RV& r, const XV& x, const AV& a,
                  const SizeType startingColumn) :
    m_r (r), m_x (x), m_a (a)
  {
    static_assert (Kokkos::Impl::is_view<RV>::value,
                   "V_Scal_Functor: RV is not a Kokkos::View.");
    static_assert (Kokkos::Impl::is_view<AV>::value,
                   "V_Scal_Functor: AV is not a Kokkos::View.");
    static_assert (Kokkos::Impl::is_view<XV>::value,
                   "V_Scal_Functor: XV is not a Kokkos::View.");
    static_assert (RV::rank == 1,
                   "V_Scal_Functor: RV is not rank 1.");
    static_assert (AV::rank == 1,
                   "V_Scal_Functor: AV is not rank 1.");
    static_assert (XV::rank == 1,
                   "V_Scal_Functor: XV is not rank 1.");

    if (startingColumn != 0) {
      m_a = Kokkos::subview (a, std::make_pair (startingColumn, static_cast<SizeType> (a.dimension_0 ())));
    }
  }

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
  typedef SizeType size_type;
  typedef Kokkos::Details::ArithTraits<typename RV::non_const_value_type> ATS;

  RV m_r;
  XV m_x;
  const typename XV::non_const_value_type m_a;

  V_Scal_Functor (const RV& r, const XV& x,
                  const typename XV::non_const_value_type& a,
                  const SizeType /* startingColumn */ ) :
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
    if (scalar_x == 1) {
      m_r(i) = m_x(i);
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
MV_Scal_Unrolled (const RMV& r, const aVector& av, const XMV& x,
                  const SizeType startingColumn, int a = 2)
{
  typedef typename XMV::execution_space execution_space;

  if (a == 0) {
    MV_Scal_Unroll_Functor<RMV, aVector, XMV, 0, UNROLL, SizeType> op (r, x, av, startingColumn);
    const SizeType numRows = x.dimension_0 ();
    Kokkos::RangePolicy<execution_space, SizeType> policy (0, numRows);
    Kokkos::parallel_for (policy, op);
    return;
  }
  if (a == -1) {
    MV_Scal_Unroll_Functor<RMV, aVector, XMV, -1, UNROLL, SizeType> op (r, x, av, startingColumn);
    const SizeType numRows = x.dimension_0 ();
    Kokkos::RangePolicy<execution_space, SizeType> policy (0, numRows);
    Kokkos::parallel_for (policy, op);
    return;
  }
  if (a == 1) {
    MV_Scal_Unroll_Functor<RMV, aVector, XMV, 1, UNROLL, SizeType> op (r, x, av, startingColumn);
    const SizeType numRows = x.dimension_0 ();
    Kokkos::RangePolicy<execution_space, SizeType> policy (0, numRows);
    Kokkos::parallel_for (policy, op);
    return;
  }

  // a arbitrary (not -1, 0, or 1)
  MV_Scal_Unroll_Functor<RMV, aVector, XMV, 2, UNROLL, SizeType> op (r, x, av, startingColumn);
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
MV_Scal_Generic (const RVector& r,
                 const aVector& av,
                 const XVector& x,
                 const SizeType startingColumn,
                 int a = 2)
{
  typedef typename XVector::execution_space execution_space;
  const SizeType numRows = x.dimension_0 ();
  Kokkos::RangePolicy<execution_space, SizeType> policy (0, numRows);

  if (a == 0) {
    MV_Scal_Functor<RVector, aVector, XVector, 0, SizeType> op (r, x, av, startingColumn);
    Kokkos::parallel_for (policy, op);
    return;
  }
  if (a == -1) {
    MV_Scal_Functor<RVector, aVector, XVector, -1, SizeType> op (r, x, av, startingColumn);
    Kokkos::parallel_for (policy, op);
    return;
  }
  if (a == 1) {
    MV_Scal_Functor<RVector, aVector, XVector, 1, SizeType> op (r, x, av, startingColumn);
    Kokkos::parallel_for (policy, op);
    return;
  }

  // a arbitrary (not -1, 0, or 1)
  MV_Scal_Functor<RVector, aVector, XVector, 2, SizeType> op (r, x, av, startingColumn);
  Kokkos::parallel_for (policy, op);
}

// Variant of MV_Scal_Generic for single vectors (1-D Views) r and x.
// As above, av is either a 1-D View (and only its first entry will be
// read), or a scalar.
template<class RV, class AV, class XV, class SizeType>
void
V_Scal_Generic (const RV& r, const AV& av, const XV& x,
                const SizeType startingColumn, int a = 2)
{
  static_assert (Kokkos::Impl::is_view<RV>::value,
                 "V_Scal_Generic: RV is not a Kokkos::View.");
  static_assert (Kokkos::Impl::is_view<XV>::value,
                 "V_Scal_Generic: XV is not a Kokkos::View.");
  static_assert (RV::rank == 1,
                 "V_Scal_Generic: RV is not rank 1.");
  static_assert (XV::rank == 1,
                 "V_Scal_Generic: XV is not rank 1.");

  typedef typename RV::execution_space execution_space;
  const SizeType numRows = x.dimension_0 ();
  Kokkos::RangePolicy<execution_space, SizeType> policy (0, numRows);

  if (a == 0) {
    V_Scal_Functor<RV, AV, XV, 0, SizeType> op (r, x, av, startingColumn);
    Kokkos::parallel_for (policy, op);
    return;
  }
  if (a == -1) {
    V_Scal_Functor<RV, AV, XV, -1, SizeType> op (r, x, av, startingColumn);
    Kokkos::parallel_for (policy, op);
    return;
  }
  if (a == 1) {
    V_Scal_Functor<RV, AV, XV, 1, SizeType> op (r, x, av, startingColumn);
    Kokkos::parallel_for (policy, op);
    return;
  }

  // a arbitrary (not -1, 0, or 1)
  V_Scal_Functor<RV, AV, XV, 2, SizeType> op (r, x, av, startingColumn);
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
  const SizeType numCols = x.dimension_1 ();

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

/// \brief Implementation of KokkosBlas::scal for (multi)vectors.
///
/// Compute any of the following:
///
/// 1. R(i,j) = a*X(i,j) for a in -1,0,1
/// 2. R(i,j) = av(j)*X(i,j)
///
/// Only the specializations for rank == 2 or rank == 1 have
/// valid and meaningful definitions.
template<class RMV, class AV, class XMV, int rank = RMV::rank>
struct Scal;

/// \brief Partial specialization of Scal for 2-D Views and 1-D View AV.
///
/// Compute any of the following:
///
/// 1. R(i,j) = a*X(i,j) for a in -1,0,1
/// 2. R(i,j) = alpha(j)*X(i,j)
template<class RMV, class AV, class XMV>
struct Scal<RMV, AV, XMV, 2>
#ifndef KOKKOSKERNELS_ETI_ONLY
{
  typedef typename XMV::size_type size_type;
  typedef Kokkos::Details::ArithTraits<typename XMV::non_const_value_type> ATA;

  static void
  scal (const RMV& R, const AV& av, const XMV& X)
  {
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
}
#endif
;

/// \brief Partial specialization of Scal for 2-D Views and scalar AV.
///
/// Compute any of the following:
///
/// 1. R(i,j) = a*X(i,j) for a in -1,0,1
/// 2. R(i,j) = alpha*X(i,j)
template<class RMV, class XMV>
struct Scal<RMV, typename XMV::non_const_value_type, XMV, 2>
#ifndef KOKKOSKERNELS_ETI_ONLY
{
  typedef typename XMV::non_const_value_type AV;
  typedef typename XMV::size_type size_type;
  typedef Kokkos::Details::ArithTraits<typename XMV::non_const_value_type> ATA;

  static void
  scal (const RMV& R, const AV& alpha, const XMV& X)
  {
    static_assert (Kokkos::Impl::is_view<RMV>::value, "KokkosBlas::Impl::"
                   "Scal<2-D, AV=scalar>: RMV is not a Kokkos::View.");
    static_assert (Kokkos::Impl::is_view<XMV>::value, "KokkosBlas::Impl::"
                   "Scal<2-D, AV=scalar>: XMV is not a Kokkos::View.");
    static_assert (RMV::rank == 2, "KokkosBlas::Impl::Scal<2-D, AV=scalar>: "
                   "RMV is not rank 2.");
    static_assert (XMV::rank == 2, "KokkosBlas::Impl::Scal<2-D, AV=scalar>: "
                   "XMV is not rank 2.");

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
}
#endif
;

/// \brief Partial specialization of Scal for scalar AV (instead of
///   1-D View) and 1-D RMV and XMV.
///
/// Compute any of the following:
///
/// 1. R(i) = a*X(i) for a in -1,0,1
/// 2. R(i) = alpha*X(i)
template<class RMV, class XMV>
struct Scal<RMV, typename RMV::non_const_value_type, XMV, 1>
#ifndef KOKKOSKERNELS_ETI_ONLY
{
  typedef typename XMV::non_const_value_type AV;
  typedef typename XMV::size_type size_type;
  typedef Kokkos::Details::ArithTraits<typename XMV::non_const_value_type> ATA;

  static void
  scal (const RMV& R, const AV& alpha, const XMV& X)
  {
    static_assert (Kokkos::Impl::is_view<RMV>::value, "KokkosBlas::Impl::"
                   "Scal<1-D>: RMV is not a Kokkos::View.");
    static_assert (Kokkos::Impl::is_view<XMV>::value, "KokkosBlas::Impl::"
                   "Scal<1-D>: XMV is not a Kokkos::View.");
    static_assert (RMV::rank == 1, "KokkosBlas::Impl::Scal<1-D>: "
                   "RMV is not rank 1.");
    static_assert (XMV::rank == 1, "KokkosBlas::Impl::Scal<1-D>: "
                   "XMV is not rank 1.");

    const size_type numRows = X.dimension_0 ();
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

    if (numRows < static_cast<size_type> (INT_MAX)) {
      typedef int index_type;
      V_Scal_Generic<RMV, AV, XMV, index_type> (R, alpha, X, a);
    }
    else {
      typedef typename XMV::size_type index_type;
      V_Scal_Generic<RMV, AV, XMV, index_type> (R, alpha, X, a);
    }
  }
}
#endif
;

//
// Macro for declaration of full specialization of
// KokkosBlas::Impl::Scal, for rank == 2 and a 1-D View of
// coefficients.  This is NOT for users!!!  All the declarations of
// full specializations go in this header file.  We may spread out
// definitions (see _DEF macro below) across one or more .cpp files.
//

#define KOKKOSBLAS1_IMPL_MV_SCAL_MULTICOEFF_DECL( SCALAR, LAYOUT, EXEC_SPACE, MEM_SPACE ) \
extern template struct Scal<Kokkos::View<SCALAR**, \
                         LAYOUT, \
                         Kokkos::Device<EXEC_SPACE, MEM_SPACE>, \
                         Kokkos::MemoryTraits<Kokkos::Unmanaged> >, \
            Kokkos::View<const SCALAR*, \
                         EXEC_SPACE::array_layout, \
                         Kokkos::Device<EXEC_SPACE, MEM_SPACE>, \
                         Kokkos::MemoryTraits<Kokkos::Unmanaged> >, \
            Kokkos::View<const SCALAR**, \
                         LAYOUT, \
                         Kokkos::Device<EXEC_SPACE, MEM_SPACE>, \
                         Kokkos::MemoryTraits<Kokkos::Unmanaged> >, \
                            2>;

//
// Macro for declaration of full specialization of
// KokkosBlas::Impl::Scal, for rank == 2 and a single scalar
// coefficient.  This is NOT for users!!!  All the declarations of
// full specializations go in this header file.  We may spread out
// definitions (see _DEF macro below) across one or more .cpp files.
//

#define KOKKOSBLAS1_IMPL_MV_SCAL_SINGLECOEFF_DECL( SCALAR, LAYOUT, EXEC_SPACE, MEM_SPACE ) \
extern template struct Scal<Kokkos::View<SCALAR**, \
                         LAYOUT, \
                         Kokkos::Device<EXEC_SPACE, MEM_SPACE>, \
                         Kokkos::MemoryTraits<Kokkos::Unmanaged> >, \
            SCALAR, \
            Kokkos::View<const SCALAR**, \
                         LAYOUT, \
                         Kokkos::Device<EXEC_SPACE, MEM_SPACE>, \
                         Kokkos::MemoryTraits<Kokkos::Unmanaged> >, \
                            2>;

//
// Macro for declaration of full specialization of
// KokkosBlas::Impl::Scal, for rank == 1 and a single scalar
// coefficient.  This is NOT for users!!!  All the declarations of
// full specializations go in this header file.  We may spread out
// definitions (see _DEF macro below) across one or more .cpp files.
//

#define KOKKOSBLAS1_IMPL_V_SCAL_SINGLECOEFF_DECL( SCALAR, LAYOUT, EXEC_SPACE, MEM_SPACE ) \
extern template struct Scal<Kokkos::View<SCALAR*, \
                         LAYOUT, \
                         Kokkos::Device<EXEC_SPACE, MEM_SPACE>, \
                         Kokkos::MemoryTraits<Kokkos::Unmanaged> >, \
            SCALAR, \
            Kokkos::View<const SCALAR*, \
                         LAYOUT, \
                         Kokkos::Device<EXEC_SPACE, MEM_SPACE>, \
                         Kokkos::MemoryTraits<Kokkos::Unmanaged> >, \
                            1>;

//
// Macro for definition of full specialization of
// KokkosBlas::Impl::Scalar for rank == 2 and a single scalar
// coefficient.  This is NOT for users!!!
//

#define KOKKOSBLAS1_IMPL_MV_SCAL_SINGLECOEFF_DEF( SCALAR, LAYOUT, EXEC_SPACE, MEM_SPACE ) \
template struct Scal<Kokkos::View<SCALAR**, \
                  LAYOUT, \
                  Kokkos::Device<EXEC_SPACE, MEM_SPACE>, \
                  Kokkos::MemoryTraits<Kokkos::Unmanaged> >, \
     SCALAR, \
     Kokkos::View<const SCALAR**, \
                  LAYOUT, \
                  Kokkos::Device<EXEC_SPACE, MEM_SPACE>, \
                  Kokkos::MemoryTraits<Kokkos::Unmanaged> >, \
     2>;

//
// Macro for definition of full specialization of
// KokkosBlas::Impl::Scalar for rank == 2 and a 1-D View of
// coefficients.  This is NOT for users!!!
//

#define KOKKOSBLAS1_IMPL_MV_SCAL_MULTICOEFF_DEF( SCALAR, LAYOUT, EXEC_SPACE, MEM_SPACE ) \
template struct Scal<Kokkos::View<SCALAR**, \
                  LAYOUT, \
                  Kokkos::Device<EXEC_SPACE, MEM_SPACE>, \
                  Kokkos::MemoryTraits<Kokkos::Unmanaged> >, \
     Kokkos::View<const SCALAR*,                         \
                  EXEC_SPACE::array_layout,                     \
                  Kokkos::Device<EXEC_SPACE, MEM_SPACE>,          \
                  Kokkos::MemoryTraits<Kokkos::Unmanaged> >, \
     Kokkos::View<const SCALAR**, \
                  LAYOUT, \
                  Kokkos::Device<EXEC_SPACE, MEM_SPACE>, \
                  Kokkos::MemoryTraits<Kokkos::Unmanaged> >, \
                     2>;

//
// Macro for definition of full specialization of
// KokkosBlas::Impl::Scalar for rank == 1 and a single scalar
// coefficient.  This is NOT for users!!!
//

#define KOKKOSBLAS1_IMPL_V_SCAL_SINGLECOEFF_DEF( SCALAR, LAYOUT, EXEC_SPACE, MEM_SPACE ) \
template struct Scal<Kokkos::View<SCALAR*, \
                  LAYOUT, \
                  Kokkos::Device<EXEC_SPACE, MEM_SPACE>, \
                  Kokkos::MemoryTraits<Kokkos::Unmanaged> >, \
     SCALAR, \
     Kokkos::View<const SCALAR*, \
                  LAYOUT, \
                  Kokkos::Device<EXEC_SPACE, MEM_SPACE>, \
                  Kokkos::MemoryTraits<Kokkos::Unmanaged> >, \
                     1>;

} // namespace Impl
} // namespace KokkosBlas

#include<generated_specializations_hpp/KokkosBlas1_impl_MV_scal_multicoeff_decl_specializations.hpp>
#include<generated_specializations_hpp/KokkosBlas1_impl_MV_scal_singlecoeff_decl_specializations.hpp>
#include<generated_specializations_hpp/KokkosBlas1_impl_V_scal_singlecoeff_decl_specializations.hpp>

#endif // KOKKOS_BLAS1_MV_IMPL_SCAL_HPP_
