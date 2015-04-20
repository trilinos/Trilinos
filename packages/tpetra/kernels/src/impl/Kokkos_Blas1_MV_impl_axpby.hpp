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
#ifndef KOKKOS_BLAS1_MV_IMPL_AXPBY_HPP_
#define KOKKOS_BLAS1_MV_IMPL_AXPBY_HPP_

#include <TpetraKernels_config.h>
#include <Kokkos_Core.hpp>
#include <Kokkos_InnerProductSpaceTraits.hpp>

namespace KokkosBlas {
namespace Impl {

//
// axpby
//

// Functor for multivectors X and Y and 1-D views a and b, that
// computes any of the following:
//
// 1. Y(i,j) = alpha*X(i,j) + beta*Y(i,j) for alpha,beta in -1,0,1
// 2. Y(i,j) = a(j)*X(i,j) + beta*Y(i,j) for beta in -1,0,1
// 3. Y(i,j) = alpha*X(i,j) + beta*Y(i,j) for alpha in -1,0,1
// 4. Y(i,j) = a(j)*X(i,j) + b(j)*Y(i,j)
//
// The template parameters scalar_x and scalar_y correspond to alpha
// resp. beta in the operation y = alpha*x + beta*y.  The values -1,
// 0, and -1 correspond to literal values of those coefficients.  The
// value 2 tells the functor to use the corresponding vector of
// coefficients.  Any literal coefficient of zero has BLAS semantics
// of ignoring the corresponding (multi)vector entry.  This does not
// apply to coefficients in the a and b vectors, if they are used.
template<class AV, class XMV, class BV, class YMV,
         int scalar_x, int scalar_y, class SizeType = typename YMV::size_type>
struct MV_Axpby_Functor
{
  typedef typename YMV::execution_space execution_space;
  typedef SizeType size_type;
  typedef Kokkos::Details::ArithTraits<typename YMV::non_const_value_type> ATS;

  const size_type numCols;
  XMV m_x;
  YMV m_y;
  AV m_a;
  BV m_b;

  MV_Axpby_Functor (const XMV& X, const YMV& Y, const AV& a, const BV& b) :
    numCols (X.dimension_1 ()), m_x (X), m_y (Y), m_a (a), m_b (b)
  {
#ifdef KOKKOS_HAVE_CXX11
    // XMV and YMV must be Kokkos::View specializations.
    static_assert (Kokkos::Impl::is_view<AV>::value, "KokkosBlas::Impl::"
                   "MV_Axpby_Functor: a is not a Kokkos::View.");
    static_assert (Kokkos::Impl::is_view<XMV>::value, "KokkosBlas::Impl::"
                   "MV_Axpby_Functor: X is not a Kokkos::View.");
    static_assert (Kokkos::Impl::is_view<BV>::value, "KokkosBlas::Impl::"
                   "MV_Axpby_Functor: b is not a Kokkos::View.");
    static_assert (Kokkos::Impl::is_view<YMV>::value, "KokkosBlas::Impl::"
                   "MV_Axpby_Functor: Y is not a Kokkos::View.");
    // YMV must be nonconst (else it can't be an output argument).
    static_assert (Kokkos::Impl::is_same<typename YMV::value_type,
                   typename YMV::non_const_value_type>::value,
                   "KokkosBlas::Impl::MV_Axpby_Functor: Y is const.  "
                   "It must be nonconst, because it is an output argument "
                   "(we have to be able to write to its entries).");
    static_assert ((int) YMV::rank == (int) XMV::rank, "KokkosBlas::Impl::MV_Axpby_Functor: "
                   "X and Y must have the same rank.");
    static_assert (YMV::rank == 2, "KokkosBlas::Impl::MV_Axpby_Functor: "
                   "XMV and YMV must have rank 2.");
    static_assert (AV::rank == 1, "KokkosBlas::Impl::MV_Axpby_Functor: "
                   "AV must have rank 1.");
    static_assert (BV::rank == 1, "KokkosBlas::Impl::MV_Axpby_Functor: "
                   "BV must have rank 1.");
#endif // KOKKOS_HAVE_CXX11
  }

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
        m_y(i,k) = ATS::zero ();
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
        m_y(i,k) = -m_y(i,k);
      }
    }
    if (scalar_x == 0 && scalar_y == 1) {
      return; // Y(i,j) := Y(i,j)
    }
    if (scalar_x == 0 && scalar_y == 2) {
#ifdef KOKKOS_HAVE_PRAGMA_IVDEP
#pragma ivdep
#endif
#ifdef KOKKOS_HAVE_PRAGMA_VECTOR
#pragma vector always
#endif
      for (size_type k = 0; k < numCols; ++k) {
        m_y(i,k) = m_b(k)*m_y(i,k);
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
        m_y(i,k) = -m_x(i,k);
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
        m_y(i,k) = -m_x(i,k) - m_y(i,k);
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
        m_y(i,k) = -m_x(i,k) + m_y(i,k);
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
        m_y(i,k) = -m_x(i,k) + m_b(k)*m_y(i,k);
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
        m_y(i,k) = m_x(i,k);
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
        m_y(i,k) = m_x(i,k) - m_y(i,k);
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
        m_y(i,k) = m_x(i,k) + m_y(i,k);
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
        m_y(i,k) = m_x(i,k) + m_b(k)*m_y(i,k);
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
        m_y(i,k) = m_a(k)*m_x(i,k);
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
        m_y(i,k) = m_a(k)*m_x(i,k) - m_y(i,k);
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
        m_y(i,k) = m_a(k)*m_x(i,k) + m_y(i,k);
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
        m_y(i,k) = m_a(k)*m_x(i,k) + m_b(k)*m_y(i,k);
      }
    }
  }
};

// Variant of MV_Axpby_Functor, where a and b are scalars.
// This functor computes any of the following:
//
// 1. Y(i,j) = alpha*X(i,j) + beta*Y(i,j) for alpha,beta in -1,0,1
// 2. Y(i,j) = a*X(i,j) + beta*Y(i,j) for beta in -1,0,1
// 3. Y(i,j) = alpha*X(i,j) + beta*Y(i,j) for alpha in -1,0,1
// 4. Y(i,j) = a*X(i,j) + b*Y(i,j)
//
// The template parameters scalar_x and scalar_y correspond to alpha
// resp. beta in the operation y = alpha*x + beta*y.  The values -1,
// 0, and -1 correspond to literal values of those coefficients.  The
// value 2 tells the functor to use the corresponding vector of
// coefficients.  Any literal coefficient of zero has BLAS semantics
// of ignoring the corresponding (multi)vector entry.  This does not
// apply to coefficients in the a and b vectors, if they are used.
//
// This version works by partial specialization on AV and BV.
// In this partial specialization, both AV and BV are scalars.
template<class XMV, class YMV, int scalar_x, int scalar_y, class SizeType>
struct MV_Axpby_Functor<typename XMV::non_const_value_type, XMV,
                        typename YMV::non_const_value_type, YMV,
                        scalar_x, scalar_y, SizeType>
{
  typedef typename YMV::execution_space execution_space;
  typedef SizeType size_type;
  typedef Kokkos::Details::ArithTraits<typename YMV::non_const_value_type> ATS;

  const size_type numCols;
  XMV m_x;
  YMV m_y;
  const typename XMV::non_const_value_type m_a;
  const typename YMV::non_const_value_type m_b;

  MV_Axpby_Functor (const XMV& X, const YMV& Y,
                    const typename XMV::non_const_value_type& a,
                    const typename YMV::non_const_value_type& b) :
    numCols (X.dimension_1 ()), m_x (X), m_y (Y), m_a (a), m_b (b)
  {
#ifdef KOKKOS_HAVE_CXX11
    static_assert (Kokkos::Impl::is_view<XMV>::value, "KokkosBlas::Impl::"
                   "MV_Axpby_Functor: X is not a Kokkos::View.");
    static_assert (Kokkos::Impl::is_view<YMV>::value, "KokkosBlas::Impl::"
                   "MV_Axpby_Functor: Y is not a Kokkos::View.");
    // YMV must be nonconst (else it can't be an output argument).
    static_assert (Kokkos::Impl::is_same<typename YMV::value_type,
                   typename YMV::non_const_value_type>::value,
                   "KokkosBlas::Impl::MV_Axpby_Functor: Y is const.  "
                   "It must be nonconst, because it is an output argument "
                   "(we have to be able to write to its entries).");
    static_assert ((int) YMV::rank == (int) XMV::rank, "KokkosBlas::Impl::"
                   "MV_Axpby_Functor: X and Y must have the same rank.");
    static_assert (YMV::rank == 2, "KokkosBlas::Impl::MV_Axpby_Functor: "
                   "XMV and YMV must have rank 2.");
#endif // KOKKOS_HAVE_CXX11
  }

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
        m_y(i,k) = ATS::zero ();
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
        m_y(i,k) = -m_y(i,k);
      }
    }
    if (scalar_x == 0 && scalar_y == 1) {
      return; // Y(i,j) := Y(i,j)
    }
    if (scalar_x == 0 && scalar_y == 2) {
#ifdef KOKKOS_HAVE_PRAGMA_IVDEP
#pragma ivdep
#endif
#ifdef KOKKOS_HAVE_PRAGMA_VECTOR
#pragma vector always
#endif
      for (size_type k = 0; k < numCols; ++k) {
        m_y(i,k) = m_b*m_y(i,k);
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
        m_y(i,k) = -m_x(i,k);
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
        m_y(i,k) = -m_x(i,k) - m_y(i,k);
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
        m_y(i,k) = -m_x(i,k) + m_y(i,k);
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
        m_y(i,k) = -m_x(i,k) + m_b*m_y(i,k);
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
        m_y(i,k) = m_x(i,k);
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
        m_y(i,k) = m_x(i,k) - m_y(i,k);
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
        m_y(i,k) = m_x(i,k) + m_y(i,k);
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
        m_y(i,k) = m_x(i,k) + m_b*m_y(i,k);
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
        m_y(i,k) = m_a*m_x(i,k);
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
        m_y(i,k) = m_a*m_x(i,k) - m_y(i,k);
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
        m_y(i,k) = m_a*m_x(i,k) + m_y(i,k);
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
        m_y(i,k) = m_a*m_x(i,k) + m_b*m_y(i,k);
      }
    }
  }
};


// Column-unrolled variant of MV_Axpby_Functor.  The number of columns
// in X and Y, UNROLL, is a compile-time constant.
template<class AV, class XMV, class BV, class YMV,
         int scalar_x, int scalar_y, int UNROLL, class SizeType>
struct MV_Axpby_Unroll_Functor
{
  typedef typename YMV::execution_space execution_space;
  typedef SizeType size_type;
  typedef Kokkos::Details::ArithTraits<typename YMV::non_const_value_type> ATS;

  XMV m_x;
  YMV m_y;
  AV m_a;
  BV m_b;

  MV_Axpby_Unroll_Functor (const XMV& x, const YMV& y, const AV& a, const BV& b) :
    m_x (x), m_y (y), m_a (a), m_b (b)
  {
#ifdef KOKKOS_HAVE_CXX11
    static_assert (Kokkos::Impl::is_view<AV>::value, "KokkosBlas::Impl::"
                   "MV_Axpby_Unroll_Functor: a is not a Kokkos::View.");
    static_assert (Kokkos::Impl::is_view<XMV>::value, "KokkosBlas::Impl::"
                   "MV_Axpby_Unroll_Functor: X is not a Kokkos::View.");
    static_assert (Kokkos::Impl::is_view<BV>::value, "KokkosBlas::Impl::"
                   "MV_Axpby_Unroll_Functor: b is not a Kokkos::View.");
    static_assert (Kokkos::Impl::is_view<YMV>::value, "KokkosBlas::Impl::"
                   "MV_Axpby_Unroll_Functor: Y is not a Kokkos::View.");
    static_assert (Kokkos::Impl::is_same<typename YMV::value_type,
                   typename YMV::non_const_value_type>::value,
                   "KokkosBlas::Impl::MV_Axpby_Unroll_Functor: Y is const.  "
                   "It must be nonconst, because it is an output argument "
                   "(we have to be able to write to its entries).");
    static_assert ((int) YMV::rank == (int) XMV::rank,
                   "KokkosBlas::Impl::MV_Axpby_Unroll_Functor: "
                   "X and Y must have the same rank.");
    static_assert (YMV::rank == 2, "KokkosBlas::Impl::MV_Axpby_Unroll_Functor: "
                   "XMV and YMV must have rank 2.");
    static_assert (AV::rank == 1, "KokkosBlas::Impl::MV_Axpby_Unroll_Functor: "
                   "AV must have rank 1.");
    static_assert (BV::rank == 1, "KokkosBlas::Impl::MV_Axpby_Unroll_Functor: "
                   "BV must have rank 1.");
#endif // KOKKOS_HAVE_CXX11
  }

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
        m_y(i,k) = ATS::zero ();
      }
    }
    if (scalar_x == 0 && scalar_y == -1) {
#ifdef KOKKOS_HAVE_PRAGMA_UNROLL
#pragma unroll
#endif
      for (int k = 0; k < UNROLL; ++k) {
        m_y(i,k) = -m_y(i,k);
      }
    }
    if (scalar_x == 0 && scalar_y == 1) {
      return; // Y(i,j) := Y(i,j)
    }
    if (scalar_x == 0 && scalar_y == 2) {
#ifdef KOKKOS_HAVE_PRAGMA_UNROLL
#pragma unroll
#endif
      for (int k = 0; k < UNROLL; ++k) {
        m_y(i,k) = m_b(k)*m_y(i,k);
      }
    }
    if (scalar_x == -1 && scalar_y == 0) {
#ifdef KOKKOS_HAVE_PRAGMA_UNROLL
#pragma unroll
#endif
      for (int k = 0; k < UNROLL; ++k) {
        m_y(i,k) = -m_x(i,k);
      }
    }
    if (scalar_x == -1 && scalar_y == -1) {
#ifdef KOKKOS_HAVE_PRAGMA_UNROLL
#pragma unroll
#endif
      for (int k = 0; k < UNROLL; ++k) {
        m_y(i,k) = -m_x(i,k) - m_y(i,k);
      }
    }
    if (scalar_x == -1 && scalar_y == 1) {
#ifdef KOKKOS_HAVE_PRAGMA_UNROLL
#pragma unroll
#endif
      for (int k = 0; k < UNROLL; ++k) {
        m_y(i,k) = -m_x(i,k) + m_y(i,k);
      }
    }
    if (scalar_x == -1 && scalar_y == 2) {
#ifdef KOKKOS_HAVE_PRAGMA_UNROLL
#pragma unroll
#endif
      for (int k = 0; k < UNROLL; ++k) {
        m_y(i,k) = -m_x(i,k) + m_b(k)*m_y(i,k);
      }
    }
    if (scalar_x == 1 && scalar_y == 0) {
#ifdef KOKKOS_HAVE_PRAGMA_UNROLL
#pragma unroll
#endif
      for (int k = 0; k < UNROLL; ++k) {
        m_y(i,k) = m_x(i,k);
      }
    }
    if (scalar_x == 1 && scalar_y == -1) {
#ifdef KOKKOS_HAVE_PRAGMA_UNROLL
#pragma unroll
#endif
      for (int k = 0; k < UNROLL; ++k) {
        m_y(i,k) = m_x(i,k) - m_y(i,k);
      }
    }
    if (scalar_x == 1 && scalar_y == 1) {
#ifdef KOKKOS_HAVE_PRAGMA_UNROLL
#pragma unroll
#endif
      for (int k = 0; k < UNROLL; ++k) {
        m_y(i,k) = m_x(i,k) + m_y(i,k);
      }
    }
    if (scalar_x == 1 && scalar_y == 2) {
#ifdef KOKKOS_HAVE_PRAGMA_UNROLL
#pragma unroll
#endif
      for (int k = 0; k < UNROLL; ++k) {
        m_y(i,k) = m_x(i,k) + m_b(k)*m_y(i,k);
      }
    }
    if (scalar_x == 2 && scalar_y == 0) {
#ifdef KOKKOS_HAVE_PRAGMA_UNROLL
#pragma unroll
#endif
      for (int k = 0; k < UNROLL; ++k) {
        m_y(i,k) = m_a(k)*m_x(i,k);
      }
    }
    if (scalar_x == 2 && scalar_y == -1) {
#ifdef KOKKOS_HAVE_PRAGMA_UNROLL
#pragma unroll
#endif
      for (int k = 0; k < UNROLL; ++k) {
        m_y(i,k) = m_a(k)*m_x(i,k) - m_y(i,k);
      }
    }
    if (scalar_x == 2 && scalar_y == 1) {
#ifdef KOKKOS_HAVE_PRAGMA_UNROLL
#pragma unroll
#endif
      for (int k = 0; k < UNROLL; ++k) {
        m_y(i,k) = m_a(k)*m_x(i,k) + m_y(i,k);
      }
    }
    if (scalar_x == 2 && scalar_y == 2) {
#ifdef KOKKOS_HAVE_PRAGMA_UNROLL
#pragma unroll
#endif
      for (int k = 0; k < UNROLL; ++k) {
        m_y(i,k) = m_a(k)*m_x(i,k) + m_b(k)*m_y(i,k);
      }
    }
  }
};


// Variant of MV_Axpby_Unroll_Functor for single coefficients (rather
// than vectors of coefficients) a and b.  The number of columns in X
// and Y, UNROLL, is a compile-time constant.
template<class XMV, class YMV, int scalar_x, int scalar_y,
         int UNROLL, class SizeType>
struct MV_Axpby_Unroll_Functor<typename XMV::non_const_value_type, XMV,
                               typename YMV::non_const_value_type, YMV,
                               scalar_x, scalar_y, UNROLL, SizeType>
{
  typedef typename YMV::execution_space execution_space;
  typedef SizeType size_type;
  typedef Kokkos::Details::ArithTraits<typename YMV::non_const_value_type> ATS;

  XMV m_x;
  YMV m_y;
  const typename XMV::non_const_value_type m_a;
  const typename YMV::non_const_value_type m_b;

  MV_Axpby_Unroll_Functor (const XMV& X, const YMV& Y,
                           const typename XMV::non_const_value_type& a,
                           const typename YMV::non_const_value_type& b) :
    m_x (X), m_y (Y), m_a (a), m_b (b)
  {
#ifdef KOKKOS_HAVE_CXX11
    static_assert (Kokkos::Impl::is_view<XMV>::value, "KokkosBlas::Impl::"
                   "MV_Axpby_Unroll_Functor: X is not a Kokkos::View.");
    static_assert (Kokkos::Impl::is_view<YMV>::value, "KokkosBlas::Impl::"
                   "MV_Axpby_Unroll_Functor: Y is not a Kokkos::View.");
    static_assert (Kokkos::Impl::is_same<typename YMV::value_type,
                   typename YMV::non_const_value_type>::value,
                   "KokkosBlas::Impl::MV_Axpby_Unroll_Functor: Y is const.  "
                   "It must be nonconst, because it is an output argument "
                   "(we have to be able to write to its entries).");
    static_assert ((int) YMV::rank == (int) XMV::rank, "KokkosBlas::Impl::MV_Axpby_Unroll_Functor: "
                   "X and Y must have the same rank.");
    static_assert (YMV::rank == 2, "KokkosBlas::Impl::MV_Axpby_Unroll_Functor: "
                   "XMV and YMV must have rank 2.");
#endif // KOKKOS_HAVE_CXX11
  }

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
        m_y(i,k) = ATS::zero ();
      }
    }
    if (scalar_x == 0 && scalar_y == -1) {
#ifdef KOKKOS_HAVE_PRAGMA_UNROLL
#pragma unroll
#endif
      for (int k = 0; k < UNROLL; ++k) {
        m_y(i,k) = -m_y(i,k);
      }
    }
    if (scalar_x == 0 && scalar_y == 1) {
      return; // Y(i,j) := Y(i,j)
    }
    if (scalar_x == 0 && scalar_y == 2) {
#ifdef KOKKOS_HAVE_PRAGMA_UNROLL
#pragma unroll
#endif
      for (int k = 0; k < UNROLL; ++k) {
        m_y(i,k) = m_b*m_y(i,k);
      }
    }
    if (scalar_x == -1 && scalar_y == 0) {
#ifdef KOKKOS_HAVE_PRAGMA_UNROLL
#pragma unroll
#endif
      for (int k = 0; k < UNROLL; ++k) {
        m_y(i,k) = -m_x(i,k);
      }
    }
    if (scalar_x == -1 && scalar_y == -1) {
#ifdef KOKKOS_HAVE_PRAGMA_UNROLL
#pragma unroll
#endif
      for (int k = 0; k < UNROLL; ++k) {
        m_y(i,k) = -m_x(i,k) - m_y(i,k);
      }
    }
    if (scalar_x == -1 && scalar_y == 1) {
#ifdef KOKKOS_HAVE_PRAGMA_UNROLL
#pragma unroll
#endif
      for (int k = 0; k < UNROLL; ++k) {
        m_y(i,k) = -m_x(i,k) + m_y(i,k);
      }
    }
    if (scalar_x == -1 && scalar_y == 2) {
#ifdef KOKKOS_HAVE_PRAGMA_UNROLL
#pragma unroll
#endif
      for (int k = 0; k < UNROLL; ++k) {
        m_y(i,k) = -m_x(i,k) + m_b*m_y(i,k);
      }
    }
    if (scalar_x == 1 && scalar_y == 0) {
#ifdef KOKKOS_HAVE_PRAGMA_UNROLL
#pragma unroll
#endif
      for (int k = 0; k < UNROLL; ++k) {
        m_y(i,k) = m_x(i,k);
      }
    }
    if (scalar_x == 1 && scalar_y == -1) {
#ifdef KOKKOS_HAVE_PRAGMA_UNROLL
#pragma unroll
#endif
      for (int k = 0; k < UNROLL; ++k) {
        m_y(i,k) = m_x(i,k) - m_y(i,k);
      }
    }
    if (scalar_x == 1 && scalar_y == 1) {
#ifdef KOKKOS_HAVE_PRAGMA_UNROLL
#pragma unroll
#endif
      for (int k = 0; k < UNROLL; ++k) {
        m_y(i,k) = m_x(i,k) + m_y(i,k);
      }
    }
    if (scalar_x == 1 && scalar_y == 2) {
#ifdef KOKKOS_HAVE_PRAGMA_UNROLL
#pragma unroll
#endif
      for (int k = 0; k < UNROLL; ++k) {
        m_y(i,k) = m_x(i,k) + m_b*m_y(i,k);
      }
    }
    if (scalar_x == 2 && scalar_y == 0) {
#ifdef KOKKOS_HAVE_PRAGMA_UNROLL
#pragma unroll
#endif
      for (int k = 0; k < UNROLL; ++k) {
        m_y(i,k) = m_a*m_x(i,k);
      }
    }
    if (scalar_x == 2 && scalar_y == -1) {
#ifdef KOKKOS_HAVE_PRAGMA_UNROLL
#pragma unroll
#endif
      for (int k = 0; k < UNROLL; ++k) {
        m_y(i,k) = m_a*m_x(i,k) - m_y(i,k);
      }
    }
    if (scalar_x == 2 && scalar_y == 1) {
#ifdef KOKKOS_HAVE_PRAGMA_UNROLL
#pragma unroll
#endif
      for (int k = 0; k < UNROLL; ++k) {
        m_y(i,k) = m_a*m_x(i,k) + m_y(i,k);
      }
    }
    if (scalar_x == 2 && scalar_y == 2) {
#ifdef KOKKOS_HAVE_PRAGMA_UNROLL
#pragma unroll
#endif
      for (int k = 0; k < UNROLL; ++k) {
        m_y(i,k) = m_a*m_x(i,k) + m_b*m_y(i,k);
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
template<class AV, class XV, class BV, class YV,
         int scalar_x, int scalar_y, class SizeType>
struct V_Axpby_Functor {
  typedef typename YV::execution_space execution_space;
  typedef SizeType size_type;
  typedef Kokkos::Details::ArithTraits<typename YV::non_const_value_type> ATS;

  XV m_x;
  YV m_y;
  AV m_a;
  BV m_b;

  V_Axpby_Functor (const XV& x, const YV& y, const AV& a, const BV& b) :
    m_x (x), m_y (y), m_a (a), m_b (b)
  {
#ifdef KOKKOS_HAVE_CXX11
    static_assert (Kokkos::Impl::is_view<XV>::value, "KokkosBlas::Impl::"
                   "V_Axpby_Functor: X is not a Kokkos::View.");
    static_assert (Kokkos::Impl::is_view<YV>::value, "KokkosBlas::Impl::"
                   "V_Axpby_Functor: Y is not a Kokkos::View.");
    static_assert (Kokkos::Impl::is_same<typename YV::value_type,
                   typename YV::non_const_value_type>::value,
                   "KokkosBlas::Impl::V_Axpby_Functor: Y is const.  "
                   "It must be nonconst, because it is an output argument "
                   "(we have to be able to write to its entries).");
    static_assert ((int) YV::rank == (int) XV::rank, "KokkosBlas::Impl::V_Axpby_Functor: "
                   "X and Y must have the same rank.");
    static_assert (YV::rank == 1, "KokkosBlas::Impl::V_Axpby_Functor: "
                   "XV and YV must have rank 1.");
#endif // KOKKOS_HAVE_CXX11
  }

  KOKKOS_INLINE_FUNCTION
  void operator() (const size_type& i) const
  {
    // scalar_x and scalar_y are compile-time constants (since they
    // are template parameters), so the compiler should evaluate these
    // branches at compile time.
    if (scalar_x == 0 && scalar_y == 0) {
      m_y(i) = ATS::zero ();
    }
    if (scalar_x == 0 && scalar_y == -1) {
      m_y(i) = -m_y(i);
    }
    if (scalar_x == 0 && scalar_y == 1) {
      return; // m_y(i) = m_y(i);
    }
    if (scalar_x == 0 && scalar_y == 2) {
      m_y(i) = m_b(0)*m_y(i);
    }
    if (scalar_x == -1 && scalar_y == 0) {
      m_y(i) = -m_x(i);
    }
    if (scalar_x == -1 && scalar_y == -1) {
      m_y(i) = -m_x(i) - m_y(i);
    }
    if (scalar_x == -1 && scalar_y == 1) {
      m_y(i) = -m_x(i) + m_y(i);
    }
    if (scalar_x == -1 && scalar_y == 2) {
      m_y(i) = -m_x(i) + m_b(0)*m_y(i);
    }
    if (scalar_x == 1 && scalar_y == 0) {
      m_y(i) = m_x(i);
    }
    if (scalar_x == 1 && scalar_y == -1) {
      m_y(i) = m_x(i) - m_y(i);
    }
    if (scalar_x == 1 && scalar_y == 1) {
      m_y(i) = m_x(i) + m_y(i);
    }
    if (scalar_x == 1 && scalar_y == 2) {
      m_y(i) = m_x(i) + m_b(0)*m_y(i);
    }
    if (scalar_x == 2 && scalar_y == 0) {
      m_y(i) = m_a(0)*m_x(i);
    }
    if (scalar_x == 2 && scalar_y == -1) {
      m_y(i) = m_a(0)*m_x(i) - m_y(i);
    }
    if (scalar_x == 2 && scalar_y == 1) {
      m_y(i) = m_a(0)*m_x(i) + m_y(i);
    }
    if (scalar_x == 2 && scalar_y == 2) {
      m_y(i) = m_a(0)*m_x(i) + m_b(0)*m_y(i);
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
template<class XV, class YV,
         int scalar_x, int scalar_y, class SizeType>
struct V_Axpby_Functor<typename XV::non_const_value_type, XV,
                       typename YV::non_const_value_type, YV,
                       scalar_x, scalar_y, SizeType> {
  typedef typename YV::execution_space execution_space;
  typedef SizeType size_type;
  typedef Kokkos::Details::ArithTraits<typename YV::non_const_value_type> ATS;

  XV m_x;
  YV m_y;
  const typename XV::non_const_value_type m_a;
  const typename YV::non_const_value_type m_b;

  V_Axpby_Functor (const XV& x, const YV& y,
                   const typename XV::non_const_value_type& a,
                   const typename YV::non_const_value_type& b) :
    m_x (x), m_y (y), m_a (a), m_b (b)
  {
#ifdef KOKKOS_HAVE_CXX11
    static_assert (Kokkos::Impl::is_view<XV>::value, "KokkosBlas::Impl::"
                   "V_Axpby_Functor: X is not a Kokkos::View.");
    static_assert (Kokkos::Impl::is_view<YV>::value, "KokkosBlas::Impl::"
                   "V_Axpby_Functor: Y is not a Kokkos::View.");
    static_assert (Kokkos::Impl::is_same<typename YV::value_type,
                   typename YV::non_const_value_type>::value,
                   "KokkosBlas::Impl::V_Axpby_Functor: R is const.  "
                   "It must be nonconst, because it is an output argument "
                   "(we have to be able to write to its entries).");
    static_assert ((int) YV::rank == (int) XV::rank, "KokkosBlas::Impl::V_Axpby_Functor: "
                   "X and Y must have the same rank.");
    static_assert (YV::rank == 1, "KokkosBlas::Impl::V_Axpby_Functor: "
                   "XV and YV must have rank 1.");
#endif // KOKKOS_HAVE_CXX11
  }

  KOKKOS_INLINE_FUNCTION
  void operator() (const size_type& i) const
  {
    // scalar_x and scalar_y are compile-time constants (since they
    // are template parameters), so the compiler should evaluate these
    // branches at compile time.
    if (scalar_x == 0 && scalar_y == 0) {
      m_y(i) = ATS::zero ();
    }
    if (scalar_x == 0 && scalar_y == -1) {
      m_y(i) = -m_y(i);
    }
    if (scalar_x == 0 && scalar_y == 1) {
      return; // m_y(i) = m_y(i);
    }
    if (scalar_x == 0 && scalar_y == 2) {
      m_y(i) = m_b*m_y(i);
    }
    if (scalar_x == -1 && scalar_y == 0) {
      m_y(i) = -m_x(i);
    }
    if (scalar_x == -1 && scalar_y == -1) {
      m_y(i) = -m_x(i) - m_y(i);
    }
    if (scalar_x == -1 && scalar_y == 1) {
      m_y(i) = -m_x(i) + m_y(i);
    }
    if (scalar_x == -1 && scalar_y == 2) {
      m_y(i) = -m_x(i) + m_b*m_y(i);
    }
    if (scalar_x == 1 && scalar_y == 0) {
      m_y(i) = m_x(i);
    }
    if (scalar_x == 1 && scalar_y == -1) {
      m_y(i) = m_x(i) - m_y(i);
    }
    if (scalar_x == 1 && scalar_y == 1) {
      m_y(i) = m_x(i) + m_y(i);
    }
    if (scalar_x == 1 && scalar_y == 2) {
      m_y(i) = m_x(i) + m_b*m_y(i);
    }
    if (scalar_x == 2 && scalar_y == 0) {
      m_y(i) = m_a*m_x(i);
    }
    if (scalar_x == 2 && scalar_y == -1) {
      m_y(i) = m_a*m_x(i) - m_y(i);
    }
    if (scalar_x == 2 && scalar_y == 1) {
      m_y(i) = m_a*m_x(i) + m_y(i);
    }
    if (scalar_x == 2 && scalar_y == 2) {
      m_y(i) = m_a*m_x(i) + m_b*m_y(i);
    }
  }
};

// Invoke the unrolled multivector functor that computes any of the
// following:
//
// 1. Y(i,j) = a*X(i,j) + b*Y(i,j) for a,b in -1,0,1
// 2. Y(i,j) = av(j)*X(i,j) + b*Y(i,j) for b in -1,0,1
// 3. Y(i,j) = a*X(i,j) + b*Y(i,j) for a in -1,0,1
// 4. Y(i,j) = av(j)*X(i,j) + bv(j)*Y(i,j)
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
template<class AV, class XMV, class BV, class YMV,
         int UNROLL, class SizeType>
void
MV_Axpby_Unrolled (const AV& av, const XMV& x,
                   const BV& bv, const YMV& y,
                   int a = 2, int b = 2)
{
#ifdef KOKKOS_HAVE_CXX11
  static_assert (Kokkos::Impl::is_view<XMV>::value, "KokkosBlas::Impl::"
                 "MV_Axpby_Unrolled: X is not a Kokkos::View.");
  static_assert (Kokkos::Impl::is_view<YMV>::value, "KokkosBlas::Impl::"
                 "MV_Axpby_Unrolled: Y is not a Kokkos::View.");
  static_assert (Kokkos::Impl::is_same<typename YMV::value_type,
                   typename YMV::non_const_value_type>::value,
                 "KokkosBlas::Impl::MV_Axpby_Unrolled: Y is const.  "
                 "It must be nonconst, because it is an output argument "
                 "(we have to be able to write to its entries).");
  static_assert ((int) YMV::rank == (int) XMV::rank, "KokkosBlas::Impl::MV_Axpby_Unrolled: "
                 "X and Y must have the same rank.");
  static_assert (YMV::rank == 2, "KokkosBlas::Impl::MV_Axpby_Unrolled: "
                 "XMV and YMV must have rank 2.");
#endif // KOKKOS_HAVE_CXX11
  typedef typename YMV::execution_space execution_space;
  const SizeType numRows = x.dimension_0 ();
  Kokkos::RangePolicy<execution_space, SizeType> policy (0, numRows);

  if (a == 0 && b == 0) {
    MV_Axpby_Unroll_Functor<AV, XMV, BV, YMV, 0, 0, UNROLL, SizeType> op (x, y, av, bv);
    Kokkos::parallel_for (policy, op);
    return;
  }
  if (a == 0 && b == -1) {
    MV_Axpby_Unroll_Functor<AV, XMV, BV, YMV, 0, -1, UNROLL, SizeType> op (x, y, av, bv);
    Kokkos::parallel_for (policy, op);
    return;
  }
  if (a == 0 && b == 1) {
    MV_Axpby_Unroll_Functor<AV, XMV, BV, YMV, 0, 1, UNROLL, SizeType> op (x, y, av, bv);
    Kokkos::parallel_for (policy, op);
    return;
  }
  if (a == 0 && b == 2) {
    MV_Axpby_Unroll_Functor<AV, XMV, BV, YMV, 0, 2, UNROLL, SizeType> op (x, y, av, bv);
    Kokkos::parallel_for (policy, op);
    return;
  }
  // a == -1
  if (a == -1 && b == 0) {
    MV_Axpby_Unroll_Functor<AV, XMV, BV, YMV, -1, 0, UNROLL, SizeType> op (x, y, av, bv);
    Kokkos::parallel_for (policy, op);
    return;
  }
  if (a == -1 && b == -1) {
    MV_Axpby_Unroll_Functor<AV, XMV, BV, YMV, -1, -1, UNROLL, SizeType> op (x, y, av, bv);
    Kokkos::parallel_for (policy, op);
    return;
  }
  if (a == -1 && b == 1) {
    MV_Axpby_Unroll_Functor<AV, XMV, BV, YMV, -1, 1, UNROLL, SizeType> op (x, y, av, bv);
    Kokkos::parallel_for (policy, op);
    return;
  }
  if (a == -1 && b == 2) {
    MV_Axpby_Unroll_Functor<AV, XMV, BV, YMV, -1, 2, UNROLL, SizeType> op (x, y, av, bv);
    Kokkos::parallel_for (policy, op);
    return;
  }
  // a == 1
  if (a == 1 && b == 0) {
    MV_Axpby_Unroll_Functor<AV, XMV, BV, YMV, 1, 0, UNROLL, SizeType> op (x, y, av, bv);
    Kokkos::parallel_for (policy, op);
    return;
  }
  if (a == 1 && b == -1) {
    MV_Axpby_Unroll_Functor<AV, XMV, BV, YMV, 1, -1, UNROLL, SizeType> op (x, y, av, bv);
    Kokkos::parallel_for (policy, op);
    return;
  }
  if (a == 1 && b == 1) {
    MV_Axpby_Unroll_Functor<AV, XMV, BV, YMV, 1, 1, UNROLL, SizeType> op (x, y, av, bv);
    Kokkos::parallel_for (policy, op);
    return;
  }
  if (a == 1 && b == 2) {
    MV_Axpby_Unroll_Functor<AV, XMV, BV, YMV, 1, 2, UNROLL, SizeType> op (x, y, av, bv);
    Kokkos::parallel_for (policy, op);
    return;
  }
  // a == 2
  if (a == 2 && b == 0) {
    MV_Axpby_Unroll_Functor<AV, XMV, BV, YMV, 2, 0, UNROLL, SizeType> op (x, y, av, bv);
    Kokkos::parallel_for (policy, op);
    return;
  }
  if (a == 2 && b == -1) {
    MV_Axpby_Unroll_Functor<AV, XMV, BV, YMV, 2, -1, UNROLL, SizeType> op (x, y, av, bv);
    Kokkos::parallel_for (policy, op);
    return;
  }
  if (a == 2 && b == 1) {
    MV_Axpby_Unroll_Functor<AV, XMV, BV, YMV, 2, 1, UNROLL, SizeType> op (x, y, av, bv);
    Kokkos::parallel_for (policy, op);
    return;
  }

  // a and b arbitrary (not -1, 0, or 1)
  MV_Axpby_Unroll_Functor<AV, XMV, BV, YMV, 2, 2, UNROLL, SizeType> op (x, y, av, bv);
  Kokkos::parallel_for (policy, op);
}

// Invoke the "generic" (not unrolled) multivector functor that
// computes any of the following:
//
// 1. Y(i,j) = a*X(i,j) + b*Y(i,j) for a,b in -1,0,1
// 2. Y(i,j) = av(j)*X(i,j) + b*Y(i,j) for b in -1,0,1
// 3. Y(i,j) = a*X(i,j) + b*Y(i,j) for a in -1,0,1
// 4. Y(i,j) = av(j)*X(i,j) + bv(j)*Y(i,j)
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
template<class AV, class XMV, class BV, class YMV, class SizeType>
void
MV_Axpby_Generic (const AV& av, const XMV& x,
                  const BV& bv, const YMV& y,
                  int a = 2, int b = 2)
{
#ifdef KOKKOS_HAVE_CXX11
  static_assert (Kokkos::Impl::is_view<XMV>::value, "KokkosBlas::Impl::"
                 "MV_Axpby_Generic: X is not a Kokkos::View.");
  static_assert (Kokkos::Impl::is_view<YMV>::value, "KokkosBlas::Impl::"
                 "MV_Axpby_Generic: Y is not a Kokkos::View.");
  static_assert (Kokkos::Impl::is_same<typename YMV::value_type,
                   typename YMV::non_const_value_type>::value,
                 "KokkosBlas::Impl::MV_Axpby_Generic: Y is const.  "
                 "It must be nonconst, because it is an output argument "
                 "(we have to be able to write to its entries).");
  static_assert ((int) YMV::rank == (int) XMV::rank, "KokkosBlas::Impl::MV_Axpby_Generic: "
                 "X and Y must have the same rank.");
  static_assert (YMV::rank == 2, "KokkosBlas::Impl::MV_Axpby_Generic: "
                 "XMV and YMV must have rank 2.");
#endif // KOKKOS_HAVE_CXX11
  typedef typename YMV::execution_space execution_space;
  const SizeType numRows = x.dimension_0 ();
  Kokkos::RangePolicy<execution_space, SizeType> policy (0, numRows);

  if (a == 0 && b == 0) {
    MV_Axpby_Functor<AV, XMV, BV, YMV, 0, 0, SizeType> op (x, y, av, bv);
    Kokkos::parallel_for (policy, op);
    return;
  }
  if (a == 0 && b == -1) {
    MV_Axpby_Functor<AV, XMV, BV, YMV, 0, -1, SizeType> op (x, y, av, bv);
    Kokkos::parallel_for (policy, op);
    return;
  }
  if (a == 0 && b == 1) {
    MV_Axpby_Functor<AV, XMV, BV, YMV, 0, 1, SizeType> op (x, y, av, bv);
    Kokkos::parallel_for (policy, op);
    return;
  }
  if (a == 0 && b == 2) {
    MV_Axpby_Functor<AV, XMV, BV, YMV, 0, 2, SizeType> op (x, y, av, bv);
    Kokkos::parallel_for (policy, op);
    return;
  }
  // a == -1
  if (a == -1 && b == 0) {
    MV_Axpby_Functor<AV, XMV, BV, YMV, -1, 0, SizeType> op (x, y, av, bv);
    Kokkos::parallel_for (policy, op);
    return;
  }
  if (a == -1 && b == -1) {
    MV_Axpby_Functor<AV, XMV, BV, YMV, -1, -1, SizeType> op (x, y, av, bv);
    Kokkos::parallel_for (policy, op);
    return;
  }
  if (a == -1 && b == 1) {
    MV_Axpby_Functor<AV, XMV, BV, YMV, -1, 1, SizeType> op (x, y, av, bv);
    Kokkos::parallel_for (policy, op);
    return;
  }
  if (a == -1 && b == 2) {
    MV_Axpby_Functor<AV, XMV, BV, YMV, -1, 2, SizeType> op (x, y, av, bv);
    Kokkos::parallel_for (policy, op);
    return;
  }
  // a == 1
  if (a == 1 && b == 0) {
    MV_Axpby_Functor<AV, XMV, BV, YMV, 1, 0, SizeType> op (x, y, av, bv);
    Kokkos::parallel_for (policy, op);
    return;
  }
  if (a == 1 && b == -1) {
    MV_Axpby_Functor<AV, XMV, BV, YMV, 1, -1, SizeType> op (x, y, av, bv);
    Kokkos::parallel_for (policy, op);
    return;
  }
  if (a == 1 && b == 1) {
    MV_Axpby_Functor<AV, XMV, BV, YMV, 1, 1, SizeType> op (x, y, av, bv);
    Kokkos::parallel_for (policy, op);
    return;
  }
  if (a == 1 && b == 2) {
    MV_Axpby_Functor<AV, XMV, BV, YMV, 1, 2, SizeType> op (x, y, av, bv);
    Kokkos::parallel_for (policy, op);
    return;
  }
  // a == 2
  if (a == 2 && b == 0) {
    MV_Axpby_Functor<AV, XMV, BV, YMV, 2, 0, SizeType> op (x, y, av, bv);
    Kokkos::parallel_for (policy, op);
    return;
  }
  if (a == 2 && b == -1) {
    MV_Axpby_Functor<AV, XMV, BV, YMV, 2, -1, SizeType> op (x, y, av, bv);
    Kokkos::parallel_for (policy, op);
    return;
  }
  if (a == 2 && b == 1) {
    MV_Axpby_Functor<AV, XMV, BV, YMV, 2, 1, SizeType> op (x, y, av, bv);
    Kokkos::parallel_for (policy, op);
    return;
  }

  // a and b arbitrary (not -1, 0, or 1)
  MV_Axpby_Functor<AV, XMV, BV, YMV, 2, 2, SizeType> op (x, y, av, bv);
  Kokkos::parallel_for (policy, op);
}

// Variant of MV_Axpby_Generic for single vectors (1-D Views) x and y.
// As above, either av and bv are both 1-D Views (and only the first
// entry of each will be read), or both av and bv are scalars.
template<class AV, class XV, class BV, class YV, class SizeType>
void
V_Axpby_Generic (const AV& av, const XV& x,
                 const BV& bv, const YV& y,
                 int a = 2, int b = 2)
{
#ifdef KOKKOS_HAVE_CXX11
  static_assert (Kokkos::Impl::is_view<XV>::value, "KokkosBlas::Impl::"
                 "V_Axpby_Generic: X is not a Kokkos::View.");
  static_assert (Kokkos::Impl::is_view<YV>::value, "KokkosBlas::Impl::"
                 "V_Axpby_Generic: Y is not a Kokkos::View.");
  static_assert (Kokkos::Impl::is_same<typename YV::value_type,
                   typename YV::non_const_value_type>::value,
                 "KokkosBlas::Impl::V_Axpby_Generic: Y is const.  "
                 "It must be nonconst, because it is an output argument "
                 "(we have to be able to write to its entries).");
  static_assert ((int) YV::rank == (int) XV::rank, "KokkosBlas::Impl::V_Axpby_Generic: "
                 "X and Y must have the same rank.");
  static_assert (YV::rank == 1, "KokkosBlas::Impl::V_Axpby_Generic: "
                 "XV and YV must have rank 1.");
#endif // KOKKOS_HAVE_CXX11

  typedef typename YV::execution_space execution_space;
  const SizeType numRows = x.dimension_0 ();
  Kokkos::RangePolicy<execution_space, SizeType> policy (0, numRows);

  if (a == 0 && b == 0) {
    V_Axpby_Functor<AV, XV, BV, YV, 0, 0, SizeType> op (x, y, av, bv);
    Kokkos::parallel_for (policy, op);
    return;
  }
  if (a == 0 && b == -1) {
    V_Axpby_Functor<AV, XV, BV, YV, 0, -1, SizeType> op (x, y, av, bv);
    Kokkos::parallel_for (policy, op);
    return;
  }
  if (a == 0 && b == 1) {
    V_Axpby_Functor<AV, XV, BV, YV, 0, 1, SizeType> op (x, y, av, bv);
    Kokkos::parallel_for (policy, op);
    return;
  }
  if (a == 0 && b == 2) {
    V_Axpby_Functor<AV, XV, BV, YV, 0, 2, SizeType> op (x, y, av, bv);
    Kokkos::parallel_for (policy, op);
    return;
  }
  // a == -1
  if (a == -1 && b == 0) {
    V_Axpby_Functor<AV, XV, BV, YV, -1, 0, SizeType> op (x, y, av, bv);
    Kokkos::parallel_for (policy, op);
    return;
  }
  if (a == -1 && b == -1) {
    V_Axpby_Functor<AV, XV, BV, YV, -1, -1, SizeType> op (x, y, av, bv);
    Kokkos::parallel_for (policy, op);
    return;
  }
  if (a == -1 && b == 1) {
    V_Axpby_Functor<AV, XV, BV, YV, -1, 1, SizeType> op (x, y, av, bv);
    Kokkos::parallel_for (policy, op);
    return;
  }
  if (a == -1 && b == 2) {
    V_Axpby_Functor<AV, XV, BV, YV, -1, 2, SizeType> op (x, y, av, bv);
    Kokkos::parallel_for (policy, op);
    return;
  }
  // a == 1
  if (a == 1 && b == 0) {
    V_Axpby_Functor<AV, XV, BV, YV, 1, 0, SizeType> op (x, y, av, bv);
    Kokkos::parallel_for (policy, op);
    return;
  }
  if (a == 1 && b == -1) {
    V_Axpby_Functor<AV, XV, BV, YV, 1, -1, SizeType> op (x, y, av, bv);
    Kokkos::parallel_for (policy, op);
    return;
  }
  if (a == 1 && b == 1) {
    V_Axpby_Functor<AV, XV, BV, YV, 1, 1, SizeType> op (x, y, av, bv);
    Kokkos::parallel_for (policy, op);
    return;
  }
  if (a == 1 && b == 2) {
    V_Axpby_Functor<AV, XV, BV, YV, 1, 2, SizeType> op (x, y, av, bv);
    Kokkos::parallel_for (policy, op);
    return;
  }
  // a == 2
  if (a == 2 && b == 0) {
    V_Axpby_Functor<AV, XV, BV, YV, 2, 0, SizeType> op (x, y, av, bv);
    Kokkos::parallel_for (policy, op);
    return;
  }
  if (a == 2 && b == -1) {
    V_Axpby_Functor<AV, XV, BV, YV, 2, -1, SizeType> op (x, y, av, bv);
    Kokkos::parallel_for (policy, op);
    return;
  }
  if (a == 2 && b == 1) {
    V_Axpby_Functor<AV, XV, BV, YV, 2, 1, SizeType> op (x, y, av, bv);
    Kokkos::parallel_for (policy, op);
    return;
  }

  // a and b arbitrary (not -1, 0, or 1)
  V_Axpby_Functor<AV, XV, BV, YV, 2, 2, SizeType> op (x, y, av, bv);
  Kokkos::parallel_for (policy, op);
}

// Compute any of the following, in a way optimized for X and Y
// being LayoutLeft:
//
// 1. Y(i,j) = a*X(i,j) + b*Y(i,j) for a,b in -1,0,1
// 2. Y(i,j) = av(j)*X(i,j) + b*Y(i,j) for b in -1,0,1
// 3. Y(i,j) = a*X(i,j) + b*Y(i,j) for a in -1,0,1
// 4. Y(i,j) = av(j)*X(i,j) + bv(j)*Y(i,j)
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
template<class AV, class XMV, class BV, class YMV, class SizeType>
void
MV_Axpby_Invoke_Left (const AV& av, const XMV& x,
                      const BV& bv, const YMV& y,
                      int a = 2, int b = 2)
{
#ifdef KOKKOS_HAVE_CXX11
    static_assert (Kokkos::Impl::is_view<XMV>::value, "KokkosBlas::Impl::MV_Axpby_Invoke_Left (MV): "
                   "X is not a Kokkos::View.");
    static_assert (Kokkos::Impl::is_view<YMV>::value, "KokkosBlas::Impl::MV_Axpby_Invoke_Left (MV): "
                   "Y is not a Kokkos::View.");
    static_assert (Kokkos::Impl::is_same<typename YMV::value_type,
                     typename YMV::non_const_value_type>::value,
                   "KokkosBlas::Impl::MV_Axpby_Invoke_Left (MV): Y is const.  "
                   "It must be nonconst, because it is an output argument "
                   "(we have to be able to write to its entries).");
    static_assert ((int) YMV::rank == (int) XMV::rank, "KokkosBlas::Impl::MV_Axpby_Invoke_Left (MV): "
                   "X and Y must have the same rank.");
    static_assert (YMV::rank == 2, "KokkosBlas::Impl::MV_Axpby_Invoke_Left (MV): "
                   "X and Y must have rank 2.");
#endif // KOKKOS_HAVE_CXX11

  const SizeType numCols = x.dimension_1 ();
  switch (numCols) {
  case 1: {
    typedef Kokkos::View<typename XMV::value_type*, typename XMV::array_layout,
      typename XMV::device_type, typename XMV::memory_traits,
      typename XMV::specialize> XV;
    typedef Kokkos::View<typename YMV::value_type*, typename YMV::array_layout,
      typename YMV::device_type, typename YMV::memory_traits,
      typename YMV::specialize> YV;

    XV x_0 = Kokkos::subview (x, Kokkos::ALL (), 0);
    YV y_0 = Kokkos::subview (y, Kokkos::ALL (), 0);
    V_Axpby_Generic<AV, XV, BV, YV, SizeType> (av, x_0, bv, y_0, a, b);
    break;
  }
  case 2:
    MV_Axpby_Unrolled<AV, XMV, BV, YMV, 2, SizeType> (av, x, bv, y, a, b);
    break;
  case 3:
    MV_Axpby_Unrolled<AV, XMV, BV, YMV, 3, SizeType> (av, x, bv, y, a, b);
    break;
  case 4:
    MV_Axpby_Unrolled<AV, XMV, BV, YMV, 4, SizeType> (av, x, bv, y, a, b);
    break;
  case 5:
    MV_Axpby_Unrolled<AV, XMV, BV, YMV, 5, SizeType> (av, x, bv, y, a, b);
    break;
  case 6:
    MV_Axpby_Unrolled<AV, XMV, BV, YMV, 6, SizeType> (av, x, bv, y, a, b);
    break;
  case 7:
    MV_Axpby_Unrolled<AV, XMV, BV, YMV, 7, SizeType> (av, x, bv, y, a, b);
    break;
  case 8:
    MV_Axpby_Unrolled<AV, XMV, BV, YMV, 8, SizeType> (av, x, bv, y, a, b);
    break;
  case 9:
    MV_Axpby_Unrolled<AV, XMV, BV, YMV, 9, SizeType> (av, x, bv, y, a, b);
    break;
  case 10:
    MV_Axpby_Unrolled<AV, XMV, BV, YMV, 10, SizeType> (av, x, bv, y, a, b);
    break;
  case 11:
    MV_Axpby_Unrolled<AV, XMV, BV, YMV, 11, SizeType> (av, x, bv, y, a, b);
    break;
  case 12:
    MV_Axpby_Unrolled<AV, XMV, BV, YMV, 12, SizeType> (av, x, bv, y, a, b);
    break;
  case 13:
    MV_Axpby_Unrolled<AV, XMV, BV, YMV, 13, SizeType> (av, x, bv, y, a, b);
    break;
  case 14:
    MV_Axpby_Unrolled<AV, XMV, BV, YMV, 14, SizeType> (av, x, bv, y, a, b);
    break;
  case 15:
    MV_Axpby_Unrolled<AV, XMV, BV, YMV, 15, SizeType> (av, x, bv, y, a, b);
    break;
  case 16:
    MV_Axpby_Unrolled<AV, XMV, BV, YMV, 16, SizeType> (av, x, bv, y, a, b);
    break;
  default:
    MV_Axpby_Generic<AV, XMV, BV, YMV, SizeType> (av, x, bv, y, a, b);
  }
}

// Compute any of the following, in a way optimized for X, Y, and R
// being LayoutRight:
//
// 1. Y(i,j) = a*X(i,j) + b*Y(i,j) for a,b in -1,0,1
// 2. Y(i,j) = av(j)*X(i,j) + b*Y(i,j) for b in -1,0,1
// 3. Y(i,j) = a*X(i,j) + b*Y(i,j) for a in -1,0,1
// 4. Y(i,j) = av(j)*X(i,j) + bv(j)*Y(i,j)
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
template<class AV, class XMV, class BV, class YMV, class SizeType>
void
MV_Axpby_Invoke_Right (const AV& av, const XMV& x,
                       const BV& bv, const YMV& y,
                       int a = 2, int b = 2)
{
#ifdef KOKKOS_HAVE_CXX11
    static_assert (Kokkos::Impl::is_view<XMV>::value,
                   "KokkosBlas::Impl::MV_Axpby_Invoke_Right (MV): "
                   "X is not a Kokkos::View.");
    static_assert (Kokkos::Impl::is_view<YMV>::value,
                   "KokkosBlas::Impl::MV_Axpby_Invoke_Right (MV): "
                   "Y is not a Kokkos::View.");
    static_assert (Kokkos::Impl::is_same<typename YMV::value_type,
                     typename YMV::non_const_value_type>::value,
                   "KokkosBlas::Impl::MV_Axpby_Invoke_Right (MV): Y is const.  "
                   "It must be nonconst, because it is an output argument "
                   "(we have to be able to write to its entries).");
    static_assert ((int) YMV::rank == (int) XMV::rank,
                   "KokkosBlas::Impl::MV_Axpby_Invoke_Right (MV): "
                   "X and Y must have the same rank.");
    static_assert (YMV::rank == 2, "KokkosBlas::Impl::MV_Axpby_Invoke_Right (MV): "
                   "X and Y must have rank 2.");
#endif // KOKKOS_HAVE_CXX11

  const SizeType numCols = x.dimension_1 ();
  if (numCols == 1) {
    typedef Kokkos::View<typename XMV::value_type*, typename XMV::array_layout,
      typename XMV::device_type, typename XMV::memory_traits,
      typename XMV::specialize> XV;
    typedef Kokkos::View<typename YMV::value_type*, typename YMV::array_layout,
      typename YMV::device_type, typename YMV::memory_traits,
      typename YMV::specialize> YV;

    XV x_0 = Kokkos::subview (x, Kokkos::ALL (), 0);
    YV y_0 = Kokkos::subview (y, Kokkos::ALL (), 0);
    V_Axpby_Generic<AV, XMV, BV, YMV, 1, SizeType> (av, x_0, bv, y_0, a, b);
  }
  else {
    MV_Axpby_Generic<AV, XMV, BV, YMV, SizeType> (av, x, bv, y, a, b);
  }
}

/// \brief Implementation of KokkosBlas::axpby for (multi)vectors.
///
/// Compute any of the following, depending on the types of the input
/// arguments of axpxy():
///
/// 1. Y(i,j) = av(j)*X(i,j) + bv(j)*Y(i,j) (if R, X, and Y are 2-D,
///    and av and bv are 1-D)
///
/// 2. Y(i,j) = av*X(i,j) + bv*Y(i,j) (if R, X, and Y are 2-D,
///    and av and bv are scalars)
///
/// 3. Y(i) = av()*X(i) + bv()*Y(i) (if R, X, and Y are 1-D, and av
///    and bv are 0-D Views (not scalars))
///
/// 4. Y(i) = av*X(i) + bv*Y(i) (if R, X, and Y are 1-D, and av and bv
///    are scalars)
///
// Any <i>scalar</i> coefficient of zero has BLAS semantics of
/// ignoring the corresponding (multi)vector entry.  This does NOT
/// apply to coefficients in av and bv vectors, if they are used.
template<class AV, class XMV, class BV, class YMV, int rank = YMV::rank>
struct Axpby {};

// Partial specialization for XMV and YMV rank-2 Views.
template<class AV, class XMV, class BV, class YMV>
struct Axpby<AV, XMV, BV, YMV, 2>
{
  typedef typename YMV::size_type size_type;

  static void
  axpby (const AV& av, const XMV& X, const BV& bv, const YMV& Y)
  {
#ifdef KOKKOS_HAVE_CXX11
    static_assert (Kokkos::Impl::is_view<XMV>::value,
                   "KokkosBlas::Impl::Axpby::axpby (MV): "
                   "X is not a Kokkos::View.");
    static_assert (Kokkos::Impl::is_view<YMV>::value,
                   "KokkosBlas::Impl::Axpby::axpby (MV): "
                   "Y is not a Kokkos::View.");
    static_assert (Kokkos::Impl::is_same<typename YMV::value_type,
                     typename YMV::non_const_value_type>::value,
                   "KokkosBlas::Impl::Axpby::axpby (MV): Y is const.  "
                   "It must be nonconst, because it is an output argument "
                   "(we have to be able to write to its entries).");
    static_assert ((int) YMV::rank == (int) XMV::rank,
                   "KokkosBlas::Impl::Axpby::axpby (MV): "
                   "X and Y must have the same rank.");
    static_assert (YMV::rank == 2, "KokkosBlas::Impl::Axpby::axpby (MV): "
                   "X and Y must have rank 2.");
#endif // KOKKOS_HAVE_CXX11

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
      MV_Axpby_Invoke_Left<AV, XMV, BV, YMV, index_type> (av, X, bv, Y, a, b);
    }
    else {
      typedef typename XMV::size_type index_type;
      MV_Axpby_Invoke_Left<AV, XMV, BV, YMV, index_type> (av, X, bv, Y, a, b);
    }
  }
};

// Partial specialization for XMV, and YMV rank-2 Views,
// and AV and BV scalars.
template<class XMV, class YMV>
struct Axpby<typename XMV::non_const_value_type, XMV,
             typename YMV::non_const_value_type, YMV, 2>
{
  typedef typename XMV::non_const_value_type AV;
  typedef typename YMV::non_const_value_type BV;
  typedef typename YMV::size_type size_type;
  typedef Kokkos::Details::ArithTraits<typename XMV::non_const_value_type> ATA;
  typedef Kokkos::Details::ArithTraits<typename YMV::non_const_value_type> ATB;

  static void
  axpby (const AV& alpha, const XMV& X, const BV& beta, const YMV& Y)
  {
#ifdef KOKKOS_HAVE_CXX11
    static_assert (Kokkos::Impl::is_view<XMV>::value,
                   "KokkosBlas::Impl::Axpby::axpby (MV): "
                   "X is not a Kokkos::View.");
    static_assert (Kokkos::Impl::is_view<YMV>::value,
                   "KokkosBlas::Impl::Axpby::axpby (MV): "
                   "Y is not a Kokkos::View.");
    static_assert (Kokkos::Impl::is_same<typename YMV::value_type,
                     typename YMV::non_const_value_type>::value,
                   "KokkosBlas::Impl::Axpby::axpby (MV): Y is const.  "
                   "It must be nonconst, because it is an output argument "
                   "(we have to be able to write to its entries).");
    static_assert ((int) YMV::rank == (int) XMV::rank,
                   "KokkosBlas::Impl::Axpby::axpby (MV): "
                   "X and Y must have the same rank.");
    static_assert (YMV::rank == 2, "KokkosBlas::Impl::Axpby::axpby (MV): "
                   "X and Y must have rank 2.");
#endif // KOKKOS_HAVE_CXX11
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
      MV_Axpby_Invoke_Left<AV, XMV, BV, YMV, index_type> (alpha, X,
                                                          beta, Y, a, b);
    }
    else {
      typedef typename XMV::size_type index_type;
      MV_Axpby_Invoke_Left<AV, XMV, BV, YMV, index_type> (alpha, X,
                                                          beta, Y, a, b);
    }
  }
};

// Partial specialization for XV and YV rank-1 Views,
// and AV and BV scalars.
template<class XV, class YV>
struct Axpby<typename XV::non_const_value_type, XV,
             typename YV::non_const_value_type, YV, 1>
{
  typedef typename XV::non_const_value_type AV;
  typedef typename YV::non_const_value_type BV;
  typedef typename YV::size_type size_type;
  typedef Kokkos::Details::ArithTraits<typename XV::non_const_value_type> ATA;
  typedef Kokkos::Details::ArithTraits<typename YV::non_const_value_type> ATB;

  static void
  axpby (const AV& alpha, const XV& X, const BV& beta, const YV& Y)
  {
#ifdef KOKKOS_HAVE_CXX11
    static_assert (Kokkos::Impl::is_view<XV>::value,
                   "KokkosBlas::Impl::Axpby::axpby (V): "
                   "X is not a Kokkos::View.");
    static_assert (Kokkos::Impl::is_view<YV>::value,
                   "KokkosBlas::Impl::Axpby::axpby (V): "
                   "Y is not a Kokkos::View.");
    static_assert (Kokkos::Impl::is_same<typename YV::value_type,
                     typename YV::non_const_value_type>::value,
                   "KokkosBlas::Impl::Axpby::axpby (V): Y is const.  "
                   "It must be nonconst, because it is an output argument "
                   "(we have to be able to write to its entries).");
    static_assert ((int) YV::rank == (int) XV::rank,
                   "KokkosBlas::Impl::Axpby::axpby (V): "
                   "X and Y must have the same rank.");
    static_assert (YV::rank == 1, "KokkosBlas::Impl::Axpby::axpby (V): "
                   "X and Y must have rank 1.");
#endif // KOKKOS_HAVE_CXX11

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
    int b = 2;
    if (beta == ATB::zero ()) {
      b = 0;
    }
    else if (beta == -ATB::one ()) {
      b = -1;
    }
    else if (beta == ATB::one ()) {
      b = 1;
    }

    if (numRows < static_cast<size_type> (INT_MAX)) {
      typedef int index_type;
      V_Axpby_Generic<typename XV::non_const_value_type, XV,
        typename YV::non_const_value_type, YV,
        index_type> (alpha, X, beta, Y, a, b);
    }
    else {
      typedef typename XV::size_type index_type;
      V_Axpby_Generic<typename XV::non_const_value_type, XV,
        typename YV::non_const_value_type, YV,
        index_type> (alpha, X, beta, Y, a, b);
    }
  }
};

//
// Declarations of full specializations of Impl::Axpby for rank == 2.
// Their definitions go in .cpp file(s) in this source directory.
//

#ifdef KOKKOS_HAVE_SERIAL
#define KOKKOSBLAS_IMPL_MV_EXEC_SPACE Kokkos::Serial
#define KOKKOSBLAS_IMPL_MV_MEM_SPACE Kokkos::HostSpace
#define KOKKOSBLAS_IMPL_MV_SCALAR double

template<>
struct Axpby<KOKKOSBLAS_IMPL_MV_SCALAR,
             Kokkos::View<const KOKKOSBLAS_IMPL_MV_SCALAR**,
                          Kokkos::LayoutLeft,
                          Kokkos::Device<KOKKOSBLAS_IMPL_MV_EXEC_SPACE, KOKKOSBLAS_IMPL_MV_MEM_SPACE>,
                          Kokkos::MemoryTraits<Kokkos::Unmanaged>,
                          Kokkos::Impl::ViewDefault>,
             KOKKOSBLAS_IMPL_MV_SCALAR,
             Kokkos::View<KOKKOSBLAS_IMPL_MV_SCALAR**,
                          Kokkos::LayoutLeft,
                          Kokkos::Device<KOKKOSBLAS_IMPL_MV_EXEC_SPACE, KOKKOSBLAS_IMPL_MV_MEM_SPACE>,
                          Kokkos::MemoryTraits<Kokkos::Unmanaged>,
                          Kokkos::Impl::ViewDefault>, 2>
{
  typedef KOKKOSBLAS_IMPL_MV_SCALAR AV;
  typedef Kokkos::View<const KOKKOSBLAS_IMPL_MV_SCALAR**,
                       Kokkos::LayoutLeft,
                       Kokkos::Device<KOKKOSBLAS_IMPL_MV_EXEC_SPACE, KOKKOSBLAS_IMPL_MV_MEM_SPACE>,
                       Kokkos::MemoryTraits<Kokkos::Unmanaged>,
                       Kokkos::Impl::ViewDefault> XMV;
  typedef KOKKOSBLAS_IMPL_MV_SCALAR BV;
  typedef Kokkos::View<KOKKOSBLAS_IMPL_MV_SCALAR**,
                       Kokkos::LayoutLeft,
                       Kokkos::Device<KOKKOSBLAS_IMPL_MV_EXEC_SPACE, KOKKOSBLAS_IMPL_MV_MEM_SPACE>,
                       Kokkos::MemoryTraits<Kokkos::Unmanaged>,
                       Kokkos::Impl::ViewDefault> YMV;
  typedef YMV::size_type size_type;
  typedef Kokkos::Details::ArithTraits<XMV::non_const_value_type> ATA;
  typedef Kokkos::Details::ArithTraits<YMV::non_const_value_type> ATB;

  static void
  axpby (const XMV::non_const_value_type& alpha, const XMV& X,
         const YMV::non_const_value_type& beta, const YMV& Y);
};

#undef KOKKOSBLAS_IMPL_MV_SCALAR
#undef KOKKOSBLAS_IMPL_MV_MEM_SPACE
#undef KOKKOSBLAS_IMPL_MV_EXEC_SPACE
#endif // KOKKOS_HAVE_SERIAL

#ifdef KOKKOS_HAVE_OPENMP
#define KOKKOSBLAS_IMPL_MV_EXEC_SPACE Kokkos::OpenMP
#define KOKKOSBLAS_IMPL_MV_MEM_SPACE Kokkos::HostSpace
#define KOKKOSBLAS_IMPL_MV_SCALAR double

template<>
struct Axpby<KOKKOSBLAS_IMPL_MV_SCALAR,
             Kokkos::View<const KOKKOSBLAS_IMPL_MV_SCALAR**,
                          Kokkos::LayoutLeft,
                          Kokkos::Device<KOKKOSBLAS_IMPL_MV_EXEC_SPACE, KOKKOSBLAS_IMPL_MV_MEM_SPACE>,
                          Kokkos::MemoryTraits<Kokkos::Unmanaged>,
                          Kokkos::Impl::ViewDefault>,
             KOKKOSBLAS_IMPL_MV_SCALAR,
             Kokkos::View<KOKKOSBLAS_IMPL_MV_SCALAR**,
                          Kokkos::LayoutLeft,
                          Kokkos::Device<KOKKOSBLAS_IMPL_MV_EXEC_SPACE, KOKKOSBLAS_IMPL_MV_MEM_SPACE>,
                          Kokkos::MemoryTraits<Kokkos::Unmanaged>,
                          Kokkos::Impl::ViewDefault>, 2>
{
  typedef KOKKOSBLAS_IMPL_MV_SCALAR AV;
  typedef Kokkos::View<const KOKKOSBLAS_IMPL_MV_SCALAR**,
                       Kokkos::LayoutLeft,
                       Kokkos::Device<KOKKOSBLAS_IMPL_MV_EXEC_SPACE, KOKKOSBLAS_IMPL_MV_MEM_SPACE>,
                       Kokkos::MemoryTraits<Kokkos::Unmanaged>,
                       Kokkos::Impl::ViewDefault> XMV;
  typedef KOKKOSBLAS_IMPL_MV_SCALAR BV;
  typedef Kokkos::View<KOKKOSBLAS_IMPL_MV_SCALAR**,
                       Kokkos::LayoutLeft,
                       Kokkos::Device<KOKKOSBLAS_IMPL_MV_EXEC_SPACE, KOKKOSBLAS_IMPL_MV_MEM_SPACE>,
                       Kokkos::MemoryTraits<Kokkos::Unmanaged>,
                       Kokkos::Impl::ViewDefault> YMV;
  typedef YMV::size_type size_type;
  typedef Kokkos::Details::ArithTraits<XMV::non_const_value_type> ATA;
  typedef Kokkos::Details::ArithTraits<YMV::non_const_value_type> ATB;

  static void
  axpby (const XMV::non_const_value_type& alpha, const XMV& X,
         const YMV::non_const_value_type& beta, const YMV& Y);
};

#undef KOKKOSBLAS_IMPL_MV_SCALAR
#undef KOKKOSBLAS_IMPL_MV_MEM_SPACE
#undef KOKKOSBLAS_IMPL_MV_EXEC_SPACE
#endif // KOKKOS_HAVE_OPENMP

#ifdef KOKKOS_HAVE_PTHREAD
#define KOKKOSBLAS_IMPL_MV_EXEC_SPACE Kokkos::Threads
#define KOKKOSBLAS_IMPL_MV_MEM_SPACE Kokkos::HostSpace
#define KOKKOSBLAS_IMPL_MV_SCALAR double

template<>
struct Axpby<KOKKOSBLAS_IMPL_MV_SCALAR,
             Kokkos::View<const KOKKOSBLAS_IMPL_MV_SCALAR**,
                          Kokkos::LayoutLeft,
                          Kokkos::Device<KOKKOSBLAS_IMPL_MV_EXEC_SPACE, KOKKOSBLAS_IMPL_MV_MEM_SPACE>,
                          Kokkos::MemoryTraits<Kokkos::Unmanaged>,
                          Kokkos::Impl::ViewDefault>,
             KOKKOSBLAS_IMPL_MV_SCALAR,
             Kokkos::View<KOKKOSBLAS_IMPL_MV_SCALAR**,
                          Kokkos::LayoutLeft,
                          Kokkos::Device<KOKKOSBLAS_IMPL_MV_EXEC_SPACE, KOKKOSBLAS_IMPL_MV_MEM_SPACE>,
                          Kokkos::MemoryTraits<Kokkos::Unmanaged>,
                          Kokkos::Impl::ViewDefault>, 2>
{
  typedef KOKKOSBLAS_IMPL_MV_SCALAR AV;
  typedef Kokkos::View<const KOKKOSBLAS_IMPL_MV_SCALAR**,
                       Kokkos::LayoutLeft,
                       Kokkos::Device<KOKKOSBLAS_IMPL_MV_EXEC_SPACE, KOKKOSBLAS_IMPL_MV_MEM_SPACE>,
                       Kokkos::MemoryTraits<Kokkos::Unmanaged>,
                       Kokkos::Impl::ViewDefault> XMV;
  typedef KOKKOSBLAS_IMPL_MV_SCALAR BV;
  typedef Kokkos::View<KOKKOSBLAS_IMPL_MV_SCALAR**,
                       Kokkos::LayoutLeft,
                       Kokkos::Device<KOKKOSBLAS_IMPL_MV_EXEC_SPACE, KOKKOSBLAS_IMPL_MV_MEM_SPACE>,
                       Kokkos::MemoryTraits<Kokkos::Unmanaged>,
                       Kokkos::Impl::ViewDefault> YMV;
  typedef YMV::size_type size_type;
  typedef Kokkos::Details::ArithTraits<XMV::non_const_value_type> ATA;
  typedef Kokkos::Details::ArithTraits<YMV::non_const_value_type> ATB;

  static void
  axpby (const XMV::non_const_value_type& alpha, const XMV& X,
         const YMV::non_const_value_type& beta, const YMV& Y);
};

#undef KOKKOSBLAS_IMPL_MV_SCALAR
#undef KOKKOSBLAS_IMPL_MV_MEM_SPACE
#undef KOKKOSBLAS_IMPL_MV_EXEC_SPACE
#endif // KOKKOS_HAVE_PTHREAD

#ifdef KOKKOS_HAVE_CUDA
#define KOKKOSBLAS_IMPL_MV_EXEC_SPACE Kokkos::Cuda
#define KOKKOSBLAS_IMPL_MV_MEM_SPACE Kokkos::CudaSpace
#define KOKKOSBLAS_IMPL_MV_SCALAR double

template<>
struct Axpby<KOKKOSBLAS_IMPL_MV_SCALAR,
             Kokkos::View<const KOKKOSBLAS_IMPL_MV_SCALAR**,
                          Kokkos::LayoutLeft,
                          Kokkos::Device<KOKKOSBLAS_IMPL_MV_EXEC_SPACE, KOKKOSBLAS_IMPL_MV_MEM_SPACE>,
                          Kokkos::MemoryTraits<Kokkos::Unmanaged>,
                          Kokkos::Impl::ViewDefault>,
             KOKKOSBLAS_IMPL_MV_SCALAR,
             Kokkos::View<KOKKOSBLAS_IMPL_MV_SCALAR**,
                          Kokkos::LayoutLeft,
                          Kokkos::Device<KOKKOSBLAS_IMPL_MV_EXEC_SPACE, KOKKOSBLAS_IMPL_MV_MEM_SPACE>,
                          Kokkos::MemoryTraits<Kokkos::Unmanaged>,
                          Kokkos::Impl::ViewDefault>, 2>
{
  typedef KOKKOSBLAS_IMPL_MV_SCALAR AV;
  typedef Kokkos::View<const KOKKOSBLAS_IMPL_MV_SCALAR**,
                       Kokkos::LayoutLeft,
                       Kokkos::Device<KOKKOSBLAS_IMPL_MV_EXEC_SPACE, KOKKOSBLAS_IMPL_MV_MEM_SPACE>,
                       Kokkos::MemoryTraits<Kokkos::Unmanaged>,
                       Kokkos::Impl::ViewDefault> XMV;
  typedef KOKKOSBLAS_IMPL_MV_SCALAR BV;
  typedef Kokkos::View<KOKKOSBLAS_IMPL_MV_SCALAR**,
                       Kokkos::LayoutLeft,
                       Kokkos::Device<KOKKOSBLAS_IMPL_MV_EXEC_SPACE, KOKKOSBLAS_IMPL_MV_MEM_SPACE>,
                       Kokkos::MemoryTraits<Kokkos::Unmanaged>,
                       Kokkos::Impl::ViewDefault> YMV;
  typedef YMV::size_type size_type;
  typedef Kokkos::Details::ArithTraits<XMV::non_const_value_type> ATA;
  typedef Kokkos::Details::ArithTraits<YMV::non_const_value_type> ATB;

  static void
  axpby (const XMV::non_const_value_type& alpha, const XMV& X,
         const YMV::non_const_value_type& beta, const YMV& Y);
};

#undef KOKKOSBLAS_IMPL_MV_SCALAR
#undef KOKKOSBLAS_IMPL_MV_MEM_SPACE
#undef KOKKOSBLAS_IMPL_MV_EXEC_SPACE
#endif // KOKKOS_HAVE_CUDA

#ifdef KOKKOS_HAVE_CUDA
#define KOKKOSBLAS_IMPL_MV_EXEC_SPACE Kokkos::Cuda
#define KOKKOSBLAS_IMPL_MV_MEM_SPACE Kokkos::CudaUVMSpace
#define KOKKOSBLAS_IMPL_MV_SCALAR double

template<>
struct Axpby<KOKKOSBLAS_IMPL_MV_SCALAR,
             Kokkos::View<const KOKKOSBLAS_IMPL_MV_SCALAR**,
                          Kokkos::LayoutLeft,
                          Kokkos::Device<KOKKOSBLAS_IMPL_MV_EXEC_SPACE, KOKKOSBLAS_IMPL_MV_MEM_SPACE>,
                          Kokkos::MemoryTraits<Kokkos::Unmanaged>,
                          Kokkos::Impl::ViewDefault>,
             KOKKOSBLAS_IMPL_MV_SCALAR,
             Kokkos::View<KOKKOSBLAS_IMPL_MV_SCALAR**,
                          Kokkos::LayoutLeft,
                          Kokkos::Device<KOKKOSBLAS_IMPL_MV_EXEC_SPACE, KOKKOSBLAS_IMPL_MV_MEM_SPACE>,
                          Kokkos::MemoryTraits<Kokkos::Unmanaged>,
                          Kokkos::Impl::ViewDefault>, 2>
{
  typedef KOKKOSBLAS_IMPL_MV_SCALAR AV;
  typedef Kokkos::View<const KOKKOSBLAS_IMPL_MV_SCALAR**,
                       Kokkos::LayoutLeft,
                       Kokkos::Device<KOKKOSBLAS_IMPL_MV_EXEC_SPACE, KOKKOSBLAS_IMPL_MV_MEM_SPACE>,
                       Kokkos::MemoryTraits<Kokkos::Unmanaged>,
                       Kokkos::Impl::ViewDefault> XMV;
  typedef KOKKOSBLAS_IMPL_MV_SCALAR BV;
  typedef Kokkos::View<KOKKOSBLAS_IMPL_MV_SCALAR**,
                       Kokkos::LayoutLeft,
                       Kokkos::Device<KOKKOSBLAS_IMPL_MV_EXEC_SPACE, KOKKOSBLAS_IMPL_MV_MEM_SPACE>,
                       Kokkos::MemoryTraits<Kokkos::Unmanaged>,
                       Kokkos::Impl::ViewDefault> YMV;
  typedef YMV::size_type size_type;
  typedef Kokkos::Details::ArithTraits<XMV::non_const_value_type> ATA;
  typedef Kokkos::Details::ArithTraits<YMV::non_const_value_type> ATB;

  static void
  axpby (const XMV::non_const_value_type& alpha, const XMV& X,
         const YMV::non_const_value_type& beta, const YMV& Y);
};

#undef KOKKOSBLAS_IMPL_MV_SCALAR
#undef KOKKOSBLAS_IMPL_MV_MEM_SPACE
#undef KOKKOSBLAS_IMPL_MV_EXEC_SPACE
#endif // KOKKOS_HAVE_CUDA

} // namespace Impl
} // namespace KokkosBlas

#endif // KOKKOS_BLAS1_MV_IMPL_AXPBY_HPP_
