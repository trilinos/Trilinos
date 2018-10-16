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
#ifndef KOKKOSBLAS1_AXPBY_MV_IMPL_HPP_
#define KOKKOSBLAS1_AXPBY_MV_IMPL_HPP_

#include<KokkosBlas1_axpby_impl.hpp>

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
struct Axpby_MV_Functor
{
  typedef typename YMV::execution_space execution_space;
  typedef SizeType size_type;
  typedef Kokkos::Details::ArithTraits<typename YMV::non_const_value_type> ATS;

  const size_type numCols;
  XMV m_x;
  YMV m_y;
  AV m_a;
  BV m_b;

  Axpby_MV_Functor (const XMV& X, const YMV& Y, const AV& a, const BV& b) :
    numCols (X.extent(1)), m_x (X), m_y (Y), m_a (a), m_b (b)
  {
    // XMV and YMV must be Kokkos::View specializations.
    static_assert (Kokkos::Impl::is_view<AV>::value, "KokkosBlas::Impl::"
                   "Axpby_MV_Functor: a is not a Kokkos::View.");
    static_assert (Kokkos::Impl::is_view<XMV>::value, "KokkosBlas::Impl::"
                   "Axpby_MV_Functor: X is not a Kokkos::View.");
    static_assert (Kokkos::Impl::is_view<BV>::value, "KokkosBlas::Impl::"
                   "Axpby_MV_Functor: b is not a Kokkos::View.");
    static_assert (Kokkos::Impl::is_view<YMV>::value, "KokkosBlas::Impl::"
                   "Axpby_MV_Functor: Y is not a Kokkos::View.");
    // YMV must be nonconst (else it can't be an output argument).
    static_assert (Kokkos::Impl::is_same<typename YMV::value_type,
                   typename YMV::non_const_value_type>::value,
                   "KokkosBlas::Impl::Axpby_MV_Functor: Y is const.  "
                   "It must be nonconst, because it is an output argument "
                   "(we have to be able to write to its entries).");
    static_assert ((int) YMV::Rank == (int) XMV::Rank, "KokkosBlas::Impl::Axpby_MV_Functor: "
                   "X and Y must have the same rank.");
    static_assert (YMV::Rank == 2, "KokkosBlas::Impl::Axpby_MV_Functor: "
                   "XMV and YMV must have rank 2.");
    static_assert (AV::Rank == 1, "KokkosBlas::Impl::Axpby_MV_Functor: "
                   "AV must have rank 1.");
    static_assert (BV::Rank == 1, "KokkosBlas::Impl::Axpby_MV_Functor: "
                   "BV must have rank 1.");
  }

  KOKKOS_INLINE_FUNCTION
  void operator() (const size_type& i) const
  {
    // scalar_x and scalar_y are compile-time constants (since they
    // are template parameters), so the compiler should evaluate these
    // branches at compile time.
    if (scalar_x == 0 && scalar_y == 0) {
#ifdef KOKKOS_ENABLE_PRAGMA_IVDEP
#pragma ivdep
#endif
#ifdef KOKKOS_ENABLE_PRAGMA_VECTOR
#pragma vector always
#endif
      for (size_type k = 0; k < numCols; ++k) {
        m_y(i,k) = ATS::zero ();
      }
    }
    if (scalar_x == 0 && scalar_y == -1) {
#ifdef KOKKOS_ENABLE_PRAGMA_IVDEP
#pragma ivdep
#endif
#ifdef KOKKOS_ENABLE_PRAGMA_VECTOR
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
#ifdef KOKKOS_ENABLE_PRAGMA_IVDEP
#pragma ivdep
#endif
#ifdef KOKKOS_ENABLE_PRAGMA_VECTOR
#pragma vector always
#endif
      for (size_type k = 0; k < numCols; ++k) {
        m_y(i,k) = m_b(k)*m_y(i,k);
      }
    }
    if (scalar_x == -1 && scalar_y == 0) {
#ifdef KOKKOS_ENABLE_PRAGMA_IVDEP
#pragma ivdep
#endif
#ifdef KOKKOS_ENABLE_PRAGMA_VECTOR
#pragma vector always
#endif
      for (size_type k = 0; k < numCols; ++k) {
        m_y(i,k) = -m_x(i,k);
      }
    }
    if (scalar_x == -1 && scalar_y == -1) {
#ifdef KOKKOS_ENABLE_PRAGMA_IVDEP
#pragma ivdep
#endif
#ifdef KOKKOS_ENABLE_PRAGMA_VECTOR
#pragma vector always
#endif
      for (size_type k = 0; k < numCols; ++k) {
        m_y(i,k) = -m_x(i,k) - m_y(i,k);
      }
    }
    if (scalar_x == -1 && scalar_y == 1) {
#ifdef KOKKOS_ENABLE_PRAGMA_IVDEP
#pragma ivdep
#endif
#ifdef KOKKOS_ENABLE_PRAGMA_VECTOR
#pragma vector always
#endif
      for (size_type k = 0; k < numCols; ++k) {
        m_y(i,k) = -m_x(i,k) + m_y(i,k);
      }
    }
    if (scalar_x == -1 && scalar_y == 2) {
#ifdef KOKKOS_ENABLE_PRAGMA_IVDEP
#pragma ivdep
#endif
#ifdef KOKKOS_ENABLE_PRAGMA_VECTOR
#pragma vector always
#endif
      for (size_type k = 0; k < numCols; ++k) {
        m_y(i,k) = -m_x(i,k) + m_b(k)*m_y(i,k);
      }
    }
    if (scalar_x == 1 && scalar_y == 0) {
#ifdef KOKKOS_ENABLE_PRAGMA_IVDEP
#pragma ivdep
#endif
#ifdef KOKKOS_ENABLE_PRAGMA_VECTOR
#pragma vector always
#endif
      for (size_type k = 0; k < numCols; ++k) {
        m_y(i,k) = m_x(i,k);
      }
    }
    if (scalar_x == 1 && scalar_y == -1) {
#ifdef KOKKOS_ENABLE_PRAGMA_IVDEP
#pragma ivdep
#endif
#ifdef KOKKOS_ENABLE_PRAGMA_VECTOR
#pragma vector always
#endif
      for (size_type k = 0; k < numCols; ++k) {
        m_y(i,k) = m_x(i,k) - m_y(i,k);
      }
    }
    if (scalar_x == 1 && scalar_y == 1) {
#ifdef KOKKOS_ENABLE_PRAGMA_IVDEP
#pragma ivdep
#endif
#ifdef KOKKOS_ENABLE_PRAGMA_VECTOR
#pragma vector always
#endif
      for (size_type k = 0; k < numCols; ++k) {
        m_y(i,k) = m_x(i,k) + m_y(i,k);
      }
    }
    if (scalar_x == 1 && scalar_y == 2) {
#ifdef KOKKOS_ENABLE_PRAGMA_IVDEP
#pragma ivdep
#endif
#ifdef KOKKOS_ENABLE_PRAGMA_VECTOR
#pragma vector always
#endif
      for (size_type k = 0; k < numCols; ++k) {
        m_y(i,k) = m_x(i,k) + m_b(k)*m_y(i,k);
      }
    }
    if (scalar_x == 2 && scalar_y == 0) {
#ifdef KOKKOS_ENABLE_PRAGMA_IVDEP
#pragma ivdep
#endif
#ifdef KOKKOS_ENABLE_PRAGMA_VECTOR
#pragma vector always
#endif
      for (size_type k = 0; k < numCols; ++k) {
        m_y(i,k) = m_a(k)*m_x(i,k);
      }
    }
    if (scalar_x == 2 && scalar_y == -1) {
#ifdef KOKKOS_ENABLE_PRAGMA_IVDEP
#pragma ivdep
#endif
#ifdef KOKKOS_ENABLE_PRAGMA_VECTOR
#pragma vector always
#endif
      for (size_type k = 0; k < numCols; ++k) {
        m_y(i,k) = m_a(k)*m_x(i,k) - m_y(i,k);
      }
    }
    if (scalar_x == 2 && scalar_y == 1) {
#ifdef KOKKOS_ENABLE_PRAGMA_IVDEP
#pragma ivdep
#endif
#ifdef KOKKOS_ENABLE_PRAGMA_VECTOR
#pragma vector always
#endif
      for (size_type k = 0; k < numCols; ++k) {
        m_y(i,k) = m_a(k)*m_x(i,k) + m_y(i,k);
      }
    }
    if (scalar_x == 2 && scalar_y == 2) {
#ifdef KOKKOS_ENABLE_PRAGMA_IVDEP
#pragma ivdep
#endif
#ifdef KOKKOS_ENABLE_PRAGMA_VECTOR
#pragma vector always
#endif
      for (size_type k = 0; k < numCols; ++k) {
        m_y(i,k) = m_a(k)*m_x(i,k) + m_b(k)*m_y(i,k);
      }
    }
  }
};

// Variant of Axpby_MV_Functor, where a and b are scalars.
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
struct Axpby_MV_Functor<typename XMV::non_const_value_type, XMV,
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

  Axpby_MV_Functor (const XMV& X, const YMV& Y,
                    const typename XMV::non_const_value_type& a,
                    const typename YMV::non_const_value_type& b) :
    numCols (X.extent(1)), m_x (X), m_y (Y), m_a (a), m_b (b)
  {
    static_assert (Kokkos::Impl::is_view<XMV>::value, "KokkosBlas::Impl::"
                   "Axpby_MV_Functor: X is not a Kokkos::View.");
    static_assert (Kokkos::Impl::is_view<YMV>::value, "KokkosBlas::Impl::"
                   "Axpby_MV_Functor: Y is not a Kokkos::View.");
    static_assert (Kokkos::Impl::is_same<typename YMV::value_type,
                   typename YMV::non_const_value_type>::value,
                   "KokkosBlas::Impl::Axpby_MV_Functor: Y is const.  "
                   "It must be nonconst, because it is an output argument "
                   "(we have to be able to write to its entries).");
    static_assert ((int) YMV::Rank == (int) XMV::Rank, "KokkosBlas::Impl::"
                   "Axpby_MV_Functor: X and Y must have the same rank.");
    static_assert (YMV::Rank == 2, "KokkosBlas::Impl::Axpby_MV_Functor: "
                   "XMV and YMV must have rank 2.");
  }

  KOKKOS_INLINE_FUNCTION
  void operator() (const size_type& i) const
  {
    // scalar_x and scalar_y are compile-time constants (since they
    // are template parameters), so the compiler should evaluate these
    // branches at compile time.
    if (scalar_x == 0 && scalar_y == 0) {
#ifdef KOKKOS_ENABLE_PRAGMA_IVDEP
#pragma ivdep
#endif
#ifdef KOKKOS_ENABLE_PRAGMA_VECTOR
#pragma vector always
#endif
      for (size_type k = 0; k < numCols; ++k) {
        m_y(i,k) = ATS::zero ();
      }
    }
    if (scalar_x == 0 && scalar_y == -1) {
#ifdef KOKKOS_ENABLE_PRAGMA_IVDEP
#pragma ivdep
#endif
#ifdef KOKKOS_ENABLE_PRAGMA_VECTOR
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
#ifdef KOKKOS_ENABLE_PRAGMA_IVDEP
#pragma ivdep
#endif
#ifdef KOKKOS_ENABLE_PRAGMA_VECTOR
#pragma vector always
#endif
      for (size_type k = 0; k < numCols; ++k) {
        m_y(i,k) = m_b*m_y(i,k);
      }
    }
    if (scalar_x == -1 && scalar_y == 0) {
#ifdef KOKKOS_ENABLE_PRAGMA_IVDEP
#pragma ivdep
#endif
#ifdef KOKKOS_ENABLE_PRAGMA_VECTOR
#pragma vector always
#endif
      for (size_type k = 0; k < numCols; ++k) {
        m_y(i,k) = -m_x(i,k);
      }
    }
    if (scalar_x == -1 && scalar_y == -1) {
#ifdef KOKKOS_ENABLE_PRAGMA_IVDEP
#pragma ivdep
#endif
#ifdef KOKKOS_ENABLE_PRAGMA_VECTOR
#pragma vector always
#endif
      for (size_type k = 0; k < numCols; ++k) {
        m_y(i,k) = -m_x(i,k) - m_y(i,k);
      }
    }
    if (scalar_x == -1 && scalar_y == 1) {
#ifdef KOKKOS_ENABLE_PRAGMA_IVDEP
#pragma ivdep
#endif
#ifdef KOKKOS_ENABLE_PRAGMA_VECTOR
#pragma vector always
#endif
      for (size_type k = 0; k < numCols; ++k) {
        m_y(i,k) = -m_x(i,k) + m_y(i,k);
      }
    }
    if (scalar_x == -1 && scalar_y == 2) {
#ifdef KOKKOS_ENABLE_PRAGMA_IVDEP
#pragma ivdep
#endif
#ifdef KOKKOS_ENABLE_PRAGMA_VECTOR
#pragma vector always
#endif
      for (size_type k = 0; k < numCols; ++k) {
        m_y(i,k) = -m_x(i,k) + m_b*m_y(i,k);
      }
    }
    if (scalar_x == 1 && scalar_y == 0) {
#ifdef KOKKOS_ENABLE_PRAGMA_IVDEP
#pragma ivdep
#endif
#ifdef KOKKOS_ENABLE_PRAGMA_VECTOR
#pragma vector always
#endif
      for (size_type k = 0; k < numCols; ++k) {
        m_y(i,k) = m_x(i,k);
      }
    }
    if (scalar_x == 1 && scalar_y == -1) {
#ifdef KOKKOS_ENABLE_PRAGMA_IVDEP
#pragma ivdep
#endif
#ifdef KOKKOS_ENABLE_PRAGMA_VECTOR
#pragma vector always
#endif
      for (size_type k = 0; k < numCols; ++k) {
        m_y(i,k) = m_x(i,k) - m_y(i,k);
      }
    }
    if (scalar_x == 1 && scalar_y == 1) {
#ifdef KOKKOS_ENABLE_PRAGMA_IVDEP
#pragma ivdep
#endif
#ifdef KOKKOS_ENABLE_PRAGMA_VECTOR
#pragma vector always
#endif
      for (size_type k = 0; k < numCols; ++k) {
        m_y(i,k) = m_x(i,k) + m_y(i,k);
      }
    }
    if (scalar_x == 1 && scalar_y == 2) {
#ifdef KOKKOS_ENABLE_PRAGMA_IVDEP
#pragma ivdep
#endif
#ifdef KOKKOS_ENABLE_PRAGMA_VECTOR
#pragma vector always
#endif
      for (size_type k = 0; k < numCols; ++k) {
        m_y(i,k) = m_x(i,k) + m_b*m_y(i,k);
      }
    }
    if (scalar_x == 2 && scalar_y == 0) {
#ifdef KOKKOS_ENABLE_PRAGMA_IVDEP
#pragma ivdep
#endif
#ifdef KOKKOS_ENABLE_PRAGMA_VECTOR
#pragma vector always
#endif
      for (size_type k = 0; k < numCols; ++k) {
        m_y(i,k) = m_a*m_x(i,k);
      }
    }
    if (scalar_x == 2 && scalar_y == -1) {
#ifdef KOKKOS_ENABLE_PRAGMA_IVDEP
#pragma ivdep
#endif
#ifdef KOKKOS_ENABLE_PRAGMA_VECTOR
#pragma vector always
#endif
      for (size_type k = 0; k < numCols; ++k) {
        m_y(i,k) = m_a*m_x(i,k) - m_y(i,k);
      }
    }
    if (scalar_x == 2 && scalar_y == 1) {
#ifdef KOKKOS_ENABLE_PRAGMA_IVDEP
#pragma ivdep
#endif
#ifdef KOKKOS_ENABLE_PRAGMA_VECTOR
#pragma vector always
#endif
      for (size_type k = 0; k < numCols; ++k) {
        m_y(i,k) = m_a*m_x(i,k) + m_y(i,k);
      }
    }
    if (scalar_x == 2 && scalar_y == 2) {
#ifdef KOKKOS_ENABLE_PRAGMA_IVDEP
#pragma ivdep
#endif
#ifdef KOKKOS_ENABLE_PRAGMA_VECTOR
#pragma vector always
#endif
      for (size_type k = 0; k < numCols; ++k) {
        m_y(i,k) = m_a*m_x(i,k) + m_b*m_y(i,k);
      }
    }
  }
};


// Column-unrolled variant of Axpby_MV_Functor.  The number of columns
// in X and Y, UNROLL, is a compile-time constant.
template<class AV, class XMV, class BV, class YMV,
         int scalar_x, int scalar_y, int UNROLL, class SizeType>
struct Axpby_MV_Unroll_Functor
{
  typedef typename YMV::execution_space execution_space;
  typedef SizeType size_type;
  typedef Kokkos::Details::ArithTraits<typename YMV::non_const_value_type> ATS;

  XMV m_x;
  YMV m_y;
  AV m_a;
  BV m_b;

  Axpby_MV_Unroll_Functor (const XMV& x, const YMV& y,
                           const AV& a, const BV& b,
                           const SizeType startingColumn) :
    m_x (x), m_y (y), m_a (a), m_b (b)
  {
    static_assert (Kokkos::Impl::is_view<AV>::value, "KokkosBlas::Impl::"
                   "Axpby_MV_Unroll_Functor: a is not a Kokkos::View.");
    static_assert (Kokkos::Impl::is_view<XMV>::value, "KokkosBlas::Impl::"
                   "Axpby_MV_Unroll_Functor: X is not a Kokkos::View.");
    static_assert (Kokkos::Impl::is_view<BV>::value, "KokkosBlas::Impl::"
                   "Axpby_MV_Unroll_Functor: b is not a Kokkos::View.");
    static_assert (Kokkos::Impl::is_view<YMV>::value, "KokkosBlas::Impl::"
                   "Axpby_MV_Unroll_Functor: Y is not a Kokkos::View.");
    static_assert (Kokkos::Impl::is_same<typename YMV::value_type,
                   typename YMV::non_const_value_type>::value,
                   "KokkosBlas::Impl::Axpby_MV_Unroll_Functor: Y is const.  "
                   "It must be nonconst, because it is an output argument "
                   "(we have to be able to write to its entries).");
    static_assert ((int) YMV::Rank == (int) XMV::Rank,
                   "KokkosBlas::Impl::Axpby_MV_Unroll_Functor: "
                   "X and Y must have the same rank.");
    static_assert (YMV::Rank == 2, "KokkosBlas::Impl::Axpby_MV_Unroll_Functor: "
                   "XMV and YMV must have rank 2.");
    static_assert (AV::Rank == 1, "KokkosBlas::Impl::Axpby_MV_Unroll_Functor: "
                   "AV must have rank 1.");
    static_assert (BV::Rank == 1, "KokkosBlas::Impl::Axpby_MV_Unroll_Functor: "
                   "BV must have rank 1.");

    if (startingColumn != 0) {
      m_a = Kokkos::subview (a, std::make_pair (startingColumn, SizeType(a.extent(0))));
      m_b = Kokkos::subview (b, std::make_pair (startingColumn, SizeType(b.extent(0))));
    }
  }

  KOKKOS_INLINE_FUNCTION
  void operator() (const size_type& i) const
  {
    // scalar_x and scalar_y are compile-time constants (since they
    // are template parameters), so the compiler should evaluate these
    // branches at compile time.

#if KOKKOSBLAS_OPTIMIZATION_LEVEL_AXPBY <= 2

    if (scalar_x == 0 && scalar_y == 0) {
#ifdef KOKKOS_ENABLE_PRAGMA_UNROLL
#pragma unroll
#endif
      for (int k = 0; k < UNROLL; ++k) {
        m_y(i,k) = ATS::zero ();
      }
    }
    if (scalar_x == 0 && scalar_y == 2) {
#ifdef KOKKOS_ENABLE_PRAGMA_UNROLL
#pragma unroll
#endif
      for (int k = 0; k < UNROLL; ++k) {
        m_y(i,k) = m_b(k)*m_y(i,k);
      }
    }
    if (scalar_x == 2 && scalar_y == 0) {
#ifdef KOKKOS_ENABLE_PRAGMA_UNROLL
#pragma unroll
#endif
      for (int k = 0; k < UNROLL; ++k) {
        m_y(i,k) = m_a(k)*m_x(i,k);
      }
    }
    if (scalar_x == 2 && scalar_y == 2) {
#ifdef KOKKOS_ENABLE_PRAGMA_UNROLL
#pragma unroll
#endif
      for (int k = 0; k < UNROLL; ++k) {
        m_y(i,k) = m_a(k)*m_x(i,k) + m_b(k)*m_y(i,k);
      }
    }

#else // KOKKOSBLAS_OPTIMIZATION_LEVEL >= 3

    if (scalar_x == 0 && scalar_y == 0) {
#ifdef KOKKOS_ENABLE_PRAGMA_UNROLL
#pragma unroll
#endif
      for (int k = 0; k < UNROLL; ++k) {
        m_y(i,k) = ATS::zero ();
      }
    }
    if (scalar_x == 0 && scalar_y == -1) {
#ifdef KOKKOS_ENABLE_PRAGMA_UNROLL
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
#ifdef KOKKOS_ENABLE_PRAGMA_UNROLL
#pragma unroll
#endif
      for (int k = 0; k < UNROLL; ++k) {
        m_y(i,k) = m_b(k)*m_y(i,k);
      }
    }
    if (scalar_x == -1 && scalar_y == 0) {
#ifdef KOKKOS_ENABLE_PRAGMA_UNROLL
#pragma unroll
#endif
      for (int k = 0; k < UNROLL; ++k) {
        m_y(i,k) = -m_x(i,k);
      }
    }
    if (scalar_x == -1 && scalar_y == -1) {
#ifdef KOKKOS_ENABLE_PRAGMA_UNROLL
#pragma unroll
#endif
      for (int k = 0; k < UNROLL; ++k) {
        m_y(i,k) = -m_x(i,k) - m_y(i,k);
      }
    }
    if (scalar_x == -1 && scalar_y == 1) {
#ifdef KOKKOS_ENABLE_PRAGMA_UNROLL
#pragma unroll
#endif
      for (int k = 0; k < UNROLL; ++k) {
        m_y(i,k) = -m_x(i,k) + m_y(i,k);
      }
    }
    if (scalar_x == -1 && scalar_y == 2) {
#ifdef KOKKOS_ENABLE_PRAGMA_UNROLL
#pragma unroll
#endif
      for (int k = 0; k < UNROLL; ++k) {
        m_y(i,k) = -m_x(i,k) + m_b(k)*m_y(i,k);
      }
    }
    if (scalar_x == 1 && scalar_y == 0) {
#ifdef KOKKOS_ENABLE_PRAGMA_UNROLL
#pragma unroll
#endif
      for (int k = 0; k < UNROLL; ++k) {
        m_y(i,k) = m_x(i,k);
      }
    }
    if (scalar_x == 1 && scalar_y == -1) {
#ifdef KOKKOS_ENABLE_PRAGMA_UNROLL
#pragma unroll
#endif
      for (int k = 0; k < UNROLL; ++k) {
        m_y(i,k) = m_x(i,k) - m_y(i,k);
      }
    }
    if (scalar_x == 1 && scalar_y == 1) {
#ifdef KOKKOS_ENABLE_PRAGMA_UNROLL
#pragma unroll
#endif
      for (int k = 0; k < UNROLL; ++k) {
        m_y(i,k) = m_x(i,k) + m_y(i,k);
      }
    }
    if (scalar_x == 1 && scalar_y == 2) {
#ifdef KOKKOS_ENABLE_PRAGMA_UNROLL
#pragma unroll
#endif
      for (int k = 0; k < UNROLL; ++k) {
        m_y(i,k) = m_x(i,k) + m_b(k)*m_y(i,k);
      }
    }
    if (scalar_x == 2 && scalar_y == 0) {
#ifdef KOKKOS_ENABLE_PRAGMA_UNROLL
#pragma unroll
#endif
      for (int k = 0; k < UNROLL; ++k) {
        m_y(i,k) = m_a(k)*m_x(i,k);
      }
    }
    if (scalar_x == 2 && scalar_y == -1) {
#ifdef KOKKOS_ENABLE_PRAGMA_UNROLL
#pragma unroll
#endif
      for (int k = 0; k < UNROLL; ++k) {
        m_y(i,k) = m_a(k)*m_x(i,k) - m_y(i,k);
      }
    }
    if (scalar_x == 2 && scalar_y == 1) {
#ifdef KOKKOS_ENABLE_PRAGMA_UNROLL
#pragma unroll
#endif
      for (int k = 0; k < UNROLL; ++k) {
        m_y(i,k) = m_a(k)*m_x(i,k) + m_y(i,k);
      }
    }
    if (scalar_x == 2 && scalar_y == 2) {
#ifdef KOKKOS_ENABLE_PRAGMA_UNROLL
#pragma unroll
#endif
      for (int k = 0; k < UNROLL; ++k) {
        m_y(i,k) = m_a(k)*m_x(i,k) + m_b(k)*m_y(i,k);
      }
    }
#endif // KOKKOSBLAS_OPTIMIZATION_LEVEL_AXPBY
  }
};


// Variant of Axpby_MV_Unroll_Functor for single coefficients (rather
// than vectors of coefficients) a and b.  The number of columns in X
// and Y, UNROLL, is a compile-time constant.
template<class XMV, class YMV, int scalar_x, int scalar_y,
         int UNROLL, class SizeType>
struct Axpby_MV_Unroll_Functor<typename XMV::non_const_value_type, XMV,
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

  Axpby_MV_Unroll_Functor (const XMV& X, const YMV& Y,
                           const typename XMV::non_const_value_type& a,
                           const typename YMV::non_const_value_type& b,
                           const SizeType /* startingColumn */) :
    m_x (X), m_y (Y), m_a (a), m_b (b)
  {
    static_assert (Kokkos::Impl::is_view<XMV>::value, "KokkosBlas::Impl::"
                   "Axpby_MV_Unroll_Functor: X is not a Kokkos::View.");
    static_assert (Kokkos::Impl::is_view<YMV>::value, "KokkosBlas::Impl::"
                   "Axpby_MV_Unroll_Functor: Y is not a Kokkos::View.");
    static_assert (Kokkos::Impl::is_same<typename YMV::value_type,
                   typename YMV::non_const_value_type>::value,
                   "KokkosBlas::Impl::Axpby_MV_Unroll_Functor: Y is const.  "
                   "It must be nonconst, because it is an output argument "
                   "(we have to be able to write to its entries).");
    static_assert ((int) YMV::Rank == (int) XMV::Rank, "KokkosBlas::Impl::"
                   "Axpby_MV_Unroll_Functor: X and Y must have the same rank.");
    static_assert (YMV::Rank == 2, "KokkosBlas::Impl::Axpby_MV_Unroll_Functor: "
                   "XMV and YMV must have rank 2.");
  }

  KOKKOS_INLINE_FUNCTION
  void operator() (const size_type& i) const
  {
    // scalar_x and scalar_y are compile-time constants (since they
    // are template parameters), so the compiler should evaluate these
    // branches at compile time.

#if KOKKOSBLAS_OPTIMIZATION_LEVEL_AXPBY <= 2

    if (scalar_x == 0 && scalar_y == 0) {
#ifdef KOKKOS_ENABLE_PRAGMA_UNROLL
#pragma unroll
#endif
      for (int k = 0; k < UNROLL; ++k) {
        m_y(i,k) = ATS::zero ();
      }
    }
    if (scalar_x == 0 && scalar_y == 2) {
#ifdef KOKKOS_ENABLE_PRAGMA_UNROLL
#pragma unroll
#endif
      for (int k = 0; k < UNROLL; ++k) {
        m_y(i,k) = m_b*m_y(i,k);
      }
    }
    if (scalar_x == 2 && scalar_y == 0) {
#ifdef KOKKOS_ENABLE_PRAGMA_UNROLL
#pragma unroll
#endif
      for (int k = 0; k < UNROLL; ++k) {
        m_y(i,k) = m_a*m_x(i,k);
      }
    }
    if (scalar_x == 2 && scalar_y == 2) {
#ifdef KOKKOS_ENABLE_PRAGMA_UNROLL
#pragma unroll
#endif
      for (int k = 0; k < UNROLL; ++k) {
        m_y(i,k) = m_a*m_x(i,k) + m_b*m_y(i,k);
      }
    }

#else // KOKKOSBLAS_OPTIMIZATION_LEVEL_AXPBY > 2

    if (scalar_x == 0 && scalar_y == 0) {
#ifdef KOKKOS_ENABLE_PRAGMA_UNROLL
#pragma unroll
#endif
      for (int k = 0; k < UNROLL; ++k) {
        m_y(i,k) = ATS::zero ();
      }
    }
    if (scalar_x == 0 && scalar_y == -1) {
#ifdef KOKKOS_ENABLE_PRAGMA_UNROLL
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
#ifdef KOKKOS_ENABLE_PRAGMA_UNROLL
#pragma unroll
#endif
      for (int k = 0; k < UNROLL; ++k) {
        m_y(i,k) = m_b*m_y(i,k);
      }
    }
    if (scalar_x == -1 && scalar_y == 0) {
#ifdef KOKKOS_ENABLE_PRAGMA_UNROLL
#pragma unroll
#endif
      for (int k = 0; k < UNROLL; ++k) {
        m_y(i,k) = -m_x(i,k);
      }
    }
    if (scalar_x == -1 && scalar_y == -1) {
#ifdef KOKKOS_ENABLE_PRAGMA_UNROLL
#pragma unroll
#endif
      for (int k = 0; k < UNROLL; ++k) {
        m_y(i,k) = -m_x(i,k) - m_y(i,k);
      }
    }
    if (scalar_x == -1 && scalar_y == 1) {
#ifdef KOKKOS_ENABLE_PRAGMA_UNROLL
#pragma unroll
#endif
      for (int k = 0; k < UNROLL; ++k) {
        m_y(i,k) = -m_x(i,k) + m_y(i,k);
      }
    }
    if (scalar_x == -1 && scalar_y == 2) {
#ifdef KOKKOS_ENABLE_PRAGMA_UNROLL
#pragma unroll
#endif
      for (int k = 0; k < UNROLL; ++k) {
        m_y(i,k) = -m_x(i,k) + m_b*m_y(i,k);
      }
    }
    if (scalar_x == 1 && scalar_y == 0) {
#ifdef KOKKOS_ENABLE_PRAGMA_UNROLL
#pragma unroll
#endif
      for (int k = 0; k < UNROLL; ++k) {
        m_y(i,k) = m_x(i,k);
      }
    }
    if (scalar_x == 1 && scalar_y == -1) {
#ifdef KOKKOS_ENABLE_PRAGMA_UNROLL
#pragma unroll
#endif
      for (int k = 0; k < UNROLL; ++k) {
        m_y(i,k) = m_x(i,k) - m_y(i,k);
      }
    }
    if (scalar_x == 1 && scalar_y == 1) {
#ifdef KOKKOS_ENABLE_PRAGMA_UNROLL
#pragma unroll
#endif
      for (int k = 0; k < UNROLL; ++k) {
        m_y(i,k) = m_x(i,k) + m_y(i,k);
      }
    }
    if (scalar_x == 1 && scalar_y == 2) {
#ifdef KOKKOS_ENABLE_PRAGMA_UNROLL
#pragma unroll
#endif
      for (int k = 0; k < UNROLL; ++k) {
        m_y(i,k) = m_x(i,k) + m_b*m_y(i,k);
      }
    }
    if (scalar_x == 2 && scalar_y == 0) {
#ifdef KOKKOS_ENABLE_PRAGMA_UNROLL
#pragma unroll
#endif
      for (int k = 0; k < UNROLL; ++k) {
        m_y(i,k) = m_a*m_x(i,k);
      }
    }
    if (scalar_x == 2 && scalar_y == -1) {
#ifdef KOKKOS_ENABLE_PRAGMA_UNROLL
#pragma unroll
#endif
      for (int k = 0; k < UNROLL; ++k) {
        m_y(i,k) = m_a*m_x(i,k) - m_y(i,k);
      }
    }
    if (scalar_x == 2 && scalar_y == 1) {
#ifdef KOKKOS_ENABLE_PRAGMA_UNROLL
#pragma unroll
#endif
      for (int k = 0; k < UNROLL; ++k) {
        m_y(i,k) = m_a*m_x(i,k) + m_y(i,k);
      }
    }
    if (scalar_x == 2 && scalar_y == 2) {
#ifdef KOKKOS_ENABLE_PRAGMA_UNROLL
#pragma unroll
#endif
      for (int k = 0; k < UNROLL; ++k) {
        m_y(i,k) = m_a*m_x(i,k) + m_b*m_y(i,k);
      }
    }

#endif // KOKKOSBLAS_OPTIMIZATION_LEVEL_AXPBY
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
Axpby_MV_Unrolled (const AV& av, const XMV& x,
                   const BV& bv, const YMV& y,
                   const SizeType startingColumn,
                   int a = 2, int b = 2)
{
  static_assert (Kokkos::Impl::is_view<XMV>::value, "KokkosBlas::Impl::"
                 "Axpby_MV_Unrolled: X is not a Kokkos::View.");
  static_assert (Kokkos::Impl::is_view<YMV>::value, "KokkosBlas::Impl::"
                 "Axpby_MV_Unrolled: Y is not a Kokkos::View.");
  static_assert (Kokkos::Impl::is_same<typename YMV::value_type,
                   typename YMV::non_const_value_type>::value,
                 "KokkosBlas::Impl::Axpby_MV_Unrolled: Y is const.  "
                 "It must be nonconst, because it is an output argument "
                 "(we have to be able to write to its entries).");
  static_assert ((int) YMV::Rank == (int) XMV::Rank, "KokkosBlas::Impl::"
                 "Axpby_MV_Unrolled: X and Y must have the same rank.");
  static_assert (YMV::Rank == 2, "KokkosBlas::Impl::Axpby_MV_Unrolled: "
                 "XMV and YMV must have rank 2.");

  typedef typename YMV::execution_space execution_space;
  const SizeType numRows = x.extent(0);
  Kokkos::RangePolicy<execution_space, SizeType> policy (0, numRows);

  if (a == 0 && b == 0) {
    Axpby_MV_Unroll_Functor<AV, XMV, BV, YMV, 0, 0, UNROLL, SizeType> op (x, y, av, bv, startingColumn);
    Kokkos::parallel_for (policy, op);
    return;
  }

#if KOKKOSBLAS_OPTIMIZATION_LEVEL_AXPBY > 2
  if (a == 0 && b == -1) {
    Axpby_MV_Unroll_Functor<AV, XMV, BV, YMV, 0, -1, UNROLL, SizeType> op (x, y, av, bv, startingColumn);
    Kokkos::parallel_for (policy, op);
    return;
  }
  if (a == 0 && b == 1) {
    Axpby_MV_Unroll_Functor<AV, XMV, BV, YMV, 0, 1, UNROLL, SizeType> op (x, y, av, bv, startingColumn);
    Kokkos::parallel_for (policy, op);
    return;
  }
#endif // KOKKOSBLAS_OPTIMIZATION_LEVEL_AXPBY

  if (a == 0 && b == 2) {
    Axpby_MV_Unroll_Functor<AV, XMV, BV, YMV, 0, 2, UNROLL, SizeType> op (x, y, av, bv, startingColumn);
    Kokkos::parallel_for (policy, op);
    return;
  }

#if KOKKOSBLAS_OPTIMIZATION_LEVEL_AXPBY > 2
  // a == -1
  if (a == -1 && b == 0) {
    Axpby_MV_Unroll_Functor<AV, XMV, BV, YMV, -1, 0, UNROLL, SizeType> op (x, y, av, bv, startingColumn);
    Kokkos::parallel_for (policy, op);
    return;
  }
  if (a == -1 && b == -1) {
    Axpby_MV_Unroll_Functor<AV, XMV, BV, YMV, -1, -1, UNROLL, SizeType> op (x, y, av, bv, startingColumn);
    Kokkos::parallel_for (policy, op);
    return;
  }
  if (a == -1 && b == 1) {
    Axpby_MV_Unroll_Functor<AV, XMV, BV, YMV, -1, 1, UNROLL, SizeType> op (x, y, av, bv, startingColumn);
    Kokkos::parallel_for (policy, op);
    return;
  }
  if (a == -1 && b == 2) {
    Axpby_MV_Unroll_Functor<AV, XMV, BV, YMV, -1, 2, UNROLL, SizeType> op (x, y, av, bv, startingColumn);
    Kokkos::parallel_for (policy, op);
    return;
  }
  // a == 1
  if (a == 1 && b == 0) {
    Axpby_MV_Unroll_Functor<AV, XMV, BV, YMV, 1, 0, UNROLL, SizeType> op (x, y, av, bv, startingColumn);
    Kokkos::parallel_for (policy, op);
    return;
  }
  if (a == 1 && b == -1) {
    Axpby_MV_Unroll_Functor<AV, XMV, BV, YMV, 1, -1, UNROLL, SizeType> op (x, y, av, bv, startingColumn);
    Kokkos::parallel_for (policy, op);
    return;
  }
  if (a == 1 && b == 1) {
    Axpby_MV_Unroll_Functor<AV, XMV, BV, YMV, 1, 1, UNROLL, SizeType> op (x, y, av, bv, startingColumn);
    Kokkos::parallel_for (policy, op);
    return;
  }
  if (a == 1 && b == 2) {
    Axpby_MV_Unroll_Functor<AV, XMV, BV, YMV, 1, 2, UNROLL, SizeType> op (x, y, av, bv, startingColumn);
    Kokkos::parallel_for (policy, op);
    return;
  }
#endif // KOKKOSBLAS_OPTIMIZATION_LEVEL_AXPBY > 2

  // a == 2
  if (a == 2 && b == 0) {
    Axpby_MV_Unroll_Functor<AV, XMV, BV, YMV, 2, 0, UNROLL, SizeType> op (x, y, av, bv, startingColumn);
    Kokkos::parallel_for (policy, op);
    return;
  }

#if KOKKOSBLAS_OPTIMIZATION_LEVEL_AXPBY > 2
  if (a == 2 && b == -1) {
    Axpby_MV_Unroll_Functor<AV, XMV, BV, YMV, 2, -1, UNROLL, SizeType> op (x, y, av, bv, startingColumn);
    Kokkos::parallel_for (policy, op);
    return;
  }
  if (a == 2 && b == 1) {
    Axpby_MV_Unroll_Functor<AV, XMV, BV, YMV, 2, 1, UNROLL, SizeType> op (x, y, av, bv, startingColumn);
    Kokkos::parallel_for (policy, op);
    return;
  }
#endif // KOKKOSBLAS_OPTIMIZATION_LEVEL_AXPBY > 2

  // a and b arbitrary (not -1, 0, or 1)
  Axpby_MV_Unroll_Functor<AV, XMV, BV, YMV, 2, 2, UNROLL, SizeType> op (x, y, av, bv, startingColumn);
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
Axpby_MV_Generic (const AV& av, const XMV& x,
                  const BV& bv, const YMV& y,
                  int a = 2, int b = 2)
{
  static_assert (Kokkos::Impl::is_view<XMV>::value, "KokkosBlas::Impl::"
                 "Axpby_MV_Generic: X is not a Kokkos::View.");
  static_assert (Kokkos::Impl::is_view<YMV>::value, "KokkosBlas::Impl::"
                 "Axpby_MV_Generic: Y is not a Kokkos::View.");
  static_assert (Kokkos::Impl::is_same<typename YMV::value_type,
                   typename YMV::non_const_value_type>::value,
                 "KokkosBlas::Impl::Axpby_MV_Generic: Y is const.  "
                 "It must be nonconst, because it is an output argument "
                 "(we have to be able to write to its entries).");
  static_assert ((int) YMV::Rank == (int) XMV::Rank, "KokkosBlas::Impl::"
                 "Axpby_MV_Generic: X and Y must have the same rank.");
  static_assert (YMV::Rank == 2, "KokkosBlas::Impl::Axpby_MV_Generic: "
                 "XMV and YMV must have rank 2.");

  typedef typename YMV::execution_space execution_space;
  const SizeType numRows = x.extent(0);
  Kokkos::RangePolicy<execution_space, SizeType> policy (0, numRows);

  if (a == 0 && b == 0) {
    Axpby_MV_Functor<AV, XMV, BV, YMV, 0, 0, SizeType> op (x, y, av, bv);
    Kokkos::parallel_for (policy, op);
    return;
  }

#if KOKKOSBLAS_OPTIMIZATION_LEVEL_AXPBY > 2
  if (a == 0 && b == -1) {
    Axpby_MV_Functor<AV, XMV, BV, YMV, 0, -1, SizeType> op (x, y, av, bv);
    Kokkos::parallel_for (policy, op);
    return;
  }
  if (a == 0 && b == 1) {
    Axpby_MV_Functor<AV, XMV, BV, YMV, 0, 1, SizeType> op (x, y, av, bv);
    Kokkos::parallel_for (policy, op);
    return;
  }
#endif // KOKKOSBLAS_OPTIMIZATION_LEVEL_AXPBY > 2

  if (a == 0 && b == 2) {
    Axpby_MV_Functor<AV, XMV, BV, YMV, 0, 2, SizeType> op (x, y, av, bv);
    Kokkos::parallel_for (policy, op);
    return;
  }

#if KOKKOSBLAS_OPTIMIZATION_LEVEL_AXPBY > 2
  // a == -1
  if (a == -1 && b == 0) {
    Axpby_MV_Functor<AV, XMV, BV, YMV, -1, 0, SizeType> op (x, y, av, bv);
    Kokkos::parallel_for (policy, op);
    return;
  }
  if (a == -1 && b == -1) {
    Axpby_MV_Functor<AV, XMV, BV, YMV, -1, -1, SizeType> op (x, y, av, bv);
    Kokkos::parallel_for (policy, op);
    return;
  }
  if (a == -1 && b == 1) {
    Axpby_MV_Functor<AV, XMV, BV, YMV, -1, 1, SizeType> op (x, y, av, bv);
    Kokkos::parallel_for (policy, op);
    return;
  }
  if (a == -1 && b == 2) {
    Axpby_MV_Functor<AV, XMV, BV, YMV, -1, 2, SizeType> op (x, y, av, bv);
    Kokkos::parallel_for (policy, op);
    return;
  }
  // a == 1
  if (a == 1 && b == 0) {
    Axpby_MV_Functor<AV, XMV, BV, YMV, 1, 0, SizeType> op (x, y, av, bv);
    Kokkos::parallel_for (policy, op);
    return;
  }
  if (a == 1 && b == -1) {
    Axpby_MV_Functor<AV, XMV, BV, YMV, 1, -1, SizeType> op (x, y, av, bv);
    Kokkos::parallel_for (policy, op);
    return;
  }
  if (a == 1 && b == 1) {
    Axpby_MV_Functor<AV, XMV, BV, YMV, 1, 1, SizeType> op (x, y, av, bv);
    Kokkos::parallel_for (policy, op);
    return;
  }
  if (a == 1 && b == 2) {
    Axpby_MV_Functor<AV, XMV, BV, YMV, 1, 2, SizeType> op (x, y, av, bv);
    Kokkos::parallel_for (policy, op);
    return;
  }
#endif // KOKKOSBLAS_OPTIMIZATION_LEVEL_AXPBY > 2

  // a == 2
  if (a == 2 && b == 0) {
    Axpby_MV_Functor<AV, XMV, BV, YMV, 2, 0, SizeType> op (x, y, av, bv);
    Kokkos::parallel_for (policy, op);
    return;
  }

#if KOKKOSBLAS_OPTIMIZATION_LEVEL_AXPBY > 2
  if (a == 2 && b == -1) {
    Axpby_MV_Functor<AV, XMV, BV, YMV, 2, -1, SizeType> op (x, y, av, bv);
    Kokkos::parallel_for (policy, op);
    return;
  }
  if (a == 2 && b == 1) {
    Axpby_MV_Functor<AV, XMV, BV, YMV, 2, 1, SizeType> op (x, y, av, bv);
    Kokkos::parallel_for (policy, op);
    return;
  }
#endif // KOKKOSBLAS_OPTIMIZATION_LEVEL_AXPBY > 2

  // a and b arbitrary (not -1, 0, or 1)
  Axpby_MV_Functor<AV, XMV, BV, YMV, 2, 2, SizeType> op (x, y, av, bv);
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
struct
Axpby_MV_Invoke_Left {

  static void run(const AV& av, const XMV& x,
                  const BV& bv, const YMV& y,
                  int a = 2, int b = 2)
{
  static_assert (Kokkos::Impl::is_view<XMV>::value, "KokkosBlas::Impl::"
                 "Axpby_MV_Invoke_Left: X is not a Kokkos::View.");
  static_assert (Kokkos::Impl::is_view<YMV>::value, "KokkosBlas::Impl::"
                 "Axpby_MV_Invoke_Left: Y is not a Kokkos::View.");
  static_assert (Kokkos::Impl::is_same<typename YMV::value_type,
                 typename YMV::non_const_value_type>::value,
                 "KokkosBlas::Impl::Axpby_MV_Invoke_Left: Y is const.  "
                 "It must be nonconst, because it is an output argument "
                 "(we have to be able to write to its entries).");
  static_assert ((int) YMV::Rank == (int) XMV::Rank, "KokkosBlas::Impl::"
                 "Axpby_MV_Invoke_Left: X and Y must have the same rank.");
  static_assert (YMV::Rank == 2, "KokkosBlas::Impl::Axpby_MV_Invoke_Left: "
                 "X and Y must have rank 2.");

  const SizeType numCols = x.extent(1);

  // Strip-mine by 8, then 4.  After that, do one column at a time.
  // We limit the number of strip-mine values in order to keep down
  // the amount of code to compile.
  SizeType j = 0;
  for ( ; j + 8 <= numCols; j += 8) {
    XMV X_cur = Kokkos::subview (x, Kokkos::ALL (), std::make_pair (j, j+8));
    YMV Y_cur = Kokkos::subview (y, Kokkos::ALL (), std::make_pair (j, j+8));

    // Passing in the starting column index lets the functor take
    // subviews of av and bv, if they are Views.  If they are scalars,
    // the functor doesn't have to do anything to them.
    Axpby_MV_Unrolled<AV, XMV, BV, YMV, 8, SizeType> (av, X_cur, bv, Y_cur, j, a, b);
  }
  for ( ; j + 4 <= numCols; j += 4) {
    XMV X_cur = Kokkos::subview (x, Kokkos::ALL (), std::make_pair (j, j+4));
    YMV Y_cur = Kokkos::subview (y, Kokkos::ALL (), std::make_pair (j, j+4));

    // Passing in the starting column index lets the functor take
    // subviews of av and bv, if they are Views.  If they are scalars,
    // the functor doesn't have to do anything to them.
    Axpby_MV_Unrolled<AV, XMV, BV, YMV, 4, SizeType> (av, X_cur, bv, Y_cur, j, a, b);
  }
  for ( ; j < numCols; ++j) {
    auto x_cur = Kokkos::subview (x, Kokkos::ALL (), j);
    auto y_cur = Kokkos::subview (y, Kokkos::ALL (), j);

    // Passing in the starting column index lets the functor take
    // subviews of av and bv, if they are Views.  If they are scalars,
    // the functor doesn't have to do anything to them.
    typedef decltype (x_cur) XV;
    typedef decltype (y_cur) YV;
    Axpby_Generic<AV, XV, BV, YV, SizeType> (av, x_cur, bv, y_cur, j, a, b);
  }
}
};

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
struct
Axpby_MV_Invoke_Right {

static void run(const AV& av, const XMV& x,
                const BV& bv, const YMV& y,
                int a = 2, int b = 2)
{
  static_assert (Kokkos::Impl::is_view<XMV>::value, "KokkosBlas::Impl::"
                 "Axpby_MV_Invoke_Right: X is not a Kokkos::View.");
  static_assert (Kokkos::Impl::is_view<YMV>::value, "KokkosBlas::Impl::"
                 "Axpby_MV_Invoke_Right: Y is not a Kokkos::View.");
  static_assert (Kokkos::Impl::is_same<typename YMV::value_type,
                 typename YMV::non_const_value_type>::value,
                 "KokkosBlas::Impl::Axpby_MV_Invoke_Right: Y is const.  "
                 "It must be nonconst, because it is an output argument "
                 "(we have to be able to write to its entries).");
  static_assert ((int) YMV::Rank == (int) XMV::Rank, "KokkosBlas::Impl::"
                 "Axpby_MV_Invoke_Right: X and Y must have the same rank.");
  static_assert (YMV::Rank == 2, "KokkosBlas::Impl::Axpby_MV_Invoke_Right: "
                 "X and Y must have rank 2.");

  const SizeType numCols = x.extent(1);
  if (numCols == 1) {
    auto x_0 = Kokkos::subview (x, Kokkos::ALL (), 0);
    auto y_0 = Kokkos::subview (y, Kokkos::ALL (), 0);
    typedef decltype (x_0) XV;
    typedef decltype (y_0) YV;
    Axpby_Generic<AV, XV, BV, YV, SizeType> (av, x_0, bv, y_0, 0, a, b);
  }
  else {
    Axpby_MV_Generic<AV, XMV, BV, YMV, SizeType> (av, x, bv, y, a, b);
  }
}
};

}
}

#endif // KOKKOSBLAS1_AXPBY_MV_IMPL_HPP_
