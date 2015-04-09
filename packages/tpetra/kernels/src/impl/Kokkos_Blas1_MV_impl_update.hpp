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
#ifndef KOKKOS_BLAS1_MV_IMPL_UPDATE_HPP_
#define KOKKOS_BLAS1_MV_IMPL_UPDATE_HPP_

#include <TpetraKernels_config.h>
#include <Kokkos_Core.hpp>
#include <Kokkos_InnerProductSpaceTraits.hpp>

namespace KokkosBlas {
namespace Impl {

//
// update
//

// Functor for multivectors X, Y, and Z, that computes
//
// Z(i,j) = alpha*X(i,j) + beta*Y(i,j) + gamma*Z(i,j)
//
// with special cases for alpha, beta, or gamma = 0.
//
// The template parameters scalar_x, scalar_y, and scalar_z correspond
// to alpha, beta, resp. gammar in the operation Z = alpha*X + beta*Y
// + gamma*Z.  The value 0 corresponds to literal values of those
// coefficients.  The value 2 tells the functor to use the
// corresponding input coefficient.  Any literal coefficient of zero
// has BLAS semantics of ignoring the corresponding (multi)vector
// entry.
template<class XMV, class YMV, class ZMV,
         int scalar_x, int scalar_y, int scalar_z,
         class SizeType = typename ZMV::size_type>
struct MV_Update_Functor
{
  typedef typename ZMV::execution_space execution_space;
  typedef SizeType                            size_type;
  typedef Kokkos::Details::ArithTraits<typename ZMV::non_const_value_type> ATS;

  const size_type numCols;
  const typename XMV::non_const_value_type alpha_;
  XMV X_;
  const typename YMV::non_const_value_type beta_;
  YMV Y_;
  const typename ZMV::non_const_value_type gamma_;
  ZMV Z_;

  MV_Update_Functor (const typename XMV::non_const_value_type& alpha, const XMV& X,
                     const typename YMV::non_const_value_type& beta, const YMV& Y,
                     const typename ZMV::non_const_value_type& gamma, const ZMV& Z) :
    numCols (X.dimension_1 ()),
    alpha_ (alpha), X_ (X),
    beta_ (beta), Y_ (Y),
    gamma_ (gamma), Z_ (Z)
  {
#ifdef KOKKOS_HAVE_CXX11
    // XMV, YMV, and ZMV must be Kokkos::View specializations.
    static_assert (Kokkos::Impl::is_view<XMV>::value, "KokkosBlas::Impl::"
                   "MV_Update_Functor: X is not a Kokkos::View.");
    static_assert (Kokkos::Impl::is_view<YMV>::value, "KokkosBlas::Impl::"
                   "MV_Update_Functor: Y is not a Kokkos::View.");
    static_assert (Kokkos::Impl::is_view<ZMV>::value, "KokkosBlas::Impl::"
                   "MV_Update_Functor: Z is not a Kokkos::View.");
    // ZMV must be nonconst (else it can't be an output argument).
    static_assert (Kokkos::Impl::is_same<typename ZMV::value_type,
                   typename ZMV::non_const_value_type>::value,
                   "KokkosBlas::Impl::MV_Update_Functor: Z is const.  "
                   "It must be nonconst, because it is an output argument "
                   "(we have to be able to write to its entries).");
    static_assert ((int) ZMV::rank == (int) XMV::rank && (int) ZMV::rank == (int) YMV::rank,
                   "KokkosBlas::Impl::MV_Update_Functor: X, Y, and Z must have the same rank.");
    static_assert (ZMV::rank == 2, "KokkosBlas::Impl::MV_Update_Functor: "
                   "XMV, YMV, and ZMV must have rank 2.");
#endif // KOKKOS_HAVE_CXX11
  }

  KOKKOS_INLINE_FUNCTION
  void operator() (const size_type& i) const
  {
    // scalar_x, scalar_y, and scalar_z are compile-time constants
    // (since they are template parameters), so the compiler should
    // evaluate these branches at compile time.
    if (scalar_x == 0) {
      if (scalar_y == 0) {
        if (scalar_z == 0) {
#ifdef KOKKOS_HAVE_PRAGMA_IVDEP
#pragma ivdep
#endif
#ifdef KOKKOS_HAVE_PRAGMA_VECTOR
#pragma vector always
#endif
          for (size_type k = 0; k < numCols; ++k) {
            Z_(i,k) = ATS::zero ();
          }
        }
        else {
#ifdef KOKKOS_HAVE_PRAGMA_IVDEP
#pragma ivdep
#endif
#ifdef KOKKOS_HAVE_PRAGMA_VECTOR
#pragma vector always
#endif
          for (size_type k = 0; k < numCols; ++k) {
            Z_(i,k) = gamma_ * Z_(i,k);
          }
        }
      }
      else {
        if (scalar_z == 0) {
#ifdef KOKKOS_HAVE_PRAGMA_IVDEP
#pragma ivdep
#endif
#ifdef KOKKOS_HAVE_PRAGMA_VECTOR
#pragma vector always
#endif
          for (size_type k = 0; k < numCols; ++k) {
            Z_(i,k) = beta_ * Y_(i,k);
          }
        }
        else {
#ifdef KOKKOS_HAVE_PRAGMA_IVDEP
#pragma ivdep
#endif
#ifdef KOKKOS_HAVE_PRAGMA_VECTOR
#pragma vector always
#endif
          for (size_type k = 0; k < numCols; ++k) {
            Z_(i,k) = beta_ * Y_(i,k) + gamma_ * Z_(i,k);
          }
        }
      }
    }
    //
    // scalar_x == 2
    //
    else {
      if (scalar_y == 0) {
        if (scalar_z == 0) {
#ifdef KOKKOS_HAVE_PRAGMA_IVDEP
#pragma ivdep
#endif
#ifdef KOKKOS_HAVE_PRAGMA_VECTOR
#pragma vector always
#endif
          for (size_type k = 0; k < numCols; ++k) {
            Z_(i,k) = alpha_ * X_(i,k);
          }
        }
        else {
#ifdef KOKKOS_HAVE_PRAGMA_IVDEP
#pragma ivdep
#endif
#ifdef KOKKOS_HAVE_PRAGMA_VECTOR
#pragma vector always
#endif
          for (size_type k = 0; k < numCols; ++k) {
            Z_(i,k) = alpha_ * X_(i,k) + gamma_ * Z_(i,k);
          }
        }
      }
      else {
        if (scalar_z == 0) {
#ifdef KOKKOS_HAVE_PRAGMA_IVDEP
#pragma ivdep
#endif
#ifdef KOKKOS_HAVE_PRAGMA_VECTOR
#pragma vector always
#endif
          for (size_type k = 0; k < numCols; ++k) {
            Z_(i,k) = alpha_ * X_(i,k) + beta_ * Y_(i,k);
          }
        }
        else {
#ifdef KOKKOS_HAVE_PRAGMA_IVDEP
#pragma ivdep
#endif
#ifdef KOKKOS_HAVE_PRAGMA_VECTOR
#pragma vector always
#endif
          for (size_type k = 0; k < numCols; ++k) {
            Z_(i,k) = alpha_ * X_(i,k) + beta_ * Y_(i,k) + gamma_ * Z_(i,k);
          }
        }
      }
    }
  }
};

// Functor for vectors X, Y, and Z, that computes
//
// Z(i) = alpha*X(i) + beta*Y(i) + gamma*Z(i)
//
// with special cases for alpha, beta, or gamma = 0.
//
// The template parameters scalar_x, scalar_y, and scalar_z correspond
// to alpha, beta, resp. gammar in the operation Z = alpha*X + beta*Y
// + gamma*Z.  The value 0 corresponds to literal values of those
// coefficients.  The value 2 tells the functor to use the
// corresponding input coefficient.  Any literal coefficient of zero
// has BLAS semantics of ignoring the corresponding vector entry.
template<class XV, class YV, class ZV,
         int scalar_x, int scalar_y, int scalar_z,
         class SizeType = typename ZV::size_type>
struct V_Update_Functor
{
  typedef typename ZV::execution_space execution_space;
  typedef SizeType                            size_type;
  typedef Kokkos::Details::ArithTraits<typename ZV::non_const_value_type> ATS;

  const size_type numCols;
  const typename XV::non_const_value_type alpha_;
  XV X_;
  const typename YV::non_const_value_type beta_;
  YV Y_;
  const typename ZV::non_const_value_type gamma_;
  ZV Z_;

  V_Update_Functor (const typename XV::non_const_value_type& alpha, const XV& X,
                    const typename YV::non_const_value_type& beta, const YV& Y,
                    const typename ZV::non_const_value_type& gamma, const ZV& Z) :
    numCols (X.dimension_1 ()),
    alpha_ (alpha), X_ (X),
    beta_ (beta), Y_ (Y),
    gamma_ (gamma), Z_ (Z)
  {
#ifdef KOKKOS_HAVE_CXX11
    // XV, YV, and ZV must be Kokkos::View specializations.
    static_assert (Kokkos::Impl::is_view<XV>::value, "KokkosBlas::Impl::"
                   "V_Update_Functor: X is not a Kokkos::View.");
    static_assert (Kokkos::Impl::is_view<YV>::value, "KokkosBlas::Impl::"
                   "V_Update_Functor: Y is not a Kokkos::View.");
    static_assert (Kokkos::Impl::is_view<ZV>::value, "KokkosBlas::Impl::"
                   "V_Update_Functor: Z is not a Kokkos::View.");
    // ZV must be nonconst (else it can't be an output argument).
    static_assert (Kokkos::Impl::is_same<typename ZV::value_type,
                   typename ZV::non_const_value_type>::value,
                   "KokkosBlas::Impl::V_Update_Functor: Z is const.  "
                   "It must be nonconst, because it is an output argument "
                   "(we have to be able to write to its entries).");
    static_assert ((int) ZV::rank == (int) XV::rank && (int) ZV::rank == (int) YV::rank,
                   "KokkosBlas::Impl::V_Update_Functor: X, Y, and Z must have the same rank.");
    static_assert (ZV::rank == 1, "KokkosBlas::Impl::V_Update_Functor: "
                   "XV, YV, and ZV must have rank 1.");
#endif // KOKKOS_HAVE_CXX11
  }

  KOKKOS_INLINE_FUNCTION
  void operator() (const size_type& i) const
  {
    // scalar_x, scalar_y, and scalar_z are compile-time constants
    // (since they are template parameters), so the compiler should
    // evaluate these branches at compile time.
    if (scalar_x == 0) {
      if (scalar_y == 0) {
        if (scalar_z == 0) {
          Z_(i) = ATS::zero ();
        }
        else {
          Z_(i) = gamma_ * Z_(i);
        }
      }
      else {
        if (scalar_z == 0) {
          Z_(i) = beta_ * Y_(i);
        }
        else {
          Z_(i) = beta_ * Y_(i) + gamma_ * Z_(i);
        }
      }
    }
    //
    // scalar_ x == 2
    //
    else {
      if (scalar_y == 0) {
        if (scalar_z == 0) {
          Z_(i) = alpha_ * X_(i);
        }
        else {
          Z_(i) = alpha_ * X_(i) + gamma_ * Z_(i);
        }
      }
      else {
        if (scalar_z == 0) {
          Z_(i) = alpha_ * X_(i) + beta_ * Y_(i);
        }
        else {
          Z_(i) = alpha_ * X_(i) + beta_ * Y_(i) + gamma_ * Z_(i);
        }
      }
    }
  }
};

// Invoke the "generic" (not unrolled) multivector functor that
// computes
//
// Z(i,j) = alpha*X(i,j) + beta*Y(i,j) + gamma*Z(i,j)
//
// with special cases for alpha, beta, or gamma = 0.
//
// a, b, and c come in as integers.  The value 0 corresponds to the
// literal values of the coefficients.  The value 2 tells the functor
// to use the corresponding coefficients: a == 2 means use alpha, b ==
// 2 means use beta, and c == 2 means use gamma.  Otherwise, the
// corresponding coefficients are ignored.
//
// Any literal coefficient of zero has BLAS semantics of ignoring the
// corresponding multivector entry.
template<class XMV, class YMV, class ZMV, class SizeType>
void
MV_Update_Generic (const typename XMV::non_const_value_type& alpha, const XMV& X,
                   const typename YMV::non_const_value_type& beta, const YMV& Y,
                   const typename ZMV::non_const_value_type& gamma, const ZMV& Z,
                   int a = 2, int b = 2, int c = 2)
{
#ifdef KOKKOS_HAVE_CXX11
    // XMV, YMV, and ZMV must be Kokkos::View specializations.
    static_assert (Kokkos::Impl::is_view<XMV>::value, "KokkosBlas::Impl::"
                   "MV_Update_Generic: X is not a Kokkos::View.");
    static_assert (Kokkos::Impl::is_view<YMV>::value, "KokkosBlas::Impl::"
                   "MV_Update_Generic: Y is not a Kokkos::View.");
    static_assert (Kokkos::Impl::is_view<ZMV>::value, "KokkosBlas::Impl::"
                   "MV_Update_Generic: Z is not a Kokkos::View.");
    // ZMV must be nonconst (else it can't be an output argument).
    static_assert (Kokkos::Impl::is_same<typename ZMV::value_type,
                   typename ZMV::non_const_value_type>::value,
                   "KokkosBlas::Impl::MV_Update_Generic: Z is const.  "
                   "It must be nonconst, because it is an output argument "
                   "(we have to be able to write to its entries).");
    static_assert ((int) ZMV::rank == (int) XMV::rank && (int) ZMV::rank == (int) YMV::rank,
                   "KokkosBlas::Impl::MV_Update_Generic: X, Y, and Z must have the same rank.");
    static_assert (ZMV::rank == 2, "KokkosBlas::Impl::MV_Update_Generic: "
                   "XMV, YMV, and ZMV must have rank 2.");
#endif // KOKKOS_HAVE_CXX11

  typedef typename XMV::execution_space execution_space;
  const SizeType numRows = X.dimension_0 ();
  Kokkos::RangePolicy<execution_space, SizeType> policy (0, numRows);

  if (a == 0) {
    if (b == 0) {
      if (c == 0) {
        MV_Update_Functor<XMV, YMV, ZMV, 0, 0, 0, SizeType> op (alpha, X, beta, Y, gamma, Z);
        Kokkos::parallel_for (policy, op);
      }
      else {
        MV_Update_Functor<XMV, YMV, ZMV, 0, 0, 2, SizeType> op (alpha, X, beta, Y, gamma, Z);
        Kokkos::parallel_for (policy, op);
      }
    }
    else {
      if (c == 0) {
        MV_Update_Functor<XMV, YMV, ZMV, 0, 2, 0, SizeType> op (alpha, X, beta, Y, gamma, Z);
        Kokkos::parallel_for (policy, op);
      }
      else {
        MV_Update_Functor<XMV, YMV, ZMV, 0, 2, 2, SizeType> op (alpha, X, beta, Y, gamma, Z);
        Kokkos::parallel_for (policy, op);
      }
    }
  }
  //
  // a == 2
  //
  else {
    if (b == 0) {
      if (c == 0) {
        MV_Update_Functor<XMV, YMV, ZMV, 2, 0, 0, SizeType> op (alpha, X, beta, Y, gamma, Z);
        Kokkos::parallel_for (policy, op);
      }
      else {
        MV_Update_Functor<XMV, YMV, ZMV, 2, 0, 2, SizeType> op (alpha, X, beta, Y, gamma, Z);
        Kokkos::parallel_for (policy, op);
      }
    }
    else {
      if (c == 0) {
        MV_Update_Functor<XMV, YMV, ZMV, 2, 2, 0, SizeType> op (alpha, X, beta, Y, gamma, Z);
        Kokkos::parallel_for (policy, op);
      }
      else {
        MV_Update_Functor<XMV, YMV, ZMV, 2, 2, 2, SizeType> op (alpha, X, beta, Y, gamma, Z);
        Kokkos::parallel_for (policy, op);
      }
    }
  }
}


// Invoke the "generic" (not unrolled) single-vector functor that
// computes
//
// Z(i) = alpha*X(i) + beta*Y(i) + gamma*Z(i)
//
// with special cases for alpha, beta, or gamma = 0.
//
// a, b, and c come in as integers.  The value 0 corresponds to the
// literal values of the coefficients.  The value 2 tells the functor
// to use the corresponding coefficients: a == 2 means use alpha, b ==
// 2 means use beta, and c == 2 means use gamma.  Otherwise, the
// corresponding coefficients are ignored.
//
// Any literal coefficient of zero has BLAS semantics of ignoring the
// corresponding vector entry.
template<class XV, class YV, class ZV, class SizeType>
void
V_Update_Generic (const typename XV::non_const_value_type& alpha, const XV& X,
                  const typename YV::non_const_value_type& beta, const YV& Y,
                  const typename ZV::non_const_value_type& gamma, const ZV& Z,
                  int a = 2, int b = 2, int c = 2)
{
#ifdef KOKKOS_HAVE_CXX11
    // XV, YV, and ZV must be Kokkos::View specializations.
    static_assert (Kokkos::Impl::is_view<XV>::value, "KokkosBlas::Impl::"
                   "V_Update_Generic: X is not a Kokkos::View.");
    static_assert (Kokkos::Impl::is_view<YV>::value, "KokkosBlas::Impl::"
                   "V_Update_Generic: Y is not a Kokkos::View.");
    static_assert (Kokkos::Impl::is_view<ZV>::value, "KokkosBlas::Impl::"
                   "V_Update_Generic: Z is not a Kokkos::View.");
    // ZV must be nonconst (else it can't be an output argument).
    static_assert (Kokkos::Impl::is_same<typename ZV::value_type,
                   typename ZV::non_const_value_type>::value,
                   "KokkosBlas::Impl::V_Update_Generic: Z is const.  "
                   "It must be nonconst, because it is an output argument "
                   "(we have to be able to write to its entries).");
    static_assert ((int) ZV::rank == (int) XV::rank && (int) ZV::rank == (int) YV::rank,
                   "KokkosBlas::Impl::V_Update_Generic: X, Y, and Z must have the same rank.");
    static_assert (ZV::rank == 1, "KokkosBlas::Impl::V_Update_Generic: "
                   "XV, YV, and ZV must have rank 1.");
#endif // KOKKOS_HAVE_CXX11

  typedef typename XV::execution_space execution_space;
  const SizeType numRows = X.dimension_0 ();
  Kokkos::RangePolicy<execution_space, SizeType> policy (0, numRows);

  if (a == 0) {
    if (b == 0) {
      if (c == 0) {
        V_Update_Functor<XV, YV, ZV, 0, 0, 0, SizeType> op (alpha, X, beta, Y, gamma, Z);
        Kokkos::parallel_for (policy, op);
      }
      else {
        V_Update_Functor<XV, YV, ZV, 0, 0, 2, SizeType> op (alpha, X, beta, Y, gamma, Z);
        Kokkos::parallel_for (policy, op);
      }
    }
    else {
      if (c == 0) {
        V_Update_Functor<XV, YV, ZV, 0, 2, 0, SizeType> op (alpha, X, beta, Y, gamma, Z);
        Kokkos::parallel_for (policy, op);
      }
      else {
        V_Update_Functor<XV, YV, ZV, 0, 2, 2, SizeType> op (alpha, X, beta, Y, gamma, Z);
        Kokkos::parallel_for (policy, op);
      }
    }
  }
  //
  // a == 2
  //
  else {
    if (b == 0) {
      if (c == 0) {
        V_Update_Functor<XV, YV, ZV, 2, 0, 0, SizeType> op (alpha, X, beta, Y, gamma, Z);
        Kokkos::parallel_for (policy, op);
      }
      else {
        V_Update_Functor<XV, YV, ZV, 2, 0, 2, SizeType> op (alpha, X, beta, Y, gamma, Z);
        Kokkos::parallel_for (policy, op);
      }
    }
    else {
      if (c == 0) {
        V_Update_Functor<XV, YV, ZV, 2, 2, 0, SizeType> op (alpha, X, beta, Y, gamma, Z);
        Kokkos::parallel_for (policy, op);
      }
      else {
        V_Update_Functor<XV, YV, ZV, 2, 2, 2, SizeType> op (alpha, X, beta, Y, gamma, Z);
        Kokkos::parallel_for (policy, op);
      }
    }
  }
}

/// \brief Implementation of KokkosBlas::update for single vectors and
///   multivectors.
///
/// Compute
///
/// Z(i,j) = alpha*X(i,j) + beta*Y(i,j) + gamma*Z(i,j),
///
/// with special cases for alpha, beta, or gamma = 0.
template<class XMV, class YMV, class ZMV, int rank = ZMV::rank>
struct Update {};

// Partial specialization for XMV, YMV, and ZMV rank-2 Views.
template<class XMV, class YMV, class ZMV>
struct Update<XMV, YMV, ZMV, 2>
{
  typedef typename XMV::size_type size_type;
  typedef Kokkos::Details::ArithTraits<typename XMV::non_const_value_type> ATA;
  typedef Kokkos::Details::ArithTraits<typename YMV::non_const_value_type> ATB;
  typedef Kokkos::Details::ArithTraits<typename ZMV::non_const_value_type> ATC;

  static void
  update (const typename XMV::non_const_value_type& alpha, const XMV& X,
          const typename YMV::non_const_value_type& beta, const YMV& Y,
          const typename ZMV::non_const_value_type& gamma, const ZMV& Z)
  {
#ifdef KOKKOS_HAVE_CXX11
    // XMV, YMV, and ZMV must be Kokkos::View specializations.
    static_assert (Kokkos::Impl::is_view<XMV>::value, "KokkosBlas::Impl::"
                   "Update<rank 2>::update: X is not a Kokkos::View.");
    static_assert (Kokkos::Impl::is_view<YMV>::value, "KokkosBlas::Impl::"
                   "Update<rank 2>::update: Y is not a Kokkos::View.");
    static_assert (Kokkos::Impl::is_view<ZMV>::value, "KokkosBlas::Impl::"
                   "Update<rank 2>::update: Z is not a Kokkos::View.");
    // ZMV must be nonconst (else it can't be an output argument).
    static_assert (Kokkos::Impl::is_same<typename ZMV::value_type,
                     typename ZMV::non_const_value_type>::value,
                   "KokkosBlas::Impl::Update<rank 2>::update: Z is const.  "
                   "It must be nonconst, because it is an output argument "
                   "(we have to be able to write to its entries).");
    static_assert ((int) ZMV::rank == (int) XMV::rank && (int) ZMV::rank == (int) YMV::rank,
                   "KokkosBlas::Impl::Update<rank 2>::update: "
                   "X, Y, and Z must have the same rank.");
    static_assert (ZMV::rank == 2, "KokkosBlas::Impl::Update<rank 2>::update: "
                   "XMV, YMV, and ZMV must have rank 2.");
#endif // KOKKOS_HAVE_CXX11

    const size_type numRows = X.dimension_0 ();
    const size_type numCols = X.dimension_1 ();
    int a = 2, b = 2, c = 2;

    if (alpha == ATA::zero ()) {
      a = 0;
    }
    else {
      a = 2;
    }
    if (beta == ATB::zero ()) {
      b = 0;
    }
    else {
      b = 2;
    }
    if (gamma == ATC::zero ()) {
      c = 0;
    }
    else {
      c = 2;
    }

    if (numRows < static_cast<size_type> (INT_MAX) &&
        numRows * numCols < static_cast<size_type> (INT_MAX)) {
      typedef int index_type;
      MV_Update_Generic<XMV, YMV, ZMV, index_type> (alpha, X, beta, Y, gamma, Z, a, b, c);
    }
    else {
      typedef typename XMV::size_type index_type;
      MV_Update_Generic<XMV, YMV, ZMV, index_type> (alpha, X, beta, Y, gamma, Z, a, b, c);
    }
  }
};

//
// Specializations for different Kokkos devices.
//

#ifdef KOKKOS_HAVE_SERIAL
#define KOKKOSBLAS_IMPL_MV_EXEC_SPACE Kokkos::Serial
#define KOKKOSBLAS_IMPL_MV_MEM_SPACE Kokkos::HostSpace

template<>
struct Update<Kokkos::View<const double**,
                           Kokkos::LayoutLeft,
                           Kokkos::Device<KOKKOSBLAS_IMPL_MV_EXEC_SPACE, KOKKOSBLAS_IMPL_MV_MEM_SPACE>,
                           Kokkos::MemoryTraits<Kokkos::Unmanaged>,
                           Kokkos::Impl::ViewDefault>,
              Kokkos::View<const double**,
                           Kokkos::LayoutLeft,
                           Kokkos::Device<KOKKOSBLAS_IMPL_MV_EXEC_SPACE, KOKKOSBLAS_IMPL_MV_MEM_SPACE>,
                           Kokkos::MemoryTraits<Kokkos::Unmanaged>,
                           Kokkos::Impl::ViewDefault>,
              Kokkos::View<double**,
                           Kokkos::LayoutLeft,
                           Kokkos::Device<KOKKOSBLAS_IMPL_MV_EXEC_SPACE, KOKKOSBLAS_IMPL_MV_MEM_SPACE>,
                           Kokkos::MemoryTraits<Kokkos::Unmanaged>,
                           Kokkos::Impl::ViewDefault>,
              2>
{
  typedef Kokkos::View<const double**,
                       Kokkos::LayoutLeft,
                       Kokkos::Device<KOKKOSBLAS_IMPL_MV_EXEC_SPACE, KOKKOSBLAS_IMPL_MV_MEM_SPACE>,
                       Kokkos::MemoryTraits<Kokkos::Unmanaged>,
                       Kokkos::Impl::ViewDefault> XMV;
  typedef Kokkos::View<const double**,
                       Kokkos::LayoutLeft,
                       Kokkos::Device<KOKKOSBLAS_IMPL_MV_EXEC_SPACE, KOKKOSBLAS_IMPL_MV_MEM_SPACE>,
                       Kokkos::MemoryTraits<Kokkos::Unmanaged>,
                       Kokkos::Impl::ViewDefault> YMV;
  typedef Kokkos::View<double**,
                       Kokkos::LayoutLeft,
                       Kokkos::Device<KOKKOSBLAS_IMPL_MV_EXEC_SPACE, KOKKOSBLAS_IMPL_MV_MEM_SPACE>,
                       Kokkos::MemoryTraits<Kokkos::Unmanaged>,
                       Kokkos::Impl::ViewDefault> ZMV;
  typedef XMV::size_type size_type;
  typedef Kokkos::Details::ArithTraits<XMV::non_const_value_type> ATA;
  typedef Kokkos::Details::ArithTraits<YMV::non_const_value_type> ATB;
  typedef Kokkos::Details::ArithTraits<ZMV::non_const_value_type> ATC;

  static void
  update (const XMV::non_const_value_type& alpha, const XMV& X,
          const YMV::non_const_value_type& beta, const YMV& Y,
          const ZMV::non_const_value_type& gamma, const ZMV& Z);
};

#undef KOKKOSBLAS_IMPL_MV_EXEC_SPACE
#undef KOKKOSBLAS_IMPL_MV_MEM_SPACE
#endif // KOKKOS_HAVE_SERIAL

#ifdef KOKKOS_HAVE_OPENMP
#define KOKKOSBLAS_IMPL_MV_EXEC_SPACE Kokkos::OpenMP
#define KOKKOSBLAS_IMPL_MV_MEM_SPACE Kokkos::HostSpace

template<>
struct Update<Kokkos::View<const double**,
                           Kokkos::LayoutLeft,
                           Kokkos::Device<KOKKOSBLAS_IMPL_MV_EXEC_SPACE, KOKKOSBLAS_IMPL_MV_MEM_SPACE>,
                           Kokkos::MemoryTraits<Kokkos::Unmanaged>,
                           Kokkos::Impl::ViewDefault>,
              Kokkos::View<const double**,
                           Kokkos::LayoutLeft,
                           Kokkos::Device<KOKKOSBLAS_IMPL_MV_EXEC_SPACE, KOKKOSBLAS_IMPL_MV_MEM_SPACE>,
                           Kokkos::MemoryTraits<Kokkos::Unmanaged>,
                           Kokkos::Impl::ViewDefault>,
              Kokkos::View<double**,
                           Kokkos::LayoutLeft,
                           Kokkos::Device<KOKKOSBLAS_IMPL_MV_EXEC_SPACE, KOKKOSBLAS_IMPL_MV_MEM_SPACE>,
                           Kokkos::MemoryTraits<Kokkos::Unmanaged>,
                           Kokkos::Impl::ViewDefault>,
              2>
{
  typedef Kokkos::View<const double**,
                       Kokkos::LayoutLeft,
                       Kokkos::Device<KOKKOSBLAS_IMPL_MV_EXEC_SPACE, KOKKOSBLAS_IMPL_MV_MEM_SPACE>,
                       Kokkos::MemoryTraits<Kokkos::Unmanaged>,
                       Kokkos::Impl::ViewDefault> XMV;
  typedef Kokkos::View<const double**,
                       Kokkos::LayoutLeft,
                       Kokkos::Device<KOKKOSBLAS_IMPL_MV_EXEC_SPACE, KOKKOSBLAS_IMPL_MV_MEM_SPACE>,
                       Kokkos::MemoryTraits<Kokkos::Unmanaged>,
                       Kokkos::Impl::ViewDefault> YMV;
  typedef Kokkos::View<double**,
                       Kokkos::LayoutLeft,
                       Kokkos::Device<KOKKOSBLAS_IMPL_MV_EXEC_SPACE, KOKKOSBLAS_IMPL_MV_MEM_SPACE>,
                       Kokkos::MemoryTraits<Kokkos::Unmanaged>,
                       Kokkos::Impl::ViewDefault> ZMV;
  typedef XMV::size_type size_type;
  typedef Kokkos::Details::ArithTraits<XMV::non_const_value_type> ATA;
  typedef Kokkos::Details::ArithTraits<YMV::non_const_value_type> ATB;
  typedef Kokkos::Details::ArithTraits<ZMV::non_const_value_type> ATC;

  static void
  update (const XMV::non_const_value_type& alpha, const XMV& X,
          const YMV::non_const_value_type& beta, const YMV& Y,
          const ZMV::non_const_value_type& gamma, const ZMV& Z);
};

#undef KOKKOSBLAS_IMPL_MV_EXEC_SPACE
#undef KOKKOSBLAS_IMPL_MV_MEM_SPACE
#endif // KOKKOS_HAVE_OPENMP

#ifdef KOKKOS_HAVE_PTHREAD
#define KOKKOSBLAS_IMPL_MV_EXEC_SPACE Kokkos::Threads
#define KOKKOSBLAS_IMPL_MV_MEM_SPACE Kokkos::HostSpace

template<>
struct Update<Kokkos::View<const double**,
                           Kokkos::LayoutLeft,
                           Kokkos::Device<KOKKOSBLAS_IMPL_MV_EXEC_SPACE, KOKKOSBLAS_IMPL_MV_MEM_SPACE>,
                           Kokkos::MemoryTraits<Kokkos::Unmanaged>,
                           Kokkos::Impl::ViewDefault>,
              Kokkos::View<const double**,
                           Kokkos::LayoutLeft,
                           Kokkos::Device<KOKKOSBLAS_IMPL_MV_EXEC_SPACE, KOKKOSBLAS_IMPL_MV_MEM_SPACE>,
                           Kokkos::MemoryTraits<Kokkos::Unmanaged>,
                           Kokkos::Impl::ViewDefault>,
              Kokkos::View<double**,
                           Kokkos::LayoutLeft,
                           Kokkos::Device<KOKKOSBLAS_IMPL_MV_EXEC_SPACE, KOKKOSBLAS_IMPL_MV_MEM_SPACE>,
                           Kokkos::MemoryTraits<Kokkos::Unmanaged>,
                           Kokkos::Impl::ViewDefault>,
              2>
{
  typedef Kokkos::View<const double**,
                       Kokkos::LayoutLeft,
                       Kokkos::Device<KOKKOSBLAS_IMPL_MV_EXEC_SPACE, KOKKOSBLAS_IMPL_MV_MEM_SPACE>,
                       Kokkos::MemoryTraits<Kokkos::Unmanaged>,
                       Kokkos::Impl::ViewDefault> XMV;
  typedef Kokkos::View<const double**,
                       Kokkos::LayoutLeft,
                       Kokkos::Device<KOKKOSBLAS_IMPL_MV_EXEC_SPACE, KOKKOSBLAS_IMPL_MV_MEM_SPACE>,
                       Kokkos::MemoryTraits<Kokkos::Unmanaged>,
                       Kokkos::Impl::ViewDefault> YMV;
  typedef Kokkos::View<double**,
                       Kokkos::LayoutLeft,
                       Kokkos::Device<KOKKOSBLAS_IMPL_MV_EXEC_SPACE, KOKKOSBLAS_IMPL_MV_MEM_SPACE>,
                       Kokkos::MemoryTraits<Kokkos::Unmanaged>,
                       Kokkos::Impl::ViewDefault> ZMV;
  typedef XMV::size_type size_type;
  typedef Kokkos::Details::ArithTraits<XMV::non_const_value_type> ATA;
  typedef Kokkos::Details::ArithTraits<YMV::non_const_value_type> ATB;
  typedef Kokkos::Details::ArithTraits<ZMV::non_const_value_type> ATC;

  static void
  update (const XMV::non_const_value_type& alpha, const XMV& X,
          const YMV::non_const_value_type& beta, const YMV& Y,
          const ZMV::non_const_value_type& gamma, const ZMV& Z);
};

#undef KOKKOSBLAS_IMPL_MV_EXEC_SPACE
#undef KOKKOSBLAS_IMPL_MV_MEM_SPACE
#endif // KOKKOS_HAVE_PTHREAD

#ifdef KOKKOS_HAVE_CUDA
#define KOKKOSBLAS_IMPL_MV_EXEC_SPACE Kokkos::Cuda
#define KOKKOSBLAS_IMPL_MV_MEM_SPACE Kokkos::CudaSpace

template<>
struct Update<Kokkos::View<const double**,
                           Kokkos::LayoutLeft,
                           Kokkos::Device<KOKKOSBLAS_IMPL_MV_EXEC_SPACE, KOKKOSBLAS_IMPL_MV_MEM_SPACE>,
                           Kokkos::MemoryTraits<Kokkos::Unmanaged>,
                           Kokkos::Impl::ViewDefault>,
              Kokkos::View<const double**,
                           Kokkos::LayoutLeft,
                           Kokkos::Device<KOKKOSBLAS_IMPL_MV_EXEC_SPACE, KOKKOSBLAS_IMPL_MV_MEM_SPACE>,
                           Kokkos::MemoryTraits<Kokkos::Unmanaged>,
                           Kokkos::Impl::ViewDefault>,
              Kokkos::View<double**,
                           Kokkos::LayoutLeft,
                           Kokkos::Device<KOKKOSBLAS_IMPL_MV_EXEC_SPACE, KOKKOSBLAS_IMPL_MV_MEM_SPACE>,
                           Kokkos::MemoryTraits<Kokkos::Unmanaged>,
                           Kokkos::Impl::ViewDefault>,
              2>
{
  typedef Kokkos::View<const double**,
                       Kokkos::LayoutLeft,
                       Kokkos::Device<KOKKOSBLAS_IMPL_MV_EXEC_SPACE, KOKKOSBLAS_IMPL_MV_MEM_SPACE>,
                       Kokkos::MemoryTraits<Kokkos::Unmanaged>,
                       Kokkos::Impl::ViewDefault> XMV;
  typedef Kokkos::View<const double**,
                       Kokkos::LayoutLeft,
                       Kokkos::Device<KOKKOSBLAS_IMPL_MV_EXEC_SPACE, KOKKOSBLAS_IMPL_MV_MEM_SPACE>,
                       Kokkos::MemoryTraits<Kokkos::Unmanaged>,
                       Kokkos::Impl::ViewDefault> YMV;
  typedef Kokkos::View<double**,
                       Kokkos::LayoutLeft,
                       Kokkos::Device<KOKKOSBLAS_IMPL_MV_EXEC_SPACE, KOKKOSBLAS_IMPL_MV_MEM_SPACE>,
                       Kokkos::MemoryTraits<Kokkos::Unmanaged>,
                       Kokkos::Impl::ViewDefault> ZMV;
  typedef XMV::size_type size_type;
  typedef Kokkos::Details::ArithTraits<XMV::non_const_value_type> ATA;
  typedef Kokkos::Details::ArithTraits<YMV::non_const_value_type> ATB;
  typedef Kokkos::Details::ArithTraits<ZMV::non_const_value_type> ATC;

  static void
  update (const XMV::non_const_value_type& alpha, const XMV& X,
          const YMV::non_const_value_type& beta, const YMV& Y,
          const ZMV::non_const_value_type& gamma, const ZMV& Z);
};

#undef KOKKOSBLAS_IMPL_MV_EXEC_SPACE
#undef KOKKOSBLAS_IMPL_MV_MEM_SPACE
#endif // KOKKOS_HAVE_CUDA

#ifdef KOKKOS_HAVE_CUDA
#define KOKKOSBLAS_IMPL_MV_EXEC_SPACE Kokkos::Cuda
#define KOKKOSBLAS_IMPL_MV_MEM_SPACE Kokkos::CudaUVMSpace

template<>
struct Update<Kokkos::View<const double**,
                           Kokkos::LayoutLeft,
                           Kokkos::Device<KOKKOSBLAS_IMPL_MV_EXEC_SPACE, KOKKOSBLAS_IMPL_MV_MEM_SPACE>,
                           Kokkos::MemoryTraits<Kokkos::Unmanaged>,
                           Kokkos::Impl::ViewDefault>,
              Kokkos::View<const double**,
                           Kokkos::LayoutLeft,
                           Kokkos::Device<KOKKOSBLAS_IMPL_MV_EXEC_SPACE, KOKKOSBLAS_IMPL_MV_MEM_SPACE>,
                           Kokkos::MemoryTraits<Kokkos::Unmanaged>,
                           Kokkos::Impl::ViewDefault>,
              Kokkos::View<double**,
                           Kokkos::LayoutLeft,
                           Kokkos::Device<KOKKOSBLAS_IMPL_MV_EXEC_SPACE, KOKKOSBLAS_IMPL_MV_MEM_SPACE>,
                           Kokkos::MemoryTraits<Kokkos::Unmanaged>,
                           Kokkos::Impl::ViewDefault>,
              2>
{
  typedef Kokkos::View<const double**,
                       Kokkos::LayoutLeft,
                       Kokkos::Device<KOKKOSBLAS_IMPL_MV_EXEC_SPACE, KOKKOSBLAS_IMPL_MV_MEM_SPACE>,
                       Kokkos::MemoryTraits<Kokkos::Unmanaged>,
                       Kokkos::Impl::ViewDefault> XMV;
  typedef Kokkos::View<const double**,
                       Kokkos::LayoutLeft,
                       Kokkos::Device<KOKKOSBLAS_IMPL_MV_EXEC_SPACE, KOKKOSBLAS_IMPL_MV_MEM_SPACE>,
                       Kokkos::MemoryTraits<Kokkos::Unmanaged>,
                       Kokkos::Impl::ViewDefault> YMV;
  typedef Kokkos::View<double**,
                       Kokkos::LayoutLeft,
                       Kokkos::Device<KOKKOSBLAS_IMPL_MV_EXEC_SPACE, KOKKOSBLAS_IMPL_MV_MEM_SPACE>,
                       Kokkos::MemoryTraits<Kokkos::Unmanaged>,
                       Kokkos::Impl::ViewDefault> ZMV;
  typedef XMV::size_type size_type;
  typedef Kokkos::Details::ArithTraits<XMV::non_const_value_type> ATA;
  typedef Kokkos::Details::ArithTraits<YMV::non_const_value_type> ATB;
  typedef Kokkos::Details::ArithTraits<ZMV::non_const_value_type> ATC;

  static void
  update (const XMV::non_const_value_type& alpha, const XMV& X,
          const YMV::non_const_value_type& beta, const YMV& Y,
          const ZMV::non_const_value_type& gamma, const ZMV& Z);
};

#undef KOKKOSBLAS_IMPL_MV_EXEC_SPACE
#undef KOKKOSBLAS_IMPL_MV_MEM_SPACE
#endif // KOKKOS_HAVE_CUDA



// Partial specialization for XV, YV, and ZV rank-1 Views.
template<class XV, class YV, class ZV>
struct Update<XV, YV, ZV, 1>
{
  typedef typename XV::size_type size_type;
  typedef Kokkos::Details::ArithTraits<typename XV::non_const_value_type> ATA;
  typedef Kokkos::Details::ArithTraits<typename YV::non_const_value_type> ATB;
  typedef Kokkos::Details::ArithTraits<typename ZV::non_const_value_type> ATC;

  static void
  update (const typename XV::non_const_value_type& alpha, const XV& X,
          const typename YV::non_const_value_type& beta, const YV& Y,
          const typename ZV::non_const_value_type& gamma, const ZV& Z)
  {
#ifdef KOKKOS_HAVE_CXX11
    // XV, YV, and ZV must be Kokkos::View specializations.
    static_assert (Kokkos::Impl::is_view<XV>::value, "KokkosBlas::Impl::"
                   "Update<rank 1>::update: X is not a Kokkos::View.");
    static_assert (Kokkos::Impl::is_view<YV>::value, "KokkosBlas::Impl::"
                   "Update<rank 1>::update: Y is not a Kokkos::View.");
    static_assert (Kokkos::Impl::is_view<ZV>::value, "KokkosBlas::Impl::"
                   "Update<rank 1>::update: Z is not a Kokkos::View.");
    // ZV must be nonconst (else it can't be an output argument).
    static_assert (Kokkos::Impl::is_same<typename ZV::value_type,
                     typename ZV::non_const_value_type>::value,
                   "KokkosBlas::Impl::Update<rank 1>::update: Z is const.  "
                   "It must be nonconst, because it is an output argument "
                   "(we have to be able to write to its entries).");
    static_assert ((int) ZV::rank == (int) XV::rank && (int) ZV::rank == (int) YV::rank,
                   "KokkosBlas::Impl::Update<rank 1>::update: "
                   "X, Y, and Z must have the same rank.");
    static_assert (ZV::rank == 1, "KokkosBlas::Impl::Update<rank 1>::update: "
                   "XV, YV, and ZV must have rank 1.");
#endif // KOKKOS_HAVE_CXX11

    const size_type numRows = X.dimension_0 ();
    const size_type numCols = X.dimension_1 ();
    int a = 2, b = 2, c = 2;

    if (alpha == ATA::zero ()) {
      a = 0;
    }
    else {
      a = 2;
    }
    if (beta == ATB::zero ()) {
      b = 0;
    }
    else {
      b = 2;
    }
    if (gamma == ATC::zero ()) {
      c = 0;
    }
    else {
      c = 2;
    }

    if (numRows < static_cast<size_type> (INT_MAX) &&
        numRows * numCols < static_cast<size_type> (INT_MAX)) {
      typedef int index_type;
      V_Update_Generic<XV, YV, ZV, index_type> (alpha, X, beta, Y, gamma, Z, a, b, c);
    }
    else {
      typedef typename XV::size_type index_type;
      V_Update_Generic<XV, YV, ZV, index_type> (alpha, X, beta, Y, gamma, Z, a, b, c);
    }
  }
};

} // namespace Impl
} // namespace KokkosBlas

#endif // KOKKOS_BLAS1_MV_IMPL_UPDATE_HPP_
