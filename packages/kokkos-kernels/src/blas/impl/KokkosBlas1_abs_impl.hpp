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
#ifndef KOKKOS_BLAS1_IMPL_ABS_HPP_
#define KOKKOS_BLAS1_IMPL_ABS_HPP_

#include <KokkosKernels_config.h>
#include <Kokkos_Core.hpp>
#include <Kokkos_ArithTraits.hpp>

namespace KokkosBlas {
namespace Impl {

//
// abs
//

// Entry-wise absolute value / magnitude: R(i,j) = abs(X(i,j)).
template<class RMV, class XMV, class SizeType = typename RMV::size_type>
struct MV_Abs_Functor
{
  typedef typename RMV::execution_space execution_space;
  typedef SizeType                            size_type;
  typedef Kokkos::Details::ArithTraits<typename XMV::non_const_value_type> ATS;

  const size_type numCols;
  RMV R_;
  XMV X_;

  MV_Abs_Functor (const RMV& R, const XMV& X) :
    numCols (X.extent(1)), R_ (R), X_ (X)
  {
    static_assert (Kokkos::Impl::is_view<RMV>::value, "KokkosBlas::Impl::"
                   "MV_Abs_Functor: RMV is not a Kokkos::View.");
    static_assert (Kokkos::Impl::is_view<XMV>::value, "KokkosBlas::Impl::"
                   "MV_Abs_Functor: XMV is not a Kokkos::View.");
    static_assert (RMV::rank == 2, "KokkosBlas::Impl::"
                   "MV_Abs_Functor: RMV is not rank 2");
    static_assert (XMV::rank == 2, "KokkosBlas::Impl::"
                   "MV_Abs_Functor: XMV is not rank 2");
  }

  KOKKOS_INLINE_FUNCTION
  void operator() (const size_type& i) const
  {
#ifdef KOKKOS_ENABLE_PRAGMA_IVDEP
#pragma ivdep
#endif
    for (size_type j = 0; j < numCols; ++j) {
      R_(i,j) = ATS::abs (X_(i,j));
    }
  }
};

// Entry-wise, in-place absolute value / magnitude: R(i,j) = abs(R(i,j)).
template<class RMV, class SizeType = typename RMV::size_type>
struct MV_AbsSelf_Functor
{
  typedef typename RMV::execution_space execution_space;
  typedef SizeType                            size_type;
  typedef Kokkos::Details::ArithTraits<typename RMV::non_const_value_type> ATS;

  const size_type numCols;
  RMV R_;

  MV_AbsSelf_Functor (const RMV& R) :
    numCols (R.extent(1)), R_ (R)
  {
    static_assert (Kokkos::Impl::is_view<RMV>::value, "KokkosBlas::Impl::"
                   "MV_Abs_Functor: RMV is not a Kokkos::View.");
    static_assert (RMV::rank == 2, "KokkosBlas::Impl::"
                   "MV_Abs_Functor: RMV is not rank 2");
  }

  KOKKOS_INLINE_FUNCTION
  void operator() (const size_type& i) const
  {
#ifdef KOKKOS_ENABLE_PRAGMA_IVDEP
#pragma ivdep
#endif
    for (size_type j = 0; j < numCols; ++j) {
      R_(i,j) = ATS::abs (R_(i,j));
    }
  }
};

// Single-vector, entry-wise absolute value / magnitude: R(i) = abs(X(i)).
template<class RV, class XV, class SizeType = typename RV::size_type>
struct V_Abs_Functor
{
  typedef typename RV::execution_space execution_space;
  typedef SizeType                            size_type;
  typedef Kokkos::Details::ArithTraits<typename XV::non_const_value_type> ATS;

  RV R_;
  XV X_;

  V_Abs_Functor (const RV& R, const XV& X) : R_ (R), X_ (X)
  {
    static_assert (Kokkos::Impl::is_view<RV>::value, "KokkosBlas::Impl::"
                   "V_Abs_Functor: RV is not a Kokkos::View.");
    static_assert (Kokkos::Impl::is_view<XV>::value, "KokkosBlas::Impl::"
                   "V_Abs_Functor: XV is not a Kokkos::View.");
    static_assert (RV::rank == 1, "KokkosBlas::Impl::"
                   "V_Abs_Functor: RV is not rank 1");
    static_assert (XV::rank == 1, "KokkosBlas::Impl::"
                   "V_Abs_Functor: XV is not rank 1");
  }

  KOKKOS_INLINE_FUNCTION
  void operator() (const size_type& i) const
  {
    R_(i) = ATS::abs (X_(i));
  }
};

// Single-vector, entry-wise, in-place absolute value / magnitude: R(i) = abs(R(i)).
template<class RV, class SizeType = typename RV::size_type>
struct V_AbsSelf_Functor
{
  typedef typename RV::execution_space execution_space;
  typedef SizeType                            size_type;
  typedef Kokkos::Details::ArithTraits<typename RV::non_const_value_type> ATS;

  RV R_;

  V_AbsSelf_Functor (const RV& R) : R_ (R)
  {
    static_assert (Kokkos::Impl::is_view<RV>::value, "KokkosBlas::Impl::"
                   "V_Abs_Functor: RV is not a Kokkos::View.");
    static_assert (RV::rank == 1, "KokkosBlas::Impl::"
                   "V_Abs_Functor: RV is not rank 1");
  }

  KOKKOS_INLINE_FUNCTION
  void operator() (const size_type& i) const
  {
    R_(i) = ATS::abs (R_(i));
  }
};

// Invoke the "generic" (not unrolled) multivector functor that
// computes entry-wise absolute value.
template<class RMV, class XMV, class SizeType>
void
MV_Abs_Generic (const RMV& R, const XMV& X)
{
  static_assert (Kokkos::Impl::is_view<RMV>::value, "KokkosBlas::Impl::"
                 "MV_Abs_Generic: RMV is not a Kokkos::View.");
  static_assert (Kokkos::Impl::is_view<XMV>::value, "KokkosBlas::Impl::"
                 "MV_Abs_Generic: XMV is not a Kokkos::View.");
  static_assert (RMV::rank == 2, "KokkosBlas::Impl::"
                 "MV_Abs_Generic: RMV is not rank 2");
  static_assert (XMV::rank == 2, "KokkosBlas::Impl::"
                 "MV_Abs_Generic: XMV is not rank 2");

  typedef typename XMV::execution_space execution_space;
  const SizeType numRows = X.extent(0);
  Kokkos::RangePolicy<execution_space, SizeType> policy (0, numRows);

  if (R == X) { // if R and X are the same (alias one another)
    MV_AbsSelf_Functor<RMV, SizeType> op (R);
    Kokkos::parallel_for (policy, op);
  }
  else {
    MV_Abs_Functor<RMV, XMV, SizeType> op (R, X);
    Kokkos::parallel_for (policy, op);
  }
}

// Variant of MV_Abs_Generic for single vectors (1-D Views) R and X.
template<class RV, class XV, class SizeType>
void
V_Abs_Generic (const RV& R, const XV& X)
{
  static_assert (Kokkos::Impl::is_view<RV>::value, "KokkosBlas::Impl::"
                 "V_Abs_Generic: RV is not a Kokkos::View.");
  static_assert (Kokkos::Impl::is_view<XV>::value, "KokkosBlas::Impl::"
                 "V_Abs_Generic: XV is not a Kokkos::View.");
  static_assert (RV::rank == 1, "KokkosBlas::Impl::"
                 "V_Abs_Generic: RV is not rank 1");
  static_assert (XV::rank == 1, "KokkosBlas::Impl::"
                 "V_Abs_Generic: XV is not rank 1");

  typedef typename XV::execution_space execution_space;
  const SizeType numRows = X.extent(0);
  Kokkos::RangePolicy<execution_space, SizeType> policy (0, numRows);

  if (R == X) { // if R and X are the same (alias one another)
    V_AbsSelf_Functor<RV, SizeType> op (R);
    Kokkos::parallel_for (policy, op);
  }
  else {
    V_Abs_Functor<RV, XV, SizeType> op (R, X);
    Kokkos::parallel_for (policy, op);
  }
}

}
}
#endif // KOKKOS_BLAS1_MV_IMPL_ABS_HPP_
