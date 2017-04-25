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
#ifndef KOKKOS_BLAS1_MV_IMPL_RECIPROCAL_HPP_
#define KOKKOS_BLAS1_MV_IMPL_RECIPROCAL_HPP_

#include <KokkosKernels_config.h>
#include <Kokkos_Core.hpp>
#include <Kokkos_ArithTraits.hpp>

namespace KokkosBlas {
namespace Impl {

// Functor that implements entry-wise reciprocal:
//
// R(i,j) = 1 / X(i,j).
template<class RMV, class XMV, class SizeType = typename RMV::size_type>
struct MV_Reciprocal_Functor
{
  typedef typename RMV::execution_space execution_space;
  typedef SizeType                            size_type;
  typedef Kokkos::Details::ArithTraits<typename XMV::non_const_value_type> ATS;

  const size_type numCols;
  RMV R_;
  XMV X_;

  MV_Reciprocal_Functor (const RMV& R, const XMV& X) :
    numCols (X.dimension_1 ()), R_ (R), X_ (X)
  {
    static_assert (Kokkos::Impl::is_view<RMV>::value, "KokkosBlas::Impl::"
                   "MV_Reciprocal_Functor: RMV is not a Kokkos::View.");
    static_assert (Kokkos::Impl::is_view<XMV>::value, "KokkosBlas::Impl::"
                   "MV_Reciprocal_Functor: XMV is not a Kokkos::View.");
    static_assert (RMV::rank == 2, "KokkosBlas::Impl::"
                   "MV_Reciprocal_Functor: RMV is not rank 2");
    static_assert (XMV::rank == 2, "KokkosBlas::Impl::"
                   "MV_Reciprocal_Functor: XMV is not rank 2");
  }

  KOKKOS_INLINE_FUNCTION
  void operator() (const size_type& i) const
  {
#ifdef KOKKOS_HAVE_PRAGMA_IVDEP
#pragma ivdep
#endif
    for (size_type j = 0; j < numCols; ++j) {
      R_(i,j) = ATS::one () / X_(i,j);
    }
  }
};

// Entry-wise, in-place reciprocal: R(i,j) = 1 / R(i,j).
template<class RMV, class SizeType = typename RMV::size_type>
struct MV_ReciprocalSelf_Functor
{
  typedef typename RMV::execution_space execution_space;
  typedef SizeType                            size_type;
  typedef Kokkos::Details::ArithTraits<typename RMV::non_const_value_type> ATS;

  const size_type numCols;
  RMV R_;

  MV_ReciprocalSelf_Functor (const RMV& R) :
    numCols (R.dimension_1 ()), R_ (R)
  {
    static_assert (Kokkos::Impl::is_view<RMV>::value, "KokkosBlas::Impl::"
                   "MV_Reciprocal_Functor: RMV is not a Kokkos::View.");
    static_assert (RMV::rank == 2, "KokkosBlas::Impl::"
                   "MV_Reciprocal_Functor: RMV is not rank 2");
  }

  KOKKOS_INLINE_FUNCTION
  void operator() (const size_type& i) const
  {
#ifdef KOKKOS_HAVE_PRAGMA_IVDEP
#pragma ivdep
#endif
    for (size_type j = 0; j < numCols; ++j) {
      R_(i,j) = ATS::one () / R_(i,j);
    }
  }
};

// Single-vector, entry-wise reciprocal: R(i) = 1 / R(i).
template<class RV, class XV, class SizeType = typename RV::size_type>
struct V_Reciprocal_Functor
{
  typedef typename RV::execution_space execution_space;
  typedef SizeType                            size_type;
  typedef Kokkos::Details::ArithTraits<typename XV::non_const_value_type> ATS;

  RV R_;
  XV X_;

  V_Reciprocal_Functor (const RV& R, const XV& X) : R_ (R), X_ (X)
  {
    static_assert (Kokkos::Impl::is_view<RV>::value, "KokkosBlas::Impl::"
                   "V_Reciprocal_Functor: RV is not a Kokkos::View.");
    static_assert (Kokkos::Impl::is_view<XV>::value, "KokkosBlas::Impl::"
                   "V_Reciprocal_Functor: XV is not a Kokkos::View.");
    static_assert (RV::rank == 1, "KokkosBlas::Impl::"
                   "V_Reciprocal_Functor: RV is not rank 1");
    static_assert (XV::rank == 1, "KokkosBlas::Impl::"
                   "V_Reciprocal_Functor: XV is not rank 1");
  }

  KOKKOS_INLINE_FUNCTION
  void operator() (const size_type& i) const
  {
    R_(i) = ATS::one () / X_(i);
  }
};

// Single-vector, entry-wise, in-place reciprocal: R(i) = 1 / R(i).
template<class RV, class SizeType = typename RV::size_type>
struct V_ReciprocalSelf_Functor
{
  typedef typename RV::execution_space execution_space;
  typedef SizeType                            size_type;
  typedef Kokkos::Details::ArithTraits<typename RV::non_const_value_type> ATS;

  RV R_;

  V_ReciprocalSelf_Functor (const RV& R) : R_ (R)
  {
    static_assert (Kokkos::Impl::is_view<RV>::value, "KokkosBlas::Impl::"
                   "V_Reciprocal_Functor: RV is not a Kokkos::View.");
    static_assert (RV::rank == 1, "KokkosBlas::Impl::"
                   "V_Reciprocal_Functor: RV is not rank 1");
  }

  KOKKOS_INLINE_FUNCTION
  void operator() (const size_type& i) const
  {
    R_(i) = ATS::one () / R_(i);
  }
};

// Invoke the "generic" (not unrolled) multivector functor that
// computes entry-wise reciprocal.
template<class RMV, class XMV, class SizeType>
void
MV_Reciprocal_Generic (const RMV& R, const XMV& X)
{
  static_assert (Kokkos::Impl::is_view<RMV>::value, "KokkosBlas::Impl::"
                 "MV_Reciprocal_Generic: RMV is not a Kokkos::View.");
  static_assert (Kokkos::Impl::is_view<XMV>::value, "KokkosBlas::Impl::"
                 "MV_Reciprocal_Generic: XMV is not a Kokkos::View.");
  static_assert (RMV::rank == 2, "KokkosBlas::Impl::"
                 "MV_Reciprocal_Generic: RMV is not rank 2");
  static_assert (XMV::rank == 2, "KokkosBlas::Impl::"
                 "MV_Reciprocal_Generic: XMV is not rank 2");

  typedef typename XMV::execution_space execution_space;
  const SizeType numRows = X.dimension_0 ();
  Kokkos::RangePolicy<execution_space, SizeType> policy (0, numRows);

  if (R == X) { // if R and X are the same (alias one another)
    MV_ReciprocalSelf_Functor<RMV, SizeType> op (R);
    Kokkos::parallel_for (policy, op);
  }
  else {
    MV_Reciprocal_Functor<RMV, XMV, SizeType> op (R, X);
    Kokkos::parallel_for (policy, op);
  }
}

// Variant of MV_Reciprocal_Generic for single vectors (1-D Views) R and X.
template<class RV, class XV, class SizeType>
void
V_Reciprocal_Generic (const RV& R, const XV& X)
{
  static_assert (Kokkos::Impl::is_view<RV>::value, "KokkosBlas::Impl::"
                 "V_Reciprocal_Generic: RV is not a Kokkos::View.");
  static_assert (Kokkos::Impl::is_view<XV>::value, "KokkosBlas::Impl::"
                 "V_Reciprocal_Generic: XV is not a Kokkos::View.");
  static_assert (RV::rank == 1, "KokkosBlas::Impl::"
                 "V_Reciprocal_Generic: RV is not rank 1");
  static_assert (XV::rank == 1, "KokkosBlas::Impl::"
                 "V_Reciprocal_Generic: XV is not rank 1");

  typedef typename XV::execution_space execution_space;
  const SizeType numRows = X.dimension_0 ();
  Kokkos::RangePolicy<execution_space, SizeType> policy (0, numRows);

  if (R == X) { // if R and X are the same (alias one another)
    V_ReciprocalSelf_Functor<RV, SizeType> op (R);
    Kokkos::parallel_for (policy, op);
  }
  else {
    V_Reciprocal_Functor<RV, XV, SizeType> op (R, X);
    Kokkos::parallel_for (policy, op);
  }
}

//! Implementation of KokkosBlas::reciprocal for (multi)vectors.
template<class RMV, class XMV, int rank = RMV::rank>
struct Reciprocal;

template<class RMV, class XMV>
struct Reciprocal<RMV, XMV, 2>
#ifndef KOKKOSKERNELS_ETI_ONLY
{
  typedef typename XMV::size_type size_type;

  static void reciprocal (const RMV& R, const XMV& X)
  {
    static_assert (Kokkos::Impl::is_view<RMV>::value, "KokkosBlas::Impl::"
                   "Reciprocal<2-D>: RMV is not a Kokkos::View.");
    static_assert (Kokkos::Impl::is_view<XMV>::value, "KokkosBlas::Impl::"
                   "Reciprocal<2-D>: XMV is not a Kokkos::View.");
    static_assert (RMV::rank == 2, "KokkosBlas::Impl::Reciprocal<2-D>: "
                   "RMV is not rank 2.");
    static_assert (XMV::rank == 2, "KokkosBlas::Impl::Reciprocal<2-D>: "
                   "XMV is not rank 2.");

    const size_type numRows = X.dimension_0 ();
    const size_type numCols = X.dimension_1 ();
    if (numRows < static_cast<size_type> (INT_MAX) &&
        numRows * numCols < static_cast<size_type> (INT_MAX)) {
      typedef int index_type;
      MV_Reciprocal_Generic<RMV, XMV, index_type> (R, X);
    }
    else {
      typedef typename XMV::size_type index_type;
      MV_Reciprocal_Generic<RMV, XMV, index_type> (R, X);
    }
  }
}
#endif
;

//! Partial specialization of Reciprocal for single vectors (1-D Views).
template<class RMV, class XMV>
struct Reciprocal<RMV, XMV, 1>
#ifndef KOKKOSKERNELS_ETI_ONLY
{
  typedef typename XMV::size_type size_type;

  static void reciprocal (const RMV& R, const XMV& X)
  {
    static_assert (Kokkos::Impl::is_view<RMV>::value, "KokkosBlas::Impl::"
                   "Reciprocal<1-D>: RMV is not a Kokkos::View.");
    static_assert (Kokkos::Impl::is_view<XMV>::value, "KokkosBlas::Impl::"
                   "Reciprocal<1-D>: XMV is not a Kokkos::View.");
    static_assert (RMV::rank == 1, "KokkosBlas::Impl::Reciprocal<1-D>: "
                   "RMV is not rank 1.");
    static_assert (XMV::rank == 1, "KokkosBlas::Impl::Reciprocal<1-D>: "
                   "XMV is not rank 1.");

    const size_type numRows = X.dimension_0 ();

    if (numRows < static_cast<size_type> (INT_MAX)) {
      typedef int index_type;
      V_Reciprocal_Generic<RMV, XMV, index_type> (R, X);
    }
    else {
      typedef typename XMV::size_type index_type;
      V_Reciprocal_Generic<RMV, XMV, index_type> (R, X);
    }
  }
}
#endif
;

//
// Macro for declaration of full specialization of
// KokkosBlas::Impl::Reciprocal for rank == 2.  This is NOT for
// users!!!  All the declarations of full specializations go in this
// header file.  We may spread out definitions (see _DEF macro below)
// across one or more .cpp files.
//

#define KOKKOSBLAS1_IMPL_MV_RECIP_DECL( SCALAR, LAYOUT, EXEC_SPACE, MEM_SPACE ) \
extern template struct Reciprocal<Kokkos::View<SCALAR**, \
                               LAYOUT, \
                               Kokkos::Device<EXEC_SPACE, MEM_SPACE>, \
                               Kokkos::MemoryTraits<Kokkos::Unmanaged> >, \
                  Kokkos::View<const SCALAR**, \
                               LAYOUT, \
                               Kokkos::Device<EXEC_SPACE, MEM_SPACE>, \
                               Kokkos::MemoryTraits<Kokkos::Unmanaged> >, \
                                  2>;

//
// Macro for declaration of full specialization of
// KokkosBlas::Impl::Reciprocal for rank == 2.  This is NOT for users!!!
//

#define KOKKOSBLAS1_IMPL_MV_RECIP_DEF( SCALAR, LAYOUT, EXEC_SPACE, MEM_SPACE ) \
template struct Reciprocal<Kokkos::View<SCALAR**, \
                        LAYOUT, \
                        Kokkos::Device<EXEC_SPACE, MEM_SPACE>, \
                        Kokkos::MemoryTraits<Kokkos::Unmanaged> >, \
           Kokkos::View<const SCALAR**, \
                        LAYOUT, \
                        Kokkos::Device<EXEC_SPACE, MEM_SPACE>, \
                        Kokkos::MemoryTraits<Kokkos::Unmanaged> >, \
                           2>;

} // namespace Impl
} // namespace KokkosBlas
#include<generated_specializations_hpp/KokkosBlas1_impl_MV_recip_decl_specializations.hpp>
#endif // KOKKOS_BLAS1_MV_IMPL_RECIPROCAL_HPP_
