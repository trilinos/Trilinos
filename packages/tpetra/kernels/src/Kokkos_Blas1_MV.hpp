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
#ifndef KOKKOS_BLAS1_MV_HPP_
#define KOKKOS_BLAS1_MV_HPP_

#include <Kokkos_Blas1_MV_impl.hpp>
#ifdef KOKKOS_HAVE_CXX11
#  include <type_traits>
#endif // KOKKOS_HAVE_CXX11

namespace KokkosBlas {

/// \brief Compute the column-wise dot products of two multivectors.
///
/// \tparam RV 1-D output View
/// \tparam XMV 2-D input View
/// \tparam YMV 2-D input View
///
/// \param dots [out] Output 1-D View to which to write results.
/// \param x [in] Input 2-D View.
/// \param y [in] Input 2-D View.
template<class RV, class XMV, class YMV>
void
dot (const RV& dots, const XMV& x, const YMV& y)
{
#ifdef KOKKOS_HAVE_CXX11
  // RV, XMV, and YMV must be Kokkos::View specializations.
  static_assert (Kokkos::Impl::is_view<RV>::value, "KokkosBlas::dot (MultiVector): "
                 "The output argument is not a Kokkos::View.");
  static_assert (Kokkos::Impl::is_view<XMV>::value, "KokkosBlas::dot (MultiVector): "
                 "The first input argument x is not a Kokkos::View.");
  static_assert (Kokkos::Impl::is_view<YMV>::value, "KokkosBlas::dot (MultiVector): "
                 "The second input argument y is not a Kokkos::View.");
  // RV must be nonconst (else it can't be an output argument).
  static_assert (Kokkos::Impl::is_same<typename RV::value_type, typename RV::non_const_value_type>::value,
                 "KokkosBlas::dot (MultiVector): The output argument is const.  "
                 "It must be nonconst, because it is an output argument "
                 "(we have to be able to write to its entries).");
  static_assert (RV::rank == 1, "KokkosBlas::dot (MultiVector): "
                 "The output argument must have rank 1.");
  static_assert (XMV::rank == 2, "KokkosBlas::dot (MultiVector): "
                 "The first input argument x must have rank 2.");
  static_assert (YMV::rank == 2, "KokkosBlas::dot (MultiVector): "
                 "The second input argument y must have rank 2.");
#else
  // We prefer to use C++11 static_assert, because it doesn't give
  // "unused typedef" warnings, like the constructs below do.
  //
  // RV, XMV, and YMV must be Kokkos::View specializations.
  typedef typename
    Kokkos::Impl::StaticAssert<Kokkos::Impl::is_view<RV>::value>::type RVIsNotView;
  typedef typename
    Kokkos::Impl::StaticAssert<Kokkos::Impl::is_view<XMV>::value>::type XMVIsNotView;
  typedef typename
    Kokkos::Impl::StaticAssert<Kokkos::Impl::is_view<YMV>::value>::type YMVIsNotView;

  // RV must be nonconst (else it can't be an output argument).
  typedef typename
    Kokkos::Impl::StaticAssert<Kokkos::Impl::is_same<typename RV::value_type,
      typename RV::non_const_value_type>::value>::type RV_is_const;

  typedef typename
    Kokkos::Impl::StaticAssert<RV::rank == 1 >::type Blas1_Dot_MultiVectorRanksDontMatch;
  typedef typename
    Kokkos::Impl::StaticAssert<XMV::rank == 2 >::type Blas1_Dot_XMultiVectorRankNot2;
  typedef typename
    Kokkos::Impl::StaticAssert<YMV::rank == 2 >::type Blas1_Dot_YMultiVectorRankNot2;
#endif // KOKKOS_HAVE_CXX11

  // Check compatibility of dimensions at run time.
  if (x.dimension_0 () != y.dimension_0 () ||
      x.dimension_1 () != y.dimension_1 () ||
      dots.dimension_0 () != x.dimension_1 ()) {
    std::ostringstream os;
    os << "KokkosBlas::dot (MultiVector): Dimensions do not match: "
       << "dots: " << dots.dimension_0 () << " x 1"
       << ", x: " << x.dimension_0 () << " x " << x.dimension_1 ()
       << ", y: " << y.dimension_0 () << " x " << y.dimension_1 ();
    Kokkos::Impl::throw_runtime_exception (os.str ());
  }

  // Any View can be assigned to an unmanaged View, and it's safe to
  // use them here.
  typedef Kokkos::View<typename RV::non_const_value_type*,
    typename RV::array_layout,
    typename RV::device_type,
    Kokkos::MemoryTraits<Kokkos::Unmanaged>,
    typename RV::specialize> RV_Internal;
  typedef Kokkos::View<typename XMV::const_value_type**,
    typename XMV::array_layout,
    typename XMV::device_type,
    Kokkos::MemoryTraits<Kokkos::Unmanaged>,
    typename XMV::specialize> XMV_Internal;
  typedef Kokkos::View<typename YMV::const_value_type**,
    typename YMV::array_layout,
    typename YMV::device_type,
    Kokkos::MemoryTraits<Kokkos::Unmanaged>,
    typename YMV::specialize> YMV_Internal;

  RV_Internal dots_i = dots;
  XMV_Internal x_i = x;
  YMV_Internal y_i = y;

  return Impl::Dot_MV<
    typename RV_Internal::value_type*,
    typename RV_Internal::array_layout,
    typename RV_Internal::device_type,
    typename RV_Internal::memory_traits,
    typename RV_Internal::specialize,
    typename XMV_Internal::value_type**,
    typename XMV_Internal::array_layout,
    typename XMV_Internal::device_type,
    typename XMV_Internal::memory_traits,
    typename XMV_Internal::specialize,
    typename YMV_Internal::value_type**,
    typename YMV_Internal::array_layout,
    typename YMV_Internal::device_type,
    typename YMV_Internal::memory_traits,
    typename YMV_Internal::specialize
    >::dot (dots_i, x_i, y_i);
}

} // namespace KokkosBlas

#endif // KOKKOS_BLAS1_MV_HPP_
