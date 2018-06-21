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

#ifndef KOKKOSBLAS1_AXPBY_HPP_
#define KOKKOSBLAS1_AXPBY_HPP_

#include<KokkosBlas1_axpby_spec.hpp>
#include<KokkosKernels_helpers.hpp>

// axpby() accepts both scalar coefficients a and b, and vector
// coefficients (apply one for each column of the input multivectors).
// This traits class helps axpby() select the correct specialization
// of AV and BV (the type of a resp. b) for invoking the
// implementation.

namespace KokkosBlas {


template<class AV, class XMV, class BV, class YMV>
void
axpby (const AV& a, const XMV& X, const BV& b, const YMV& Y)
{
  static_assert (Kokkos::Impl::is_view<XMV>::value, "KokkosBlas::axpby: "
                 "X is not a Kokkos::View.");
  static_assert (Kokkos::Impl::is_view<YMV>::value, "KokkosBlas::axpby: "
                 "Y is not a Kokkos::View.");
  static_assert (Kokkos::Impl::is_same<typename YMV::value_type,
                 typename YMV::non_const_value_type>::value,
                 "KokkosBlas::axpby: Y is const.  It must be nonconst, "
                 "because it is an output argument "
                 "(we must be able to write to its entries).");
  static_assert (int(YMV::Rank) == int(XMV::Rank), "KokkosBlas::axpby: "
                 "X and Y must have the same rank.");
  static_assert (YMV::Rank == 1 || YMV::Rank == 2, "KokkosBlas::axpby: "
                 "XMV and YMV must either have rank 1 or rank 2.");

  // Check compatibility of dimensions at run time.
  if (X.extent(0) != Y.extent(0) ||
      X.extent(1) != Y.extent(1)) {
    std::ostringstream os;
    os << "KokkosBlas::axpby: Dimensions of X and Y do not match: "
       << "X: " << X.extent(0) << " x " << X.extent(1)
       << ", Y: " << Y.extent(0) << " x " << Y.extent(1);
    Kokkos::Impl::throw_runtime_exception (os.str ());
  }

  // Create unmanaged versions of the input Views.  XMV and YMV may be
  // rank 1 or rank 2.  AV and BV may be either rank-1 Views, or
  // scalar values.
  typedef typename KokkosKernels::Impl::GetUnifiedScalarViewType<
    AV, XMV, true>::type AV_Internal;
  typedef Kokkos::View<
    typename Kokkos::Impl::if_c<
      XMV::Rank == 1,
      typename XMV::const_value_type*,
      typename XMV::const_value_type** >::type,
    typename KokkosKernels::Impl::GetUnifiedLayout<XMV>::array_layout,
    typename XMV::device_type,
    Kokkos::MemoryTraits<Kokkos::Unmanaged> > XMV_Internal;
  typedef typename KokkosKernels::Impl::GetUnifiedScalarViewType<
    BV, YMV, true>::type BV_Internal;
  typedef Kokkos::View<
    typename Kokkos::Impl::if_c<
      YMV::Rank == 1,
      typename YMV::non_const_value_type*,
      typename YMV::non_const_value_type** >::type,
    typename KokkosKernels::Impl::GetUnifiedLayout<YMV>::array_layout,
    typename YMV::device_type,
    Kokkos::MemoryTraits<Kokkos::Unmanaged> > YMV_Internal;

  AV_Internal  a_internal = a;
  XMV_Internal X_internal = X;
  BV_Internal  b_internal = b;
  YMV_Internal Y_internal = Y;

  Impl::Axpby<AV_Internal, XMV_Internal, BV_Internal,
    YMV_Internal>::axpby (a_internal, X_internal, b_internal, Y_internal);
}

template<class AV, class XMV, class YMV>
void
axpy (const AV& a, const XMV& X, const YMV& Y)
{
  axpby(a,X,Kokkos::Details::ArithTraits<typename YMV::non_const_value_type>::one(),Y);
}

} // KokkosBlas

#endif
