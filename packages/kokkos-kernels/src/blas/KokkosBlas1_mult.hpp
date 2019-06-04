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

#ifndef KOKKOSBLAS1_MULT_HPP_
#define KOKKOSBLAS1_MULT_HPP_

#include<KokkosBlas1_mult_spec.hpp>
#include<KokkosKernels_helpers.hpp>

namespace KokkosBlas {

template<class YMV, class AV, class XMV>
void
mult (typename YMV::const_value_type& gamma,
      const YMV& Y,
      typename AV::const_value_type& alpha,
      const AV& A,
      const XMV& X)
{
  static_assert (Kokkos::Impl::is_view<YMV>::value, "KokkosBlas::mult: "
                 "Y is not a Kokkos::View.");
  static_assert (Kokkos::Impl::is_view<AV>::value, "KokkosBlas::mult: "
                 "A is not a Kokkos::View.");
  static_assert (Kokkos::Impl::is_same<typename YMV::value_type,
                   typename YMV::non_const_value_type>::value,
                 "KokkosBlas::mult: Y is const.  "
                 "It must be nonconst, because it is an output argument "
                 "(we have to be able to write to its entries).");
  static_assert ((XMV::rank == 1 && YMV::rank == 1) ||
                 (XMV::rank == 2 && YMV::rank == 2),
                 "KokkosBlas::mult: Y and X must be either both rank 1, "
                 "or both rank 2.");
  static_assert (AV::rank == 1, "KokkosBlas::mult: A must have rank 1.");

  // Check compatibility of dimensions at run time.
  if (Y.extent(0) != A.extent(0) ||
      Y.extent(0) != X.extent(0) ||
      Y.extent(1) != X.extent(1)) {
    std::ostringstream os;
    os << "KokkosBlas::mult: Dimensions do not match: "
       << "Y: " << Y.extent(0) << " x " << Y.extent(1)
       << ", A: " << A.extent(0) << " x " << A.extent(0)
       << ", X: " << X.extent(0) << " x " << X.extent(1);
    Kokkos::Impl::throw_runtime_exception (os.str ());
  }

  // Create unmanaged versions of the input Views.
  typedef Kokkos::View<
    typename Kokkos::Impl::if_c<
      YMV::rank == 1,
      typename YMV::non_const_value_type*,
      typename YMV::non_const_value_type** >::type,
    typename KokkosKernels::Impl::GetUnifiedLayout<YMV>::array_layout,
    typename YMV::device_type,
    Kokkos::MemoryTraits<Kokkos::Unmanaged> > YMV_Internal;
  typedef Kokkos::View<
    typename AV::const_value_type*,
    typename KokkosKernels::Impl::GetUnifiedLayout<AV>::array_layout,
    typename AV::device_type,
    Kokkos::MemoryTraits<Kokkos::Unmanaged> > AV_Internal;
  typedef Kokkos::View<
    typename Kokkos::Impl::if_c<
      XMV::rank == 1,
      typename XMV::const_value_type*,
      typename XMV::const_value_type** >::type,
    typename KokkosKernels::Impl::GetUnifiedLayout<XMV>::array_layout,
    typename XMV::device_type,
    Kokkos::MemoryTraits<Kokkos::Unmanaged> > XMV_Internal;

  YMV_Internal Y_internal = Y;
  AV_Internal A_internal = A;
  XMV_Internal X_internal = X;

  Impl::Mult<YMV_Internal, AV_Internal, XMV_Internal>::mult (gamma, Y_internal, alpha,
                                                             A_internal, X_internal);
}


}

#endif // KOKKOSBLAS1_MULT_HPP_

