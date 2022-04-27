/*
//@HEADER
// ************************************************************************
//
//                        Kokkos v. 3.0
//       Copyright (2020) National Technology & Engineering
//               Solutions of Sandia, LLC (NTESS).
//
// Under the terms of Contract DE-NA0003525 with NTESS,
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
// THIS SOFTWARE IS PROVIDED BY NTESS "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL NTESS OR THE
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

#ifndef KOKKOSBLAS1_SCAL_HPP_
#define KOKKOSBLAS1_SCAL_HPP_

#include <KokkosBlas1_scal_spec.hpp>
#include <KokkosKernels_helpers.hpp>
#include <KokkosKernels_Error.hpp>

namespace KokkosBlas {

template <class RMV, class AV, class XMV>
void scal(const RMV& R, const AV& a, const XMV& X) {
  static_assert(Kokkos::is_view<RMV>::value,
                "KokkosBlas::scal: "
                "R is not a Kokkos::View.");
  static_assert(Kokkos::is_view<XMV>::value,
                "KokkosBlas::scal: "
                "X is not a Kokkos::View.");
  static_assert(std::is_same<typename RMV::value_type,
                             typename RMV::non_const_value_type>::value,
                "KokkosBlas::scal: R is const.  "
                "It must be nonconst, because it is an output argument "
                "(we have to be able to write to its entries).");
  static_assert((int)RMV::rank == (int)XMV::rank,
                "KokkosBlas::scal: "
                "R and X must have the same rank.");
  static_assert(RMV::rank == 1 || RMV::rank == 2,
                "KokkosBlas::scal: "
                "RMV and XMV must either have rank 1 or rank 2.");

  // Check compatibility of dimensions at run time.
  if (X.extent(0) != R.extent(0) || X.extent(1) != R.extent(1)) {
    std::ostringstream os;
    os << "KokkosBlas::scal: Dimensions of R and X do not match: "
       << "R: " << R.extent(0) << " x " << R.extent(1) << ", X: " << X.extent(0)
       << " x " << X.extent(1);
    KokkosKernels::Impl::throw_runtime_exception(os.str());
  }

  using UnifiedRLayout =
      typename KokkosKernels::Impl::GetUnifiedLayout<RMV>::array_layout;
  using UnifiedXLayout =
      typename KokkosKernels::Impl::GetUnifiedLayoutPreferring<
          XMV, UnifiedRLayout>::array_layout;

  // Create unmanaged versions of the input Views.  RMV and XMV may be
  // rank 1 or rank 2.  AV may be either a rank-1 View, or a scalar
  // value.
  typedef Kokkos::View<typename RMV::non_const_data_type, UnifiedRLayout,
                       typename RMV::device_type,
                       Kokkos::MemoryTraits<Kokkos::Unmanaged> >
      RMV_Internal;
  typedef Kokkos::View<typename XMV::const_data_type, UnifiedXLayout,
                       typename XMV::device_type,
                       Kokkos::MemoryTraits<Kokkos::Unmanaged> >
      XMV_Internal;
  typedef typename KokkosKernels::Impl::GetUnifiedScalarViewType<
      AV, XMV_Internal, true>::type AV_Internal;

  RMV_Internal R_internal = R;
  AV_Internal a_internal  = a;
  XMV_Internal X_internal = X;

  Impl::Scal<RMV_Internal, AV_Internal, XMV_Internal>::scal(
      R_internal, a_internal, X_internal);
}

}  // namespace KokkosBlas

#endif
