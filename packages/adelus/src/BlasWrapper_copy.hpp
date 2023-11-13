/*
//@HEADER
// ************************************************************************
//
//                        Adelus v. 1.0
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
// 3. Neither the name of NTESS nor the names of the contributors may be
// used to endorse or promote products derived from this software without
// specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY NTESS "AS IS" AND ANY EXPRESS OR IMPLIED
// WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF
// MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED.
// IN NO EVENT SHALL NTESS OR THE CONTRIBUTORS BE LIABLE FOR ANY DIRECT,
// INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
// (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR 
// SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
// HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT,
// STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING
// IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE 
// POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Vinh Dang (vqdang@sandia.gov)
//                    Joseph Kotulski (jdkotul@sandia.gov)
//                    Siva Rajamanickam (srajama@sandia.gov)
//
// ************************************************************************
//@HEADER
*/

#ifndef BLASWRAPPER_COPY_HPP_
#define BLASWRAPPER_COPY_HPP_

#include <stdexcept>

#include<BlasWrapper_copy_spec.hpp>
#include<KokkosKernels_helpers.hpp>

namespace BlasWrapper {

/// \brief Y(i,j) = X(i,j)
///
/// Copy each entry in X into the corresponding entry in Y.
///
/// \tparam YMV 1-D or 2-D Kokkos::View specialization.
/// \tparam XMV 1-D or 2-D Kokkos::View specialization.  It must have
///   the same rank as YMV, and its entries must be assignable to
///   those of YMV.
template<class XMV, class YMV>
void
copy (const XMV& X, const YMV& Y)
{
  static_assert (Kokkos::is_view<XMV>::value, "BlasWrapper::copy: "
                 "X is not a Kokkos::View.");
  static_assert (Kokkos::is_view<YMV>::value, "BlasWrapper::copy: "
                 "Y is not a Kokkos::View.");
  static_assert (std::is_same<typename YMV::value_type,
                 typename YMV::non_const_value_type>::value,
                 "BlasWrapper::copy: Y is const.  "
                 "It must be nonconst, because it is an output argument "
                 "(we have to be able to write to its entries).");
  static_assert (int(XMV::rank) == int(YMV::rank), "BlasWrapper::copy: "
                 "X and Y must have the same rank.");
  static_assert (YMV::rank == 1 || YMV::rank == 2, "BlasWrapper::copy: "
                 "YMV and XMV must either have rank 1 or rank 2.");

  // Check compatibility of dimensions at run time.
  if (X.extent(0) != Y.extent(0) ||
      X.extent(1) != Y.extent(1)) {
    std::ostringstream os;
    os << "BlasWrapper::copy (MV): Dimensions of Y and X do not match: "
       << "Y: " << Y.extent(0) << " x " << Y.extent(1)
       << ", X: " << X.extent(0) << " x " << X.extent(1);
    throw std::runtime_error (os.str ());
  }

  // Create unmanaged versions of the input Views.  RMV and XMV may be
  // rank 1 or rank 2.
  typedef Kokkos::View<
    typename std::conditional<
      YMV::rank == 1,
      typename YMV::non_const_value_type*,
      typename YMV::non_const_value_type** >::type,
    typename KokkosKernels::Impl::GetUnifiedLayout<YMV>::array_layout,
    typename YMV::device_type,
    Kokkos::MemoryTraits<Kokkos::Unmanaged> > YMV_Internal;
  typedef Kokkos::View<
    typename std::conditional<
      XMV::rank == 1,
      typename XMV::const_value_type*,
      typename XMV::const_value_type** >::type,
    typename KokkosKernels::Impl::GetUnifiedLayout<XMV>::array_layout,
    typename XMV::device_type,
    Kokkos::MemoryTraits<Kokkos::Unmanaged> > XMV_Internal;

  YMV_Internal Y_internal = Y;
  XMV_Internal X_internal = X;

  Impl::Copy<XMV_Internal, YMV_Internal>::copy (X_internal, Y_internal);
}
}

#endif // BLASWRAPPER_COPY_HPP_
