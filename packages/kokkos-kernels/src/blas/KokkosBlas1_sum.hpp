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

#ifndef KOKKOSBLAS1_SUM_HPP_
#define KOKKOSBLAS1_SUM_HPP_

#include<KokkosBlas1_sum_spec.hpp>
#include<KokkosKernels_helpers.hpp>

namespace KokkosBlas {

/// \brief Return the sum of the vector x.
///
/// \tparam XVector Type of the first vector x; a 1-D Kokkos::View.
///
/// \param x [in] Input 1-D View.
///
/// \return The sum product result; a single value.
template<class XVector>
typename XVector::non_const_value_type
sum (const XVector& x)
{
  static_assert (Kokkos::Impl::is_view<XVector>::value,
                 "KokkosBlas::sum: XVector must be a Kokkos::View.");
  static_assert (XVector::rank == 1, "KokkosBlas::sum: "
                 "Both Vector inputs must have rank 1.");

  typedef Kokkos::View<typename XVector::const_value_type*,
    typename KokkosKernels::Impl::GetUnifiedLayout<XVector>::array_layout,
    typename XVector::device_type,
    Kokkos::MemoryTraits<Kokkos::Unmanaged> > XVector_Internal;

  typedef Kokkos::View<
    typename XVector::non_const_value_type,
    Kokkos::LayoutLeft,
    Kokkos::HostSpace,
    Kokkos::MemoryTraits<Kokkos::Unmanaged> > RVector_Internal;

  typename XVector::non_const_value_type result;
  RVector_Internal R = RVector_Internal(&result);
  XVector_Internal X = x;

  Impl::Sum<RVector_Internal,XVector_Internal>::sum (R,X);
  Kokkos::fence();
  return result;
}

/// \brief R(j) = sum(X(i,j))
///
/// Replace each entry in R with the sumolute value (magnitude) of the
/// corresponding entry in X.
///
/// \tparam RMV 1-D or 2-D Kokkos::View specialization.
/// \tparam XMV 1-D or 2-D Kokkos::View specialization.  It must have
///   the same rank as RMV, and its entries must be assignable to
///   those of RMV.
template<class RV, class XMV>
void
sum (const RV& R, const XMV& X,
      typename std::enable_if<Kokkos::Impl::is_view<RV>::value, int>::type = 0)
{
  static_assert (Kokkos::Impl::is_view<RV>::value, "KokkosBlas::sum: "
                 "R is not a Kokkos::View.");
  static_assert (Kokkos::Impl::is_view<XMV>::value, "KokkosBlas::sum: "
                 "X is not a Kokkos::View.");
  static_assert (Kokkos::Impl::is_same<typename RV::value_type,
                 typename RV::non_const_value_type>::value,
                 "KokkosBlas::sum: R is const.  "
                 "It must be nonconst, because it is an output argument "
                 "(we have to be able to write to its entries).");
  static_assert (((RV::rank == 0) && (XMV::rank == 1)) ||
                 ((RV::rank == 1) && (XMV::rank == 2)), "KokkosBlas::sum: "
                 "RV and XMV must either have rank 0 and 1 or rank 1 and 2.");

  // Check compatibility of dimensions at run time.
  if (X.extent(1) != R.extent(0)) {
    std::ostringstream os;
    os << "KokkosBlas::sum (MV): Dimensions of R and X do not match: "
       << "R: " << R.extent(0)
       << ", X: " << X.extent(0) << " x " << X.extent(1);
    Kokkos::Impl::throw_runtime_exception (os.str ());
  }

  // Create unmanaged versions of the input Views.  RV and XMV may be
  // rank 1 or rank 2.
  typedef Kokkos::View<
    typename Kokkos::Impl::if_c<
      RV::rank == 0,
      typename RV::non_const_value_type,
      typename RV::non_const_value_type* >::type,
    typename KokkosKernels::Impl::GetUnifiedLayout<RV>::array_layout,
    typename RV::device_type,
    Kokkos::MemoryTraits<Kokkos::Unmanaged> > RV_Internal;
  typedef Kokkos::View<
    typename Kokkos::Impl::if_c<
      XMV::rank == 1,
      typename XMV::const_value_type*,
      typename XMV::const_value_type** >::type,
    typename KokkosKernels::Impl::GetUnifiedLayout<XMV>::array_layout,
    typename XMV::device_type,
    Kokkos::MemoryTraits<Kokkos::Unmanaged> > XMV_Internal;

  RV_Internal R_internal = R;
  XMV_Internal X_internal = X;

  Impl::Sum<RV_Internal, XMV_Internal>::sum (R_internal, X_internal);
}

}

#endif // KOKKOSBLAS1_SUM_HPP_

