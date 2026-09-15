// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

#ifndef KOKKOSBLAS1_IAMAX_HPP_
#define KOKKOSBLAS1_IAMAX_HPP_

#include <KokkosBlas1_iamax_spec.hpp>
#include <KokkosKernels_helpers.hpp>
#include <KokkosKernels_Error.hpp>

namespace KokkosBlas {

/// \brief Return the (smallest) index of the element of the maximum magnitude
/// of the vector x.
///
/// \tparam execution_space a Kokkos execution space where the kernel will run.
/// \tparam XVector Type of the first vector x; a 1-D Kokkos::View.
///
/// \param space [in] execution space instance where the kernel will run.
/// \param x [in] Input 1-D View.
///
/// \return The (smallest) index of the element of the maximum magnitude; a
/// single value.
///         Note: Returned index is 1-based for compatibility with Fortran.
template <class execution_space, class XVector,
          typename std::enable_if<Kokkos::is_execution_space_v<execution_space>, int>::type = 0>
typename XVector::size_type iamax(const execution_space& space, const XVector& x) {
  static_assert(Kokkos::is_execution_space_v<execution_space>,
                "KokkosBlas::iamax: execution_space must be a valid Kokkos "
                "execution space");
  static_assert(Kokkos::is_view<XVector>::value, "KokkosBlas::iamax: XVector must be a Kokkos::View.");
  static_assert(Kokkos::SpaceAccessibility<execution_space, typename XVector::memory_space>::accessible,
                "KokkosBlas::iamax: XVector must be accessible from execution_space");
  static_assert(XVector::rank == 1,
                "KokkosBlas::iamax: "
                "Both Vector inputs must have rank 1.");

  using index_type = typename XVector::size_type;

  using XVector_Internal = Kokkos::View<typename XVector::const_value_type*,
                                        typename KokkosKernels::Impl::GetUnifiedLayout<XVector>::array_layout,
                                        execution_space, Kokkos::MemoryTraits<Kokkos::Unmanaged> >;

  using layout_t = typename XVector_Internal::array_layout;

  using RVector_Internal =
      Kokkos::View<index_type, layout_t, Kokkos::HostSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged> >;

  index_type result;
  RVector_Internal R = RVector_Internal(&result, layout_t());
  XVector_Internal X = x;

  Impl::Iamax<execution_space, RVector_Internal, XVector_Internal>::iamax(space, R, X);
  space.fence();
  return result;
}

/// \brief Return the (smallest) index of the element of the maximum magnitude
/// of the vector x.
///
/// The kernel is executed in the default stream/queue associated
/// with the execution space of XVector.
///
/// \tparam XVector Type of the first vector x; a 1-D Kokkos::View.
///
/// \param x [in] Input 1-D View.
///
/// \return The (smallest) index of the element of the maximum magnitude; a
/// single value.
///         Note: Returned index is 1-based for compatibility with Fortran.
template <class XVector>
typename XVector::size_type iamax(const XVector& x) {
  return iamax(typename XVector::execution_space{}, x);
}

/// \brief R(j) = iamax(X(i,j))
///
/// Replace each entry in R with the (smallest) index of the element of the
/// maximum magnitude of the corresponding entry in X.
/// This function is non-blocking and thread-safe.
///
/// \tparam RMV 0-D or 1-D Kokkos::View specialization.
/// \tparam XMV 1-D or 2-D Kokkos::View specialization.
///
/// \param space [in] execution space instance where the kernel will run.
/// \param R [out] Output View (rank 0 or 1) containing the results.
/// \param X [in] Input View (rank 1 or 2).
///
/// Note for TPL cuBLAS: When TPL cuBLAS iamax is used and returns result to a
/// view, RMV must be 0-D view and XMV must be 1-D view.
template <class execution_space, class RV, class XMV>
void iamax(const execution_space& space, const RV& R, const XMV& X,
           typename std::enable_if<Kokkos::is_view<RV>::value, int>::type = 0) {
  static_assert(Kokkos::is_execution_space_v<execution_space>,
                "KokkosBlas::iamax: execution_space must be a valid Kokkos "
                "execution space.");
  static_assert(Kokkos::is_view<RV>::value,
                "KokkosBlas::iamax: "
                "R is not a Kokkos::View.");
  static_assert(Kokkos::is_view<XMV>::value,
                "KokkosBlas::iamax: "
                "X is not a Kokkos::View.");
  static_assert(Kokkos::SpaceAccessibility<execution_space, typename XMV::memory_space>::accessible,
                "KokkosBlas::iamax: XMV must be accessible from execution_space.");
  static_assert(std::is_same<typename RV::value_type, typename RV::non_const_value_type>::value,
                "KokkosBlas::iamax: R is const.  "
                "It must be nonconst, because it is an output argument "
                "(we have to be able to write to its entries).");
  static_assert(((RV::rank == 0) && (XMV::rank == 1)) || ((RV::rank == 1) && (XMV::rank == 2)),
                "KokkosBlas::iamax: "
                "RV and XMV must either have rank 0 and 1 or rank 1 and 2.");

  using index_type = typename XMV::size_type;
  static_assert(std::is_same<typename RV::value_type, index_type>::value,
                "KokkosBlas::iamax: R must have the type of"
                "the Xvectors size_type it is an output argument "
                "(we have to be able to write to its entries).");

  // Check compatibility of dimensions at run time.
  if constexpr (XMV::rank == 2 && RV::rank == 1) {
    if (X.extent(1) != R.extent(0)) {
      std::ostringstream os;
      os << "KokkosBlas::iamax (MV): Dimensions of R and X do not match: "
         << "R: " << R.extent(0) << ", X: " << X.extent(0) << " x " << X.extent(1);
      KokkosKernels::Impl::throw_runtime_exception(os.str());
    }
  }

  using UnifiedXLayout  = typename KokkosKernels::Impl::GetUnifiedLayout<XMV>::array_layout;
  using UnifiedRVLayout = typename KokkosKernels::Impl::GetUnifiedLayoutPreferring<RV, UnifiedXLayout>::array_layout;
  using UnifiedRVDevice =
      std::conditional_t<Kokkos::SpaceAccessibility<Kokkos::HostSpace, typename RV::memory_space>::assignable,
                         Kokkos::HostSpace, execution_space>;

  // Create unmanaged versions of the input Views.  RV may be rank 0 or rank 1.
  // XMV may be rank 1 or rank 2.
  using RV_Internal  = Kokkos::View<typename RV::non_const_data_type, UnifiedRVLayout, UnifiedRVDevice,
                                   Kokkos::MemoryTraits<Kokkos::Unmanaged> >;
  using XMV_Internal = Kokkos::View<typename XMV::const_data_type, UnifiedXLayout, execution_space,
                                    Kokkos::MemoryTraits<Kokkos::Unmanaged> >;

  RV_Internal R_internal  = KokkosKernels::Impl::unificationCast<RV_Internal>(R);
  XMV_Internal X_internal = KokkosKernels::Impl::unificationCast<XMV_Internal>(X);

  Impl::Iamax<execution_space, RV_Internal, XMV_Internal>::iamax(space, R_internal, X_internal);
}

/// \brief R(j) = iamax(X(i,j))
///
/// Replace each entry in R with the (smallest) index of the element of the
/// maximum magnitude of the corresponding entry in X.
/// This function is non-blocking and thread-safe.
/// The kernel is executed in the default stream/queue associated
/// with the execution space of XVector.
///
/// \tparam RMV 0-D or 1-D Kokkos::View specialization.
/// \tparam XMV 1-D or 2-D Kokkos::View specialization.
///
/// Note for TPL cuBLAS: When TPL cuBLAS iamax is used and returns result to a
/// view, RMV must be 0-D view and XMV must be 1-D view.
template <class RV, class XMV>
void iamax(const RV& R, const XMV& X, typename std::enable_if<Kokkos::is_view<RV>::value, int>::type = 0) {
  iamax(typename XMV::execution_space{}, R, X);
}

}  // namespace KokkosBlas

#endif  // KOKKOSBLAS1_IAMAX_HPP_
