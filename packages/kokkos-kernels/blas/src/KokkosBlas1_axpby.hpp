// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

#ifndef KOKKOSBLAS1_AXPBY_HPP
#define KOKKOSBLAS1_AXPBY_HPP

#include <KokkosBlas1_axpby_spec.hpp>
#include <KokkosBlas_serial_axpy.hpp>
#include <KokkosKernels_helpers.hpp>
#include <KokkosBlas1_axpby_unification.hpp>
#include <KokkosKernels_Error.hpp>

// axpby() accepts both scalar coefficients a and b, and vector
// coefficients (apply one for each column of the input multivectors).
// This traits class helps axpby() select the correct specialization
// of AV (type of 'a') and BV (type of 'b') for invoking the
// implementation.

namespace KokkosBlas {

/// \brief Computes Y := a*X + b*Y
///
/// This function is non-blocking and thread-safe.
///
/// \tparam execution_space The type of execution space where the kernel
///                         will run.
/// \tparam AV              Scalar or 0-D or 1-D Kokkos::View.
/// \tparam XMV             1-D Kokkos::View or 2-D Kokkos::View. It
///                         must have the same rank as YMV.
/// \tparam BV              Scalar or 0-D or 1-D Kokkos::View.
/// \tparam YMV             1-D or 2-D Kokkos::View.
///
/// \param exec_space [in] The execution space instance on which the kernel
///                        will run.
/// \param a          [in] Input of type AV:
///                        - scaling parameter for 1-D or 2-D X,
///                        - scaling parameters for 2-D X.
/// \param X          [in] View of type XMV. It must have the same
///                        extent(s) as Y.
/// \param b          [in] input of type BV:
///                        - scaling parameter for 1-D or 2-D Y,
///                        - scaling parameters for 2-D Y.
/// \param Y      [in/out] View of type YMV in which the results will be
///                        stored.
template <class execution_space, class AV, class XMV, class BV, class YMV>
void axpby(const execution_space& exec_space, const AV& a, const XMV& X, const BV& b, const YMV& Y) {
  static_assert(Kokkos::is_execution_space_v<execution_space>,
                "KokkosBlas::axpby: execution_space must be a valid Kokkos "
                "execution space");
  static_assert(Kokkos::is_view_v<XMV>, "KokkosBlas::axpby: X must be a Kokkos::View.");
  static_assert(Kokkos::SpaceAccessibility<execution_space, typename XMV::memory_space>::accessible,
                "KokkosBlas::axpby: X must be accessible from execution_space.");
  static_assert(Kokkos::is_view_v<YMV>, "KokkosBlas::axpby: Y must be a Kokkos::View.");
  static_assert(Kokkos::SpaceAccessibility<execution_space, typename YMV::memory_space>::accessible,
                "KokkosBlas::axpby: Y must be accessible from execution_space.");
  static_assert(std::is_same<typename YMV::value_type, typename YMV::non_const_value_type>::value,
                "KokkosBlas::axpby: Y is const.  "
                "It must be nonconst, because it is an output argument "
                "(we have to be able to write to its entries).");
  static_assert((int)XMV::rank == (int)YMV::rank,
                "KokkosBlas::axpby: "
                "X and Y must have the same rank");
  static_assert((int)XMV::rank == 1 || (int)XMV::rank == 2,
                "KokkosBlas::axpby: X and Y must be either rank 1 or rank 2");
  // Check compatibility of dimensions at run time. We already checked X and Y have matching ranks.
  if constexpr ((int)XMV::rank == 1) {
    if (X.extent(0) != Y.extent(0)) {
      std::ostringstream os;
      os << "KokkosBlas::axpby: Dimensions of X and Y do not match: "
         << "X: " << X.extent(0) << ", Y: " << Y.extent(0);
      KokkosKernels::Impl::throw_runtime_exception(os.str());
    }
  } else if constexpr ((int)XMV::rank == 2) {
    if (X.extent(0) != Y.extent(0) || X.extent(1) != Y.extent(1)) {
      std::ostringstream os;
      os << "KokkosBlas::axpby: Dimensions of X and Y do not match: "
         << "X: " << X.extent(0) << " x " << X.extent(1) << ", Y: " << Y.extent(0) << " x " << Y.extent(1);
      KokkosKernels::Impl::throw_runtime_exception(os.str());
    }
  }
  // Check coefficients (both types and runtime extents).
  // They can be either rank-0 in host or device space, OR rank-1 in device space.
  if constexpr (Kokkos::is_view_v<AV>) {
    if constexpr ((int)AV::rank == 1) {
      static_assert(Kokkos::SpaceAccessibility<execution_space, typename AV::memory_space>::accessible,
                    "KokkosBlas::axpby: when a is a rank-1 View, it must be accessible from execution_space.");
      // And check extent at runtime
      if (a.extent(0) != size_t(1) && a.extent(0) != X.extent(1)) {
        std::ostringstream os;
        os << "KokkosBlas::axpby: Dimensions of a and X do not match: "
           << "a has extent " << a.extent(0) << " but X has " << X.extent(1) << " columns";
        KokkosKernels::Impl::throw_runtime_exception(os.str());
      }
    }
  }
  if constexpr (Kokkos::is_view_v<BV>) {
    if constexpr ((int)BV::rank == 1) {
      static_assert(Kokkos::SpaceAccessibility<execution_space, typename BV::memory_space>::accessible,
                    "KokkosBlas::axpby: when b is a rank-1 View, it must be accessible from execution_space.");
      // And check extent at runtime
      if (b.extent(0) != size_t(1) && b.extent(0) != Y.extent(1)) {
        std::ostringstream os;
        os << "KokkosBlas::axpby: Dimensions of b and Y do not match: "
           << "b has extent " << b.extent(0) << " but Y has " << Y.extent(1) << " columns";
        KokkosKernels::Impl::throw_runtime_exception(os.str());
      }
    }
  }
  // Get unified types (all preferring the layout of the InternalTypeX)
  using PreferredLayout = typename KokkosKernels::Impl::GetUnifiedLayout<XMV>::array_layout;
  using InternalTypeX   = Kokkos::View<typename XMV::const_data_type, PreferredLayout, typename XMV::device_type,
                                     Kokkos::MemoryTraits<Kokkos::Unmanaged> >;
  using InternalTypeY =
      Kokkos::View<typename YMV::non_const_data_type,
                   typename KokkosKernels::Impl::GetUnifiedLayoutPreferring<YMV, PreferredLayout>::array_layout,
                   typename YMV::device_type, Kokkos::MemoryTraits<Kokkos::Unmanaged> >;
  InternalTypeX internal_X = X;
  InternalTypeY internal_Y = Y;
  // Run the rank-1 or rank-2 cases, with the appropriate coefficient unification
  if constexpr ((int)XMV::rank == 1) {
    using InternalTypeA =
        typename KokkosBlas::Impl::UnifiedAxpbyCoeff<AV, typename XMV::non_const_value_type, PreferredLayout>::type;
    using InternalTypeB =
        typename KokkosBlas::Impl::UnifiedAxpbyCoeff<BV, typename YMV::non_const_value_type, PreferredLayout>::type;
    InternalTypeA internal_a = KokkosBlas::Impl::unifyAxpbyCoeff<InternalTypeA>(a);
    InternalTypeB internal_b = KokkosBlas::Impl::unifyAxpbyCoeff<InternalTypeB>(b);
    Impl::Axpby<execution_space, InternalTypeA, InternalTypeX, InternalTypeB, InternalTypeY>::axpby(
        exec_space, internal_a, internal_X, internal_b, internal_Y);
  } else {
    // Note BMK: axpby_mv implementation has special support for rank-1, extent 1 coefficients (a and/or b)
    // even when X and Y have >1 column. This case is treated as if the coefficient were rank-0.
    // If TPL support is added for axpby_mv, we need to handle this case at runtime because TPLs likely won't support
    // it.
    using InternalTypeA =
        typename KokkosBlas::Impl::UnifiedAxpbyMvCoeff<AV, typename XMV::non_const_value_type, PreferredLayout>::type;
    using InternalTypeB =
        typename KokkosBlas::Impl::UnifiedAxpbyMvCoeff<BV, typename YMV::non_const_value_type, PreferredLayout>::type;
    InternalTypeA internal_a = KokkosBlas::Impl::unifyAxpbyMvCoeff<InternalTypeA>(a);
    InternalTypeB internal_b = KokkosBlas::Impl::unifyAxpbyMvCoeff<InternalTypeB>(b);
    Impl::Axpby<execution_space, InternalTypeA, InternalTypeX, InternalTypeB, InternalTypeY>::axpby(
        exec_space, internal_a, internal_X, internal_b, internal_Y);
  }
}

/// \brief Computes Y := a*X + b*Y
///
/// This function is non-blocking and thread-safe.
/// The kernel is executed in the default stream/queue
/// associated with the execution space of XMV.
///
/// \tparam AV  Scalar or 0-D Kokkos::View or 1-D Kokkos::View.
/// \tparam XMV 1-D Kokkos::View or 2-D Kokkos::View. It must
///             have the same rank as YMV.
/// \tparam BV  Scalar or 0-D Kokkos::View or 1-D Kokkos::View.
/// \tparam YMV 1-D Kokkos::View or 2-D Kokkos::View.
///
/// \param a     [in] Input of type AV:
///                   - scaling parameter for 1-D or 2-D X,
///                   - scaling parameters for 2-D X.
/// \param X     [in] View of type XMV. It must have the same
///                   extent(s) as Y.
/// \param b     [in] input of type BV:
///                   - scaling parameter for 1-D or 2-D Y,
///                   - scaling parameters for 2-D Y.
/// \param Y [in/out] View of type YMV in which the results will be
///                   stored.
template <class AV, class XMV, class BV, class YMV>
void axpby(const AV& a, const XMV& X, const BV& b, const YMV& Y) {
  axpby(typename XMV::execution_space{}, a, X, b, Y);
}

/// \brief Computes Y := a*X + Y
///
/// This function is non-blocking and thread-safe.
///
/// \tparam execution_space The type of execution space where the kernel
///                         will run.
/// \tparam AV              Scalar or 0-D or 1-D Kokkos::View.
/// \tparam XMV             1-D or 2-D Kokkos::View. It must have the
///                         the same rank as YMV.
/// \tparam YMV             1-D or 2-D Kokkos::View.
///
/// \param exec_space [in] The execution space instance on which the kernel
///                        will run.
/// \param a          [in] Input of type AV:
///                        - scaling parameter for 1-D or 2-D X,
///                        - scaling parameters for 2-D X.
/// \param X          [in] View of type XMV. It must have the same
///                        extent(s) as Y.
/// \param Y      [in/out] View of type YMV in which the results will be
///                        stored.
template <class execution_space, class AV, class XMV, class YMV>
void axpy(const execution_space& exec_space, const AV& a, const XMV& X, const YMV& Y) {
  axpby(exec_space, a, X, KokkosKernels::ArithTraits<typename YMV::non_const_value_type>::one(), Y);
}

/// \brief Computes Y := a*X + Y
///
/// This function is non-blocking and thread-safe.
/// The kernel is executed in the default stream/queue
/// associated with the execution space of XMV.
///
/// \tparam AV  Scalar or 0-D Kokkos::View or 1-D Kokkos::View.
/// \tparam XMV 1-D Kokkos::View or 2-D Kokkos::View. It must
///             have the same rank as YMV.
/// \tparam YMV 1-D Kokkos::View or 2-D Kokkos::View.
///
/// \param a     [in] Input of type AV:
///                   - scaling parameter for 1-D or 2-D X,
///                   - scaling parameters for 2-D X.
/// \param X     [in] View of type XMV. It must have the same
///                   extent(s) as Y.
/// \param Y [in/out] View of type YMV in which the results will be
///                   stored.
template <class AV, class XMV, class YMV>
void axpy(const AV& a, const XMV& X, const YMV& Y) {
  axpy(typename XMV::execution_space{}, a, X, Y);
}

///
/// Serial axpy on device
///
template <class scalar_type, class XMV, class YMV>
KOKKOS_FUNCTION void serial_axpy(const scalar_type alpha, const XMV X, YMV Y) {
  static_assert(Kokkos::is_view<XMV>::value, "KokkosBlas::serial_axpy: XMV is not a Kokkos::View");
  static_assert(Kokkos::is_view<YMV>::value, "KokkosBlas::serial_axpy: YMV is not a Kokkos::View");
  static_assert(XMV::rank == 1 || XMV::rank == 2, "KokkosBlas::serial_axpy: XMV must have rank 1 or 2.");
  static_assert(XMV::rank == YMV::rank, "KokkosBlas::serial_axpy: XMV and YMV must have the same rank.");

#ifndef NDEBUG
  if (X.extent(0) != Y.extent(0) || X.extent(1) != Y.extent(1)) {
    Kokkos::abort("KokkosBlas::serial_axpy: X and Y dimensions do not match");
  }
#endif  // NDEBUG
  if constexpr (XMV::rank() == 1)
    return Impl::serial_axpy(X.extent(0), alpha, X.data(), Y.data(), X.stride(0), Y.stride(0));
  else
    return Impl::serial_axpy_mv(X.extent(0), X.extent(1), alpha, X.data(), Y.data(), X.stride(0), X.stride(1),
                                Y.stride(0), Y.stride(1));
}

}  // namespace KokkosBlas

#endif  // KOKKOSBLAS1_AXPBY_HPP
