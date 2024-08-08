//@HEADER
// ************************************************************************
//
//                        Kokkos v. 4.0
//       Copyright (2022) National Technology & Engineering
//               Solutions of Sandia, LLC (NTESS).
//
// Under the terms of Contract DE-NA0003525 with NTESS,
// the U.S. Government retains certain rights in this software.
//
// Part of Kokkos, under the Apache License v2.0 with LLVM Exceptions.
// See https://kokkos.org/LICENSE for license information.
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
//@HEADER

#ifndef KOKKOSBLAS1_AXPBY_HPP_
#define KOKKOSBLAS1_AXPBY_HPP_

#if (KOKKOSKERNELS_DEBUG_LEVEL > 0)
#include <iostream>
#endif  // KOKKOSKERNELS_DEBUG_LEVEL

#include <KokkosBlas1_axpby_spec.hpp>
#include <KokkosBlas_serial_axpy.hpp>
#include <KokkosKernels_helpers.hpp>
#include <KokkosKernels_Error.hpp>
#include <KokkosBlas1_axpby_unification_attempt_traits.hpp>

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
  using AxpbyTraits   = Impl::AxpbyUnificationAttemptTraits<execution_space, AV, XMV, BV, YMV>;
  using InternalTypeA = typename AxpbyTraits::InternalTypeA;
  using InternalTypeX = typename AxpbyTraits::InternalTypeX;
  using InternalTypeB = typename AxpbyTraits::InternalTypeB;
  using InternalTypeY = typename AxpbyTraits::InternalTypeY;

  // **********************************************************************
  // Perform compile time checks and run time checks.
  // **********************************************************************
  AxpbyTraits::performChecks(a, X, b, Y);
#if (KOKKOSKERNELS_DEBUG_LEVEL > 1)
  AxpbyTraits::printInformation(std::cout, "axpby(), unif information");
#endif  // KOKKOSKERNELS_DEBUG_LEVEL

  // **********************************************************************
  // Call Impl::Axpby<...>::axpby(...)
  // **********************************************************************
  InternalTypeX internal_X = X;
  InternalTypeY internal_Y = Y;

  if constexpr (AxpbyTraits::internalTypesAB_bothScalars) {
    // ********************************************************************
    // The unification logic applies the following general rules:
    // 1) In a 'onHost' case, it makes the internal types for 'a' and 'b'
    //    to be both scalars (hence the name 'internalTypesAB_bothScalars')
    // 2) In a 'onDevice' case, it makes the internal types for 'a' and 'b'
    //    to be Kokkos views. For performance reasons in Trilinos, the only
    //    exception for this rule is when the input types for both 'a' and
    //    'b' are already scalars, in which case the internal types for 'a'
    //    and 'b' become scalars as well, eventually changing precision in
    //    order to match the precisions of 'X' and 'Y'.
    // ********************************************************************
    if constexpr (AxpbyTraits::a_is_scalar && AxpbyTraits::b_is_scalar && AxpbyTraits::onDevice) {
      // ******************************************************************
      // We are in the exception situation for rule 2
      // ******************************************************************
      InternalTypeA internal_a(a);
      InternalTypeA internal_b(b);

      Impl::Axpby<execution_space, InternalTypeA, InternalTypeX, InternalTypeB, InternalTypeY>::axpby(
          exec_space, internal_a, internal_X, internal_b, internal_Y);
    } else {
      // ******************************************************************
      // We are in rule 1, that is, we are in a 'onHost' case now
      // ******************************************************************
      InternalTypeA internal_a(Impl::getScalarValueFromVariableAtHost<AV, Impl::typeRank<AV>()>::getValue(a));
      InternalTypeB internal_b(Impl::getScalarValueFromVariableAtHost<BV, Impl::typeRank<BV>()>::getValue(b));

      Impl::Axpby<execution_space, InternalTypeA, InternalTypeX, InternalTypeB, InternalTypeY>::axpby(
          exec_space, internal_a, internal_X, internal_b, internal_Y);
    }
  } else if constexpr (AxpbyTraits::internalTypesAB_bothViews) {
    constexpr bool internalLayoutA_isStride(std::is_same_v<typename InternalTypeA::array_layout, Kokkos::LayoutStride>);
    constexpr bool internalLayoutB_isStride(std::is_same_v<typename InternalTypeB::array_layout, Kokkos::LayoutStride>);

    const size_t numScalarsA(Impl::getAmountOfScalarsInCoefficient(a));
    const size_t numScalarsB(Impl::getAmountOfScalarsInCoefficient(b));

    const size_t strideA(Impl::getStrideInCoefficient(a));
    const size_t strideB(Impl::getStrideInCoefficient(b));

    Kokkos::LayoutStride layoutStrideA{numScalarsA, strideA};
    Kokkos::LayoutStride layoutStrideB{numScalarsB, strideB};

    InternalTypeA internal_a;
    InternalTypeB internal_b;

    if constexpr (internalLayoutA_isStride) {
      // ******************************************************************
      // Prepare internal_a
      // ******************************************************************
      typename AxpbyTraits::InternalTypeA_managed managed_a("managed_a", layoutStrideA);
      if constexpr (AxpbyTraits::atInputLayoutA_isStride) {
        Kokkos::deep_copy(managed_a, a);
      } else {
        Impl::populateRank1Stride1ViewWithScalarOrNonStrideView(a, managed_a);
      }
      internal_a = managed_a;

      if constexpr (internalLayoutB_isStride) {
        // ****************************************************************
        // Prepare internal_b
        // ****************************************************************
        typename AxpbyTraits::InternalTypeB_managed managed_b("managed_b", layoutStrideB);
        if constexpr (AxpbyTraits::atInputLayoutB_isStride) {
          Kokkos::deep_copy(managed_b, b);
        } else {
          Impl::populateRank1Stride1ViewWithScalarOrNonStrideView(b, managed_b);
        }
        internal_b = managed_b;

        // ****************************************************************
        // Call Impl::Axpby<...>::axpby(...)
        // ****************************************************************
        Impl::Axpby<execution_space, InternalTypeA, InternalTypeX, InternalTypeB, InternalTypeY>::axpby(
            exec_space, internal_a, internal_X, internal_b, internal_Y);
      } else {
        // ****************************************************************
        // Prepare internal_b
        // ****************************************************************
        typename AxpbyTraits::InternalTypeB_managed managed_b("managed_b", numScalarsB);
        if constexpr (AxpbyTraits::atInputLayoutB_isStride) {
          Kokkos::deep_copy(managed_b, b);
        } else {
          Impl::populateRank1Stride1ViewWithScalarOrNonStrideView(b, managed_b);
        }
        internal_b = managed_b;

        // ****************************************************************
        // Call Impl::Axpby<...>::axpby(...)
        // ****************************************************************
        Impl::Axpby<execution_space, InternalTypeA, InternalTypeX, InternalTypeB, InternalTypeY>::axpby(
            exec_space, internal_a, internal_X, internal_b, internal_Y);
      }
    } else {
      // ******************************************************************
      // Prepare internal_a
      // ******************************************************************
      typename AxpbyTraits::InternalTypeA_managed managed_a("managed_a", numScalarsA);
      if constexpr (AxpbyTraits::atInputLayoutA_isStride) {
        Kokkos::deep_copy(managed_a, a);
      } else {
        Impl::populateRank1Stride1ViewWithScalarOrNonStrideView(a, managed_a);
      }
      internal_a = managed_a;

      if constexpr (internalLayoutB_isStride) {
        // ****************************************************************
        // Prepare internal_b
        // ****************************************************************
        typename AxpbyTraits::InternalTypeB_managed managed_b("managed_b", layoutStrideB);
        if constexpr (AxpbyTraits::atInputLayoutB_isStride) {
          Kokkos::deep_copy(managed_b, b);
        } else {
          Impl::populateRank1Stride1ViewWithScalarOrNonStrideView(b, managed_b);
        }
        internal_b = managed_b;

        // ****************************************************************
        // Call Impl::Axpby<...>::axpby(...)
        // ****************************************************************
        Impl::Axpby<execution_space, InternalTypeA, InternalTypeX, InternalTypeB, InternalTypeY>::axpby(
            exec_space, internal_a, internal_X, internal_b, internal_Y);
      } else {
        // ****************************************************************
        // Prepare internal_b
        // ****************************************************************
        typename AxpbyTraits::InternalTypeB_managed managed_b("managed_b", numScalarsB);
        if constexpr (AxpbyTraits::atInputLayoutB_isStride) {
          Kokkos::deep_copy(managed_b, b);
        } else {
          Impl::populateRank1Stride1ViewWithScalarOrNonStrideView(b, managed_b);
        }
        internal_b = managed_b;

        // ****************************************************************
        // Call Impl::Axpby<...>::axpby(...)
        // ****************************************************************
        Impl::Axpby<execution_space, InternalTypeA, InternalTypeX, InternalTypeB, InternalTypeY>::axpby(
            exec_space, internal_a, internal_X, internal_b, internal_Y);
      }
    }
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
  axpby(exec_space, a, X, Kokkos::ArithTraits<typename YMV::non_const_value_type>::one(), Y);
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
#if (KOKKOSKERNELS_DEBUG_LEVEL > 0)
  static_assert(Kokkos::is_view<XMV>::value, "KokkosBlas::serial_axpy: XMV is not a Kokkos::View");
  static_assert(Kokkos::is_view<YMV>::value, "KokkosBlas::serial_axpy: YMV is not a Kokkos::View");
  static_assert(XMV::rank == 1 || XMV::rank == 2, "KokkosBlas::serial_axpy: XMV must have rank 1 or 2.");
  static_assert(XMV::rank == YMV::rank, "KokkosBlas::serial_axpy: XMV and YMV must have the same rank.");

  if (X.extent(0) != Y.extent(0) || X.extent(1) != Y.extent(1)) {
    Kokkos::abort("KokkosBlas::serial_axpy: X and Y dimensions do not match");
  }
#endif  // KOKKOSKERNELS_DEBUG_LEVEL

  return Impl::serial_axpy_mv(X.extent(0), X.extent(1), alpha, X.data(), Y.data(), X.stride_0(), X.stride_1(),
                              Y.stride_0(), Y.stride_1());
}

}  // namespace KokkosBlas

#endif
